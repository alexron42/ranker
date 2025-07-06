{
  library(dplyr)
  library(readr)
  library(purrr)
  library(beepr)
  library(parallel)
  library(RobustRankAggreg)
  library(Matrix)
  library(tibble)
  library(ggplot2)
  library(tidyr)
}

compare_ranks <- function(
    csv_file,
    sample_sizes = 1000,
    cores = parallel::detectCores() - 2
) {
  cat("Starting compare_ranks on ", csv_file, "\n"); flush.console()
  cl <- makeCluster(cores)
  on.exit({ cat("Stopping cluster\n"); flush.console(); stopCluster(cl) }, add = TRUE)
  clusterEvalQ(cl, { library(Matrix); library(RobustRankAggreg) })
  
  get_steady_state <- function(P, tol = 1e-10, max_iter = 1000) {
    v <- rep(1 / nrow(P), nrow(P))
    for (i in seq_len(max_iter)) {
      v_new <- as.numeric(v %*% P)
      if (max(abs(v_new - v)) < tol) return(v_new)
      v <- v_new
    }
    warning("steady-state did not converge"); v
  }
  
  validate_transition_matrix <- function(P, tol = 1e-6) {
    if (!is.numeric(P) || nrow(P) != ncol(P)) stop("Invalid P")
    if (anyNA(P) || any(P < 0)) stop("P contains NA/negative")
    if (any(abs(rowSums(P) - 1) > tol)) stop("Rows do not sum to 1")
    TRUE
  }
  
  run_markov_generic <- function(mat, combs, strategy, params = list()) {
    nf <- sweep(mat, 2, colSums(mat), "/")
    comp_fun <- switch(
      strategy,
      baseline = function(i) {
        D <- sweep(nf, 2, nf[i, ], `-`)
        (rowSums(D > 0) + 0.5 * rowSums(D == 0)) / ncol(nf)
      },
      mc1 = function(i) {
        d <- sweep(nf, 2, nf[i, ], `-`)
        rowMeans(1/(1+exp(-params$sf*d)), na.rm = TRUE)
      },
      mc2 = function(i) {
        D <- sweep(nf, 2, nf[i, ], `-`)
        ((rowSums(D > 0) + 0.5 * rowSums(D == 0)) / ncol(nf)) + params$alpha
      },
      mc3 = function(i) {
        d <- sweep(nf, 2, nf[i, ], `-`)
        rowMeans(1/(1+exp(-params$w*d)), na.rm = TRUE)
      },
      bmarc = function(i) {
        sapply(seq_len(nrow(mat)), function(j)
          if(i!=j) sum(mat[i,] > mat[j,]) / ncol(mat) else 0
        )
      },
      bmarcw = function(i) {
        sapply(seq_len(nrow(mat)), function(j)
          mean(1/(1+exp(-2*(mat[i,]-mat[j,]))))
        )
      },
      stop("Unknown strategy: ", strategy)
    )
    P <- do.call(rbind, parLapply(cl, seq_len(nrow(nf)), comp_fun))
    diag(P) <- 0
    P <- sweep(P, 1, rowSums(P), "/")
    validate_transition_matrix(P)
    ss <- get_steady_state(P)
    tibble(Combination = combs, Score = ss, Rank = rank(ss, ties.method = "first"))
  }
  
  run_borda <- function(mat, combs) {
    scores <- rowSums(apply(mat, 2, rank))
    tibble(Combination = combs, Score = scores, Rank = rank(-scores, ties.method = "first"))
  }
  
  run_rankboost <- function(mat, combs) {
    lr <- 0.05; nfeat <- ncol(mat); wts <- rep(1/nfeat, nfeat)
    for (i in seq_len(500)) {
      prev <- wts
      sc <- mat %*% wts
      cors <- apply(mat, 2, function(x) cor(x, sc, method = "kendall", use = "complete.obs"))
      wts <- pmax(wts + lr * cors, 0); wts <- wts / sum(wts)
      if (max(abs(wts - prev)) < 1e-4) break
    }
    final <- mat %*% wts
    tibble(Combination = combs, Score = as.numeric(final), Rank = rank(-as.numeric(final), ties.method = "first"))
  }
  
  run_bayesmallows <- function(mat, combs) {
    mat2 <- t(mat)
    rank_mat <- t(apply(mat2, 1, function(x) rank(-x, ties.method = "average")))
    bm_data  <- BayesMallows::setup_rank_data(rank_mat)
    model_opts   <- BayesMallows::set_model_options(metric = "kendall")
    compute_opts <- BayesMallows::set_compute_options(nmc = 5000, burnin = 1000)
    priors       <- BayesMallows::set_priors()
    inits        <- BayesMallows::set_initial_values()
    progress     <- BayesMallows::set_progress_report(verbose = FALSE)
    model <- BayesMallows::compute_mallows(
      data            = bm_data,
      model_options   = model_opts,
      compute_options = compute_opts,
      priors          = priors,
      initial_values  = inits,
      progress_report = progress
    )
    burnin(model) <- compute_opts$burnin
    cp     <- BayesMallows::compute_consensus(model, type = "CP")
    scores <- cp$cumprob; scores[is.na(scores)] <- 0
    ranks  <- rank(-scores, ties.method = "first")
    tibble(Combination = combs, Score = scores, Rank = ranks)
  }
  
  run_birra <- function(mat, combs) {
    bins <- 5; iters <- 10; n <- nrow(mat); p <- ncol(mat)
    data_r <- apply(-mat, 2, rank) / n
    bdat <- ceiling(data_r * bins)
    guess <- rowMeans(mat)
    BFmat <- matrix(0, bins, p)
    for (iter in seq_len(iters)) {
      if (iter > 1 && cor(rank(guess), prev) > 1-1e-15) break
      prev <- guess
      ord <- order(guess)
      guess[ord[1:floor(0.05*n)]] <- 1
      guess[ord[(floor(0.05*n)+1):n]] <- 0
      for (j in seq_len(p)) for (b in seq_len(bins)) {
        sel <- which(bdat[, j] <= b)
        BFmat[b, j] <- log((sum(guess[sel])+1)/(sum(!guess[sel])+1)/(0.05/0.95))
      }
      BFmat <- apply(BFmat, 2, stats::smooth)
      BFmat <- apply(BFmat, 2, function(x) rev(cummax(rev(x))))
      guess <- rank(-rowSums(sapply(seq_len(p), function(j) BFmat[bdat[,j], j]), na.rm = TRUE))
    }
    inv_rank <- rank(-rowSums(sapply(seq_len(p), function(j) BFmat[bdat[,j], j]), na.rm = TRUE))
    inv_score <- length(inv_rank) + 1L - inv_rank
    tibble(Combination = combs, Score = inv_score, Rank = rank(-inv_score, ties.method = "first"))
  }
  
  run_cemc_k <- function(mat, combs) {
    rk <- apply(mat, 2, rank, ties.method = "average")
    kendals <- sapply(seq_len(ncol(rk)), function(k) {
      vals <- apply(rk, 1, function(r) sort(r)[k])
      mean(apply(rk, 2, function(fr) cor(rank(vals), fr, method = "kendall", use = "complete.obs")), na.rm = TRUE)
    })
    best <- which.max(kendals)
    vals <- apply(rk, 1, function(r) sort(r)[best])
    tibble(Combination = combs, Score = vals, Rank = rank(-vals, ties.method = "first"))
  }
  
  run_cemc_s <- function(mat, combs) {
    rk <- apply(mat, 2, rank, ties.method = "average")
    spears <- sapply(seq_len(ncol(rk)), function(u) {
      sth <- apply(rk, 1, function(r) mean(sort(r)[1:u]))
      mean(cor(sth, rk, method = "spearman"), na.rm = TRUE)
    })
    best <- which.max(spears)
    sth_final <- apply(rk, 1, function(r) mean(sort(r)[1:best]))
    tibble(Combination = combs, Score = sth_final, Rank = rank(-sth_final, ties.method = "first"))
  }
  
  run_rra <- function(mat, combs) {
    lists <- lapply(seq_len(ncol(mat)), function(i) combs[order(-mat[, i])])
    res <- aggregateRanks(glist = lists, N = length(combs), method = "RRA")
    tibble(Combination = res$Name, Score = res$Score, adjP = p.adjust(res$Score, method = "fdr"), Rank = rank(res$Score, ties.method = "first"))
  }
  
  methods <- list(
    Borda        = run_borda,
    Markov       = function(m,c) run_markov_generic(m,c,"baseline"),
    MC1          = function(m,c) run_markov_generic(m,c,"mc1", list(sf=nrow(m)^0.8)),
    MC2          = function(m,c) run_markov_generic(m,c,"mc2", list(alpha=0.01)),
    MC3          = function(m,c) run_markov_generic(m,c,"mc3", list(w=nrow(m)^0.8)),
    BMARC        = function(m,c) run_markov_generic(m,c,"bmarc"),
    BMARCW       = function(m,c) run_markov_generic(m,c,"bmarcw"),
    CEMC_k       = run_cemc_k,
    CEMC_s       = run_cemc_s,
    RankBoost    = run_rankboost,
    #BayesMallows = run_bayesmallows
    BIRRA        = run_birra,
    RRA          = run_rra
  )
  
  data <- read_csv(csv_file, show_col_types = FALSE)
  data$SynergyRank <- rank(-data$Synergy, ties.method = "average")
  folder <- tools::file_path_sans_ext(basename(csv_file))
  dir.create(folder, showWarnings = FALSE)
  
  data <- data %>% select(-Synergy)
  metadata_cols <- c("Drug1", "Drug2", "Cell_Line", "Tissue", "SynergyRank")
  feats <- setdiff(names(data), metadata_cols)
  data$Combination <- apply(data[feats], 1, paste, collapse = "_")
  
  all_rt <- list()
  set.seed(42)
  for (n in sample_sizes) {
    cat("--- Running sample size: ", n, "\n"); flush.console()
    sub <- if (n == nrow(data)) data else slice_sample(data, n = n)
    mat <- as.matrix(sub[feats])
    combs <- sub$Combination
    rt <- tibble(Method = names(methods), Time = NA_real_)
    out <- vector("list", length(methods)); names(out) <- names(methods)
    
    for (m in names(methods)) {
      cat(" Method: ", m, "\n"); flush.console()
      t0 <- system.time({
        out[[m]] <- tryCatch(
          methods[[m]](mat, combs),
          error = function(e) {
            warning(sprintf("%s failed: %s", m, e$message))
            tibble(Combination = combs, Score = NA_real_, Rank = NA_integer_)
          }
        )
      })["elapsed"]
      rt <- rt %>% mutate(Time = if_else(Method == m, as.numeric(t0), Time))
      cat("  -> ", m, " done in ", round(as.numeric(t0),3), "s\n"); flush.console()
    }
    
    all_rt[[as.character(n)]] <- rt
    
    valid <- keep(out, ~ !all(is.na(.x$Rank)))
    combined <- reduce(
      map2(valid, names(valid), ~ .x %>%
             select(Combination, Rank) %>%
             rename(!!.y := Rank)),
      inner_join, by = "Combination"
    )
    
    combined <- inner_join(
      sub %>% select(Combination, Drug1, Drug2, Cell_Line, Tissue, SynergyRank),
      combined,
      by = "Combination"
    ) %>%
      select(Drug1, Drug2, Cell_Line, Tissue, SynergyRank, everything(), -Combination)
    
    write_csv(combined, file.path(folder, paste0("Combined_Ranks_", n, ".csv")))
    
    markov_methods <- c("Markov","MC1","MC2","MC3","BMARC","BMARCW")
    for(m in markov_methods){
      rho <- cor(combined[[m]], combined$Borda, method = "spearman", use = "complete.obs")
      if(!is.na(rho) && rho < 0){
        combined[[m]] <- max(combined[[m]], na.rm=TRUE) + 1L - combined[[m]]
      }
    }
    
    Mr <- as.matrix(combined %>% select(-Drug1, -Drug2, -Cell_Line, -Tissue))
    corr_s <- cor(Mr, method = "spearman")
    write_csv(as.data.frame(corr_s) %>% rownames_to_column("Method1"),
              file.path(folder, paste0("Spearman_Rank_Corr_", n, ".csv")))
  }
  
  result <- bind_rows(all_rt, .id = "SampleSize") %>%
    pivot_wider(names_from = SampleSize, values_from = Time)
  write_csv(result, file.path(folder, "Runtimes_All.csv"))
  beep()
  cat("Done\n")
  invisible(result)
}

# full_n <- nrow(read_csv("Generic.csv", show_col_types = FALSE))
# 
# compare_ranks(
#   "Generic.csv",
#   sample_sizes = full_n
# )

#Run Merck 2-5
{fold_files <- paste0("Merck_Fold", 1:5, ".csv")
purrr::walk(fold_files, function(f) {
  full_n <- nrow(readr::read_csv(f, show_col_types = FALSE))
  compare_ranks(f, sample_sizes = full_n)
})}


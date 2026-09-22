# z-scores from all simulations


# Load libraries ----------------------------------------------------------
# check libraries carfully. important on the server/cluster
.libPaths("/home/tietje/mamba/envs/r453/lib/R/library")
rm(list=(ls())) # clear workspace
library(data.table)
library(parallel)
library(RhpcBLASctl)
print(paste("System time: ", Sys.time()))





# Bash input ###############################
(args <- commandArgs(trailingOnly=TRUE))
ncores       = as.numeric(args[1])
blas_thred   = as.numeric(args[2])
omp_thred    = as.numeric(args[3])
pair_select  = args[4]        # "all" or a subset flag 
pairs_n      = as.numeric(args[5])

# restrict this to avoid thread oversubscription
blas_set_num_threads(blas_thred)
omp_set_num_threads(omp_thred)


# load_pairs --------------------------------------------------------------
# Build the pair list directly from pairwise-analysis output so this
# server job can run before the main analysis script.
result_files <- dir("data/counterfactuals2/results", full.names = TRUE,
                    pattern = "results_.*chunk")
result_list <- lapply(result_files, readRDS)
result_dt <- rbindlist(result_list, fill = TRUE, use.names = TRUE)
combos <- unique(result_dt[, .(dataset_id, species_pair)])
rm(result_files, result_list, result_dt)
# subset?
if (pair_select != "all") {
  combos <- combos[1:pairs_n, ]
}

# function_prep -----------------------------------------------------------
n_shuffle <- 100
n_cols    <- 100




## ---- simple timestamped logging to stdout (captured by bash's redirect) ----
log_msg <- function(...) {
  cat(sprintf("[%s] [pid:%d] %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              Sys.getpid(),
              paste0(..., collapse = "")))
  flush(stdout())
}

col_cor <- function(x, y) {
  xc  <- sweep(x, 2, colMeans(x))
  yc  <- sweep(y, 2, colMeans(y))
  num <- colSums(xc * yc)
  den <- sqrt(colSums(xc^2) * colSums(yc^2))
  num / den
}

process_pair <- function(i) {
  atl  <- combos$dataset_id[i]
  pair <- combos$species_pair[i]
  t0   <- Sys.time()
  
  result <- tryCatch({
    sp1 <- sub("\\|.*", "", pair)
    sp2 <- sub(".*\\|", "", pair)
    
    sim1_list <- readRDS(sprintf("data/uncertainty/simulated_%s_%s.rds", sp1, atl))
    sim2_list <- readRDS(sprintf("data/uncertainty/simulated_%s_%s.rds", sp2, atl))
    
    n_t     <- length(sim1_list)
    t_idx   <- unique(c(1, n_t))
    t_label <- if (length(t_idx) == 1) "only" else c("first", "last")
    
    out <- vector("list", length(t_idx))
    for (j in seq_along(t_idx)) {
      t    <- t_idx[j]
      sim1 <- sim1_list[[t]]
      sim2 <- sim2_list[[t]]
      
      rank1 <- apply(sim1, 2, rank)
      rank2 <- apply(sim2, 2, rank)
      obs_cor <- col_cor(rank1, rank2)
      
      cors <- matrix(NA_real_, n_shuffle, n_cols)
      for (k in seq_len(n_shuffle)) {
        s1 <- apply(rank1, 2, sample)
        s2 <- apply(rank2, 2, sample)
        cors[k, ] <- col_cor(s1, s2)
      }
      
      z_score <- (obs_cor - colMeans(cors)) / apply(cors, 2, sd)
      
      out[[j]] <- list(dataset_id = atl, species_pair = pair,
                       timestep = t, timestep_label = t_label[j],
                       z_scores = z_score)
    }
    out
  }, error = function(e) {
    log_msg(sprintf("ERROR pair %d (%s | %s): %s", i, atl, pair, conditionMessage(e)))
    NULL
  })
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.null(result)) {
    log_msg(sprintf("OK pair %d/%d (%s | %s) - %.1f sec", i, nrow(combos), atl, pair, elapsed))
  }
  
  result
}

log_msg(sprintf("Starting: %d pairs, %d shuffles, mc.cores=%d, blas_threads=%d, omp_threads=%d",
                nrow(combos), n_shuffle, ncores, blas_thred, omp_thred))
run_start <- Sys.time()

results_nested <- mclapply(seq_len(nrow(combos)), process_pair, mc.cores = ncores)

n_failed <- sum(sapply(results_nested, is.null))
log_msg(sprintf("Finished: %d/%d succeeded, %d failed - total time %.1f min",
                nrow(combos) - n_failed, nrow(combos), n_failed,
                as.numeric(difftime(Sys.time(), run_start, units = "mins"))))

results <- unlist(results_nested[!sapply(results_nested, is.null)], recursive = FALSE)
saveRDS(results, "data/z_scores.rds")
Sys.time()









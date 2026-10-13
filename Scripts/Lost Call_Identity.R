# =======================================
# Lost Call Individual Identity Analysis
# White-faced capuchins (Cebus imitator)
# Chapoy et al., 2026
# ======================================

#Packages
suppressPackageStartupMessages({
  library(MASS)        # lda
  library(dplyr)
  library(readr)
  library(stringr)
})

# Configuration ----------
CSV_PATH     <- "individual_analyses.csv"
GROUP_COL <- "individual"
EXCLUDE_COLS <- c("rating", "Sequence position", "Sequence call length")
FEATURE_LIST <- NULL
dir.create("results", showWarnings = FALSE) #create output folder
set.seed(42) #reproducibility

# Helper functions####
chance_level <- function(y_factor) {
  # chance level for balanced classes
  1 / nlevels(y_factor)
}

balanced_accuracy <- function(y_true, y_pred) {
  tbl <- table(y_true, y_pred)
  lv <- sort(unique(c(rownames(tbl), colnames(tbl))))
  tbl <- tbl[lv, lv, drop = FALSE]
  rec <- sapply(lv, function(k) {
    denom <- sum(tbl[k, , drop = FALSE])
    if (denom == 0) return(NA_real_)
    tbl[k, k] / denom
  })
  mean(rec, na.rm = TRUE)
}

per_class_metrics <- function(y_true, y_pred) {
  tbl <- table(True = y_true, Pred = y_pred)
  classes <- sort(unique(c(rownames(tbl), colnames(tbl))))
  
  out <- lapply(classes, function(k) {
    tp <- if (k %in% rownames(tbl) && k %in% colnames(tbl)) tbl[k, k] else 0
    fn <- if (k %in% rownames(tbl)) sum(tbl[k, , drop = FALSE]) - tp else 0
    fp <- if (k %in% colnames(tbl)) sum(tbl[, k, drop = FALSE]) - tp else 0
    tn <- sum(tbl) - tp - fn - fp
    
    precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
    recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
    f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA_real_
    else 2 * precision * recall / (precision + recall)
    
    data.frame(
      class = as.character(k),
      n_true = if (k %in% rownames(tbl)) sum(tbl[k, , drop = FALSE]) else 0,
      precision = precision,
      recall = recall,
      f1 = f1,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(out)
}


# Load & prep ####
df <- read_csv(CSV_PATH)

if (!(GROUP_COL %in% names(df))) stop(sprintf("Grouping column '%s' not found.", GROUP_COL))
df[[GROUP_COL]] <- as.factor(df[[GROUP_COL]])

# If explicit list provided, use it 
if (!is.null(FEATURE_LIST)) {
  feature_cols <- setdiff(intersect(FEATURE_LIST, names(df)), EXCLUDE_COLS)
} else {
  #default
  num_cols <- names(df)[sapply(df, is.numeric)]
  feature_cols <- setdiff(num_cols, EXCLUDE_COLS)
}
feature_cols <- setdiff(feature_cols, GROUP_COL)
if (!length(feature_cols)) stop("No numeric feature columns after exclusions.")

# Keep complete cases
keep <- complete.cases(df[, c(feature_cols, GROUP_COL), drop = FALSE])
if (!all(keep)) message(sprintf("Dropping %d rows with NA in features/group.", sum(!keep)))
dat <- df[keep, c(feature_cols, GROUP_COL), drop = FALSE]
y   <- droplevels(dat[[GROUP_COL]])
X   <- dat[, feature_cols, drop = FALSE]

# Scale features
X_scaled_only <- as.data.frame(scale(X))
X_scaled <- X_scaled_only
X_scaled[[GROUP_COL]] <- y

cat(sprintf("N=%d, classes=%d, chance=%.1f%%, features=%d\n",
            nrow(dat), nlevels(y), 100*chance_level(y), length(feature_cols)))
cat("Features: ", paste(feature_cols, collapse=", "), "\n", sep="")

# MANOVA (Wilks' λ) ####
# Rank check (should equal #features)
Xmat <- as.matrix(X_scaled[, feature_cols, drop = FALSE])
if (qr(Xmat)$rank < ncol(Xmat)) {
  warning("Design still rank-deficient; consider dropping another derived/constant feature.")
}

man <- manova(as.matrix(X_scaled[, feature_cols, drop = FALSE]) ~ y)

wilks <- try(summary(man, test = "Wilks"), silent = TRUE)
if (inherits(wilks, "try-error")) {
  message("Wilks' λ failed (rank deficiency). Falling back to Pillai's trace.")
  wilks <- summary(man, test = "Pillai")
}

capture.output(wilks, file = "results/manova_test.txt")
cat("Saved: results/manova_test.txt\n")

# LDA (eigenvalues + accuracies) ####
lda_fit <- lda(x = X_scaled[, feature_cols, drop = FALSE], grouping = y)

# Eigenvalues
eigvals <- lda_fit$svd^2
cc      <- lda_fit$svd / sqrt(1 + lda_fit$svd^2)
pct     <- eigvals / sum(eigvals)

eig_df <- data.frame(
  function_index = seq_along(eigvals),
  eigenvalue = eigvals,
  canonical_correlation = cc,
  pct_variance = pct,
  cum_pct_variance = cumsum(pct)
)
write_csv(eig_df, "results/lda_eigenvalues.csv")
cat("Saved: results/lda_eigenvalues.csv\n")
head(eig_df, 10)

# Resubstitution (inflated)
y_hat_fit <- predict(lda_fit)$class
acc_resub <- mean(y_hat_fit == y)
bal_resub <- balanced_accuracy(y, y_hat_fit)
resub_metrics <- per_class_metrics(y, y_hat_fit)
write_csv(resub_metrics, "results/resub_per_class_metrics.csv")
cat(sprintf("Resubstitution accuracy: %.3f | Balanced: %.3f\n", acc_resub, bal_resub))
cat("Saved: results/resub_per_class_metrics.csv\n")

# LOOCV (honest)
lda_cv <- lda(x = X_scaled[, feature_cols, drop = FALSE], grouping = y, CV = TRUE)
y_hat_cv <- lda_cv$class
acc_cv <- mean(y_hat_cv == y)
bal_cv <- balanced_accuracy(y, y_hat_cv)
cv_metrics <- per_class_metrics(y, y_hat_cv)
write_csv(cv_metrics, "results/loocv_per_class_metrics.csv")
cat(sprintf("LOOCV accuracy: %.3f | Balanced: %.3f\n", acc_cv, bal_cv))
cat("Saved: results/loocv_per_class_metrics.csv\n")

cv_metrics %>%
  arrange(recall) %>%
  select(class, n_true, recall, precision, f1) %>%
  as.data.frame() %>%
  print()


#bootstrap CI for LOOCV accuracy ####
B <- 2000
set.seed(42)

boot_bal <- numeric(B)
boot_raw <- numeric(B)

n <- nrow(X_scaled)

for (b in seq_len(B)) {
  
  idx <- sample(seq_len(n), n, replace = TRUE)
  
  X_b <- X_scaled[idx, feature_cols, drop = FALSE]
  y_b <- y[idx]
  
  fit_b <- try(lda(x = X_b, grouping = y_b, CV = TRUE), silent = TRUE)
  
  if (inherits(fit_b, "try-error")) {
    boot_raw[b] <- NA_real_
    boot_bal[b] <- NA_real_
    next
  }
  
  yhat_b <- fit_b$class
  boot_raw[b] <- mean(yhat_b == y_b)
  boot_bal[b] <- balanced_accuracy(y_b, yhat_b)
}

boot_bal_ok <- boot_bal[!is.na(boot_bal)]
boot_raw_ok <- boot_raw[!is.na(boot_raw)]

ci_bal <- quantile(boot_bal_ok, c(0.025, 0.975))
ci_raw <- quantile(boot_raw_ok, c(0.025, 0.975))

ci_df <- data.frame(
  metric = c("raw", "balanced"),
  observed = c(acc_cv, bal_cv),
  ci_lower = c(ci_raw[1], ci_bal[1]),
  ci_upper = c(ci_raw[2], ci_bal[2])
)

print(ci_df)

#summary ####
cat("\n=== SUMMARY ===\n")
cat("Groups: ", paste(levels(y), collapse=", "), "\n", sep="")
cat("Features (", length(feature_cols), "): ", paste(feature_cols, collapse=", "), "\n", sep="")
cat("Chance level: ", sprintf("%.1f%%", 100 / nlevels(y)), "\n", sep="")
cat("Resubstitution accuracy: ", sprintf("%.1f%%", 100*acc_resub),
    " (balanced ", sprintf("%.1f%%", 100*bal_resub), ")\n", sep="")
cat("LOOCV accuracy:          ", sprintf("%.1f%%", 100*acc_cv),
    " (balanced ", sprintf("%.1f%%", 100*bal_cv), ")\n", sep="")
cat("Top eigenvalues:\n")
print(head(eig_df, 5))
cat("\nMANOVA results saved to results/manova_test.txt\n")
cat(sprintf("Individuals: %d | Calls: %d\n", nlevels(y), length(y)))


# Permutation test (pDFA) for non-independence####
# Tests whether LOOCV LDA accuracy exceeds chance under label permutation

N_PERM    <- 1000 # number of permutations
BLOCK_COL <- NULL 

# helper functions
permute_within_blocks <- function(labels, blocks) {
  labs_chr <- as.character(labels)
  out <- labs_chr
  for (b in levels(blocks)) {
    idx <- which(blocks == b)
    if (length(idx) > 1) {
      out[idx] <- sample(labs_chr[idx], length(idx), replace = FALSE)
    }
  }
  factor(out, levels = levels(labels))
}

compute_loocv_acc <- function(Xs, y_factor) {
  fit_cv <- try(
    lda(x = Xs, grouping = y_factor, CV = TRUE),
    silent = TRUE
  )
  if (inherits(fit_cv, "try-error")) {
    return(c(raw = NA_real_, bal = NA_real_))
  }
  y_hat <- fit_cv$class
  c(
    raw = mean(y_hat == y_factor),
    bal = balanced_accuracy(y_factor, y_hat)
  )
}

# observed accuracies 
obs_raw <- acc_cv
obs_bal <- bal_cv

# optional blocking vector
block_vec <- NULL
if (!is.null(BLOCK_COL)) {
  block_vec <- droplevels(as.factor(df[[BLOCK_COL]][keep]))
  if (length(block_vec) != length(y)) {
    stop("Blocking vector length mismatch after NA filtering.")
  }
}

# permutation loop
null_raw <- numeric(N_PERM)
null_bal <- numeric(N_PERM)

for (i in seq_len(N_PERM)) {
  y_perm <- if (is.null(BLOCK_COL)) {
    factor(
      sample(as.character(y), length(y), replace = FALSE),
      levels = levels(y)
    )
  } else {
    permute_within_blocks(y, block_vec)
  }
  
  accs <- compute_loocv_acc(
    X_scaled[, feature_cols, drop = FALSE],
    y_perm
  )
  
  null_raw[i] <- accs["raw"]
  null_bal[i] <- accs["bal"]
}

#remove failed permutations
null_raw_ok <- null_raw[!is.na(null_raw)]
null_bal_ok <- null_bal[!is.na(null_bal)]

# p-values (add-one correction)
p_raw <- (sum(null_raw_ok >= obs_raw) + 1) / (length(null_raw_ok) + 1)
p_bal <- (sum(null_bal_ok >= obs_bal) + 1) / (length(null_bal_ok) + 1)

# save outputs
pd_df <- data.frame(
  metric   = c("raw", "balanced"),
  observed = c(obs_raw, obs_bal),
  p_value  = c(p_raw, p_bal)
)

write_csv(
  data.frame(null_raw = null_raw, null_bal = null_bal),
  "results/pDFA_null_accuracies.csv"
)

write_csv(
  pd_df,
  "results/pDFA_summary.csv"
)

#console summary 
cat("\n=== Permutation test (pDFA-style) ===\n")
if (is.null(BLOCK_COL)) {
  cat("Permutation scheme: simple label permutation\n")
} else {
  cat(sprintf(
    "Permutation scheme: labels permuted within blocks of '%s'\n",
    BLOCK_COL
  ))
}
cat(sprintf(
  "Observed LOOCV accuracy: raw = %.3f, balanced = %.3f\n",
  obs_raw, obs_bal
))
cat(sprintf(
  "Permutation p-values:   raw = %.4f, balanced = %.4f\n",
  p_raw, p_bal
))
cat(sprintf(
  "Permutations retained:  raw = %d, balanced = %d (of %d)\n",
  length(null_raw_ok), length(null_bal_ok), N_PERM
))
cat("Saved: results/pDFA_null_accuracies.csv, results/pDFA_summary.csv\n")



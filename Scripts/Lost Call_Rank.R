# =======================================
# Lost Call Rank Analysis
# White-faced capuchins (Cebus imitator)
# Chapoy et al., 2026
# ======================================

## Packages#### 
suppressPackageStartupMessages({
  library(MASS)        # lda
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(performance)
  library(car)
  library(nlme)
  library(see)
  library(DHARMa)
  library(tidyr)
  library(effects)
})

set.seed(42)


# data cleaning####
df_rank <- read_csv("data/rank_lmm_dataset.csv")
unique(df_rank$rank_clean)

df_rank$rank_clean <- factor(df_rank$rank_clean, levels = c("subordinate", "alpha"))
contrasts(df_rank$rank_clean)

df_rank %>%
  distinct(individual, rank_clean) %>%
  count(rank_clean)

#Formant dispersion model#### 
# quick distribution checks
hist(df_rank$formant_dispersion, breaks = 30)
qqnorm(df_rank$formant_dispersion); qqline(df_rank$formant_dispersion, col="red")

# sample size
df_rank %>% count(individual, rank_clean)

# fit formant dispersion model
m_disp <- lmer(formant_dispersion ~ rank_clean + (1 | individual), data = df_rank)
summary(m_disp)
drop1(m_disp, test = "Chisq")   # optional, but consistent
emmeans(m_disp, ~ rank_clean)
r2(m_disp)

m1_ml <- lmer(formant_dispersion ~ rank_clean + (1 | individual),
              data = df_rank, REML = FALSE)

m2_ml <- lmer(formant_dispersion ~ rank_clean + Site + (1 | individual),
              data = df_rank, REML = FALSE)

AIC(m1_ml, m2_ml)
anova(m1_ml, m2_ml)

# diagnostics for formant dispersion
sim_disp <- simulateResiduals(m_disp, n = 1000)
plot(sim_disp)
testUniformity(sim_disp)
testDispersion(sim_disp)
testOutliers(sim_disp)
plot(allEffects(m_disp))

# Perform Kolmogorov-Smirnov (KS) test
ks_test_disp <- testUniformity(sim_disp)
# Print KS test result
print(ks_test_disp) #p=0.4811

# Perform dispersion test
dispersion_test_disp <- testDispersion(sim_disp)
# Print dispersion test result
print(dispersion_test_disp) #p=0.726

#calculate means for each ind. 
df_ind_means_disp <- df_rank %>%
  group_by(individual, rank_clean) %>%
  summarise(
    mean_disp = mean(formant_dispersion, na.rm = TRUE),
    n_calls = n(),
    .groups = "drop"
  )

# estimated marginal means from LMM
library(emmeans)
emm_disp <- emmeans(m_disp, ~ rank_clean)
emm_df_disp <- as.data.frame(emm_disp)

#plot- FINAL PLOT IN PAPER!!
p_disp <- ggplot(df_ind_means_disp, aes(rank_clean, mean_disp)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.7, size = 2) +
  geom_pointrange(
    data = emm_df_disp,
    aes(
      x = rank_clean,
      y = emmean,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    inherit.aes = FALSE,
    size = 1.1,
    color = "black"
  ) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_discrete(labels = c("subordinate" = "Subordinate", "alpha" = "Alpha")) +
  theme_minimal() +
  labs(
    y = "Individual mean formant dispersion (Hz)",
    x = NULL,
    title = "Formant dispersion"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
  
#Fundamental Frequency model####
# quick distribution checks
hist(df_rank$F0, breaks = 30)
shapiro.test(df_rank$F0)
qqnorm(df_rank$F0); qqline(df_rank$F0, col = "red")

# Log-transform F0 because residual diagnostics for raw F0 indicated
# departures from model assumptions.
df_rank$logF0 <- log(df_rank$F0)

## Final model
m_f0_logn <- lmer(logF0 ~ rank_clean + (1 | individual), data = df_rank)
summary(m_f0_logn)
drop1(m_f0_logn, test = "Chisq")

# Estimated marginal means, back-transformed to Hz
emmeans(m_f0_logn, ~ rank_clean, type = "response")

m1_f0_ml <- lmer(logF0 ~ rank_clean + (1 | individual),
                 data = df_rank, REML = FALSE)

m2_f0_ml <- lmer(logF0 ~ rank_clean + Site + (1 | individual),
                 data = df_rank, REML = FALSE)

AIC(m1_f0_ml, m2_f0_ml)
anova(m1_f0_ml, m2_f0_ml)

# Diagnostics
sim_f0 <- simulateResiduals(m_f0_logn, n = 1000)
plot(sim_f0)

ks_test_f0 <- testUniformity(sim_f0)
print(ks_test_f0)

dispersion_test_f0 <- testDispersion(sim_f0)
print(dispersion_test_f0)

outlier_test_f0 <- testOutliers(sim_f0)
print(outlier_test_f0)
res <- resid(m_f0_logn)
fitted <- fitted(m_f0_logn)


#Model R2
r2_f0 <- r2(m_f0_logn)
r2_f0

# calculate individual means for plotting
df_ind_means_f0 <- df_rank %>%
  group_by(individual, rank_clean) %>%
  summarise(
    mean_F0 = mean(F0, na.rm = TRUE),
    n_calls = n(),
    .groups = "drop"
  )

# estimated marginal means from logF0 model
emm_f0 <- emmeans(m_f0_logn, ~ rank_clean)
emm_df_f0 <- as.data.frame(emm_f0)

# back transform to Hz
emm_df_f0 <- emm_df_f0 %>%
  mutate(
    emmean_bt = exp(emmean),
    lower.CL_bt = exp(lower.CL),
    upper.CL_bt = exp(upper.CL)
  )

# plot
p_f0 <- ggplot(df_ind_means_f0, aes(rank_clean, mean_F0)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.7, size = 2) +
  geom_pointrange(
    data = emm_df_f0,
    aes(
      x = rank_clean,
      y = emmean_bt,
      ymin = lower.CL_bt,
      ymax = upper.CL_bt
    ),
    inherit.aes = FALSE,
    size = 1.1,
    color = "black"
  ) +
  scale_x_discrete(labels = c("subordinate" = "Subordinate", "alpha" = "Alpha")) +
  theme_minimal() +
  labs(
    y = "Individual mean fundamental frequency (Hz)",
    x = NULL,
    title = "Fundamental frequency"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# Call Duration Model ####

# quick distribution checks
hist(df_rank$call_length, breaks = 30)
qqnorm(df_rank$call_length); qqline(df_rank$call_length, col = "red")

# final call duration model
m_len <- lmer(call_length ~ rank_clean + (1 | individual), data = df_rank)
summary(m_len)
drop1(m_len, test = "Chisq")   # optional, but consistent
emmeans(m_len, ~ rank_clean)
r2(m_len)

m1_len_ml <- lmer(call_length ~ rank_clean + (1 | individual),
                  data = df_rank, REML = FALSE)

m2_len_ml <- lmer(call_length ~ rank_clean + Site + (1 | individual),
                  data = df_rank, REML = FALSE)

AIC(m1_len_ml, m2_len_ml)
anova(m1_len_ml, m2_len_ml)

# diagnostics for call length
sim_len <- simulateResiduals(m_len, n = 1000)
plot(sim_len)
testUniformity(sim_len)
testDispersion(sim_len)
testOutliers(sim_len)

# Perform Kolmogorov-Smirnov (KS) test
ks_test_len <- testUniformity(sim_len)
# Print KS test result
print(ks_test_len) #p=0.222

# Perform dispersion test
dispersion_test_len <- testDispersion(sim_len)
# Print dispersion test result
print(dispersion_test_len) #p=0.668

# calculate individual means for plotting
df_ind_means_len <- df_rank %>%
  group_by(individual, rank_clean) %>%
  summarise(
    mean_len = mean(call_length, na.rm = TRUE),
    n_calls = n(),
    .groups = "drop"
  )

# estimated marginal means from LMM
emm_len <- emmeans(m_len, ~ rank_clean)
emm_df_len <- as.data.frame(emm_len)

# plot
p_len <- ggplot(df_ind_means_len, aes(rank_clean, mean_len)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.7, size = 2) +
  geom_pointrange(
    data = emm_df_len,
    aes(
      x = rank_clean,
      y = emmean,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    inherit.aes = FALSE,
    size = 1.1,
    color = "black"
  ) +
  scale_x_discrete(labels = c("subordinate" = "Subordinate", "alpha" = "Alpha")) +
  theme_minimal() +
  labs(
    y = "Individual mean call length (s)",
    x = NULL,
    title = "Call length"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )


#Patch plots together####
library(patchwork)

final_plot <- (p_disp | p_f0 | p_len) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0, 1)
  )
#save plot
ggsave(
  "Figure_rank_effects.pdf",
  plot = final_plot,
  width = 10,
  height = 4,
  units = "in"
)

ggsave(
  "Figure_rank_effects.png",
  plot = final_plot,
  width = 10,
  height = 4,
  units = "in",
  dpi = 600
)


#pDFA for dominance rank (LOIO LDA + individual-level permutations)####

# 1) Prep data
df_pdfa <- df_rank

feat_cols <- c(
  "F0","sd_pitch","min_pitch","max_pitch",
  "pitch_instability","amplitude_instability",
  "F1","F2",
  "formant_dispersion","harmonic-noise-ratio",
  "call_length")

# keep complete cases for the multivariate model
df_pdfa <- df_pdfa %>%
  drop_na(all_of(feat_cols)) %>%
  mutate(rank_clean = droplevels(factor(rank_clean)))

#LOIO LDA function (with fold-wise standardization)
loio_lda_metrics <- function(dat, feat_cols) {
  
  inds <- unique(dat$individual)
  
  all_truth <- character(0)
  all_pred  <- character(0)
  
  for (i in inds) {
    
    train <- dat[dat$individual != i, , drop = FALSE]
    test  <- dat[dat$individual == i, , drop = FALSE]
    
    x_train <- as.matrix(train[, feat_cols, drop = FALSE])
    x_test  <- as.matrix(test[,  feat_cols, drop = FALSE])
    
    # fold-wise z standardization (fit on train, apply to test)
    x_train_scaled <- scale(x_train)
    
    # handle any zero-variance columns in this fold
    sc <- attr(x_train_scaled, "scaled:scale")
    keep <- is.finite(sc) & sc > 0
    if (!all(keep)) {
      x_train_scaled <- x_train_scaled[, keep, drop = FALSE]
      x_test <- x_test[, keep, drop = FALSE]
      # recompute scaling with the reduced set
      x_train_scaled <- scale(x_train_scaled)
    }
    
    x_test_scaled <- scale(
      x_test,
      center = attr(x_train_scaled, "scaled:center"),
      scale  = attr(x_train_scaled, "scaled:scale")
    )
    
    fit  <- lda(x = x_train_scaled, grouping = train$rank_clean)
    pred <- predict(fit, newdata = x_test_scaled)$class
    
    all_truth <- c(all_truth, as.character(test$rank_clean))
    all_pred  <- c(all_pred,  as.character(pred))
  }
  
  all_truth <- factor(all_truth, levels = levels(dat$rank_clean))
  all_pred  <- factor(all_pred,  levels = levels(dat$rank_clean))
  
  cm <- table(all_truth, all_pred)
  
  acc <- sum(diag(cm)) / sum(cm)
  recalls <- diag(cm) / rowSums(cm)
  bal_acc <- mean(recalls, na.rm = TRUE)
  
  list(cm = cm, acc = acc, bal_acc = bal_acc, recalls = recalls)
}

#Observed performance
obs <- loio_lda_metrics(df_pdfa, feat_cols)

cat("\nObserved confusion matrix:\n")
print(obs$cm)
cat(sprintf("\nObserved LOIO accuracy: %.3f\n", obs$acc))
cat(sprintf("Observed LOIO balanced accuracy: %.3f\n", obs$bal_acc))
cat("\nClass recalls:\n")
print(obs$recalls)

#Permutation test (permute rank at the individual level)
set.seed(42)
B <- 1000 

# map each individual to its observed rank
ind_rank <- df_pdfa %>% distinct(individual, rank_clean)

perm_acc <- numeric(B)
perm_bal <- numeric(B)

for (b in seq_len(B)) {
  
  # permute ranks across individuals
  permuted <- sample(ind_rank$rank_clean, size = nrow(ind_rank), replace = FALSE)
  
  ind_rank_perm <- ind_rank %>%
    mutate(rank_perm = permuted) %>%
    select(individual, rank_perm)
  
  # assign permuted rank back to all calls from that individual
  dat_perm <- df_pdfa %>%
    left_join(ind_rank_perm, by = "individual") %>%
    mutate(rank_clean = factor(rank_perm, levels = levels(df_pdfa$rank_clean))) %>%
    select(-rank_perm)
  
  res <- loio_lda_metrics(dat_perm, feat_cols)
  perm_acc[b] <- res$acc
  perm_bal[b] <- res$bal_acc
}

# empirical p values (add-one correction)

p_acc <- (sum(perm_acc >= obs$acc) + 1) / (B + 1)
p_bal <- (sum(perm_bal >= obs$bal_acc) + 1) / (B + 1)

cat(sprintf("\nPermutation p value (accuracy): %.4f\n", p_acc))
cat(sprintf("Permutation p value (balanced accuracy): %.4f\n", p_bal))




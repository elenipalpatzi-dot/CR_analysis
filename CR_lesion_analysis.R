
# Analysis code is for:
# "Lesion burden limits the association of cognitive reserve 
#  with functional stroke outcome"
# Palpatzis, Gallucci, Sperber, Hakim & Umarova 2026
# R 4.5.0
#
# Data not included (due to ethical restrictions).
# Contact corresponding author for access.



library(tidyverse)
library(MASS)
library(pROC)
library(logistf)
library(ggeffects)
library(ResourceSelection)
library(car)


#  1. Data prep 

# log-transform lesion volume, rescale age, prepare outcomes
data$log_lesion <- log(data$Volume_roi_cL + 1)
data$age_decade <- data$age / 10
data$nihss_1 <- data$nihss_24h + 1      
data$mrs_dichot <- ifelse(data$mrs_fu >= 3, 1, 0)

# tertiles on full sample
data$lesion_t3 <- cut(data$Volume_roi_cL,
                    breaks = quantile(data$Volume_roi_cL, c(0, 1/3, 2/3, 1)),
                    labels = c("Small", "Medium", "Large"),
                    include.lowest = TRUE)

# functional outcome sample, need mrs + prestroke_mrs
data_mrs <- data %>% 
  filter(!is.na(mrs_fu), !is.na(prestroke_mrs),
         !record_id %in% recurrent_ids)


# stroke severity sample, need NIHSS
data_nihss <- data %>% filter(!is.na(nihss_24h))

cat("Full sample:", nrow(data), "\n")
cat("mRS sample:", nrow(data_mrs), "| events:", sum(data_mrs$mrs_dichot), "\n")
cat("NIHSS sample:", nrow(data_nihss), "\n")

# mean-centering the variables that go into interactions
data_mrs$yoe_c <- data_mrs$yoe - mean(data_mrs$yoe)
data_mrs$log_lesion_c <- data_mrs$log_lesion - mean(data_mrs$log_lesion)
data_mrs$criq_work_10 <- data_mrs$criq_work / 10

data_nihss$yoe_c <- data_nihss$yoe - mean(data_nihss$yoe)
data_nihss$log_lesion_c <- data_nihss$log_lesion - mean(data_nihss$log_lesion)
data_nihss$criq_work_10 <- data_nihss$criq_work / 10
data_nihss$criq_work_10_c <- data_nihss$criq_work_10 - mean(data_nihss$criq_work_10)

cat("Centring values — YoE:", round(mean(data_mrs$yoe), 1),
    "| log_lesion:", round(mean(data_mrs$log_lesion), 2),
    "| CRIq-Work:", round(mean(data_mrs$criq_work), 1), "\n")



# 2. Descriptives 
cat("Age:", median(data$age), "[", quantile(data$age, .25), "-", quantile(data$age, .75), "]\n")
cat("Female:", sum(data$biol_sex == "Female"), "(", round(mean(data$biol_sex == "Female") * 100, 1), "%)\n")
cat("YoE:", median(data$yoe), "[", quantile(data$yoe, .25), "-", quantile(data$yoe, .75), "]\n")
cat("CRIq-Work:", median(data$criq_work), "[", quantile(data$criq_work, .25), "-", quantile(data$criq_work, .75), "]\n")
cat("Lesion cL:", round(median(data$Volume_roi_cL), 2), "\n")
cat("NIHSS 24h:", median(data$nihss_24h, na.rm = TRUE), "\n")

# CR proxy correlations 
cor.test(data_mrs$yoe, data_mrs$criq_work, method = "pearson")
cor.test(data_mrs$yoe, data_mrs$criq_work, method = "kendall")
cor.test(data_mrs$yoe, data_mrs$prestroke_mrs, method = "spearman")
cor.test(data_mrs$criq_work, data_mrs$prestroke_mrs, method = "spearman")

# SVD prevalence 
cat("WMH (Fazekas >= 1):", sum(data$fazekas >= 1, na.rm = TRUE), "/", sum(!is.na(data$fazekas)), "\n")
cat("Lacunes:", sum(data$lacunes > 0, na.rm = TRUE), "/", sum(!is.na(data$lacunes)), "\n")
cat("Microbleeds:", sum(data$cmb > 0, na.rm = TRUE), "/", sum(!is.na(data$cmb)), "\n")



#  3. Primary interaction models (Table 2) 

# all four CR x lesion interactions: 2 proxies x 2 outcomes
# mRS models adjust for age, sex, prestroke mRS
# NIHSS models adjust for age, sex (no prestroke mRS)

m_edu_mrs <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
                   prestroke_mrs + yoe_c + yoe_c:log_lesion_c, 
                 family = binomial, data = data_mrs)

m_work_mrs <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
                    prestroke_mrs + criq_work_10 + criq_work_10:log_lesion_c, 
                  family = binomial, data = data_mrs)

m_edu_nihss <- glm(nihss_1 ~ age_decade + log_lesion_c + biol_sex + 
                     yoe_c + yoe_c:log_lesion_c, 
                   family = Gamma(link = "log"), data = data_nihss)

m_work_nihss <- glm(nihss_1 ~ age_decade + log_lesion_c + biol_sex + 
                      criq_work_10_c + criq_work_10_c:log_lesion_c, 
                    family = Gamma(link = "log"), data = data_nihss)

# results
for (m in list(m_edu_mrs, m_work_mrs, m_edu_nihss, m_work_nihss)) {
  cat("\n")
  print(round(summary(m)$coefficients, 4))
  print(round(exp(cbind(OR = coef(m), confint(m))), 3))
}

# FDR across the 4 interaction tests
p_raw <- c(
  summary(m_edu_mrs)$coefficients["yoe_c:log_lesion_c", 4],
  summary(m_work_mrs)$coefficients["log_lesion_c:criq_work_10", 4],
  summary(m_edu_nihss)$coefficients["yoe_c:log_lesion_c", 4],
  summary(m_work_nihss)$coefficients["criq_work_10_c:log_lesion_c", 4]
)
cat("\nRaw p:", round(p_raw, 4), "\n")
cat("FDR p:", round(p.adjust(p_raw, method = "BH"), 4), "\n")



# 4. Stratified analyses (Tables A.4, A.5) 

# within each tertile: CR proxy effects without interaction term
# mRS models: age, sex, prestroke mRS, CR proxy (no lesion volume)
# NIHSS models: age, sex, CR proxy

cat("\n--- mRS by tertile ---\n")
for (q in c("Small", "Medium", "Large")) {
  d <- data_mrs %>% filter(lesion_t3 == q)
  cat("\n", q, "| n:", nrow(d), "| events:", sum(d$mrs_dichot), "\n")
  
  m_e <- glm(mrs_dichot ~ age_decade + biol_sex + prestroke_mrs + yoe_c, 
             family = binomial, data = d)
  cat("  Edu OR:", round(exp(coef(m_e)["yoe_c"]), 3),
      round(exp(confint(m_e, "yoe_c")), 3),
      "p:", round(summary(m_e)$coefficients["yoe_c", 4], 4), "\n")
  
  d$criq_work_10 <- d$criq_work / 10
  m_w <- glm(mrs_dichot ~ age_decade + biol_sex + prestroke_mrs + criq_work_10, 
             family = binomial, data = d)
  cat("  Work/10 OR:", round(exp(coef(m_w)["criq_work_10"]), 3),
      round(exp(confint(m_w, "criq_work_10")), 3),
      "p:", round(summary(m_w)$coefficients["criq_work_10", 4], 4), "\n")
}

cat("\n--- NIHSS by tertile ---\n")
for (q in c("Small", "Medium", "Large")) {
  d <- data_nihss %>% filter(lesion_t3 == q)
  cat("\n", q, "| n:", nrow(d), "\n")
  
  m_e <- glm(nihss_1 ~ age_decade + biol_sex + yoe, 
             family = Gamma(link = "log"), data = d)
  cat("  Edu MR:", round(exp(coef(m_e)["yoe"]), 3),
      "p:", round(summary(m_e)$coefficients["yoe", 4], 4), "\n")
  
  d$criq_work_10 <- d$criq_work / 10
  m_w <- glm(nihss_1 ~ age_decade + biol_sex + criq_work_10, 
             family = Gamma(link = "log"), data = d)
  cat("  Work/10 MR:", round(exp(coef(m_w)["criq_work_10"]), 3),
      "p:", round(summary(m_w)$coefficients["criq_work_10", 4], 4), "\n")
}


# ordinal mRS within each tertile
cat("\n--- Ordinal mRS by tertile ---\n")
for (q in c("Small", "Medium", "Large")) {
  d <- data_mrs %>% filter(lesion_t3 == q)
  m <- polr(factor(mrs_fu) ~ age_decade + biol_sex + prestroke_mrs + yoe_c, 
            data = d, Hess = TRUE)
  ct <- coef(summary(m))
  p <- pnorm(abs(ct["yoe_c", "t value"]), lower.tail = FALSE) * 2
  cat(q, "| OR:", round(exp(ct["yoe_c", "Value"]), 3), "| p:", round(p, 4), "\n")
}


# 5.  Sensitivity analyses 

# Excluding extreme lesion volumes (> 3 SD above mean)
upper <- mean(data_mrs$Volume_roi_cL) + 3 * sd(data_mrs$Volume_roi_cL)
data_trim <- data_mrs %>% filter(Volume_roi_cL <= upper)
cat("\n--- 3 SD exclusion: dropped", nrow(data_mrs) - nrow(data_trim), "---\n")

for (q in c("Small", "Medium", "Large")) {
  d <- data_trim %>% filter(lesion_t3 == q)
  m <- glm(mrs_dichot ~ age_decade + biol_sex + prestroke_mrs + yoe_c, 
           family = binomial, data = d)
  cat(q, "n:", nrow(d), "events:", sum(d$mrs_dichot),
      "OR:", round(exp(coef(m)["yoe_c"]), 3),
      "p:", round(summary(m)$coefficients["yoe_c", 4], 4), "\n")
}

#  Median split
med <- median(data_mrs$Volume_roi_cL)
data_mrs$les_med <- ifelse(data_mrs$Volume_roi_cL <= med, "Smaller", "Larger")
cat("\n--- Median split (median:", round(med, 2), "cL) ---\n")

for (q in c("Smaller", "Larger")) {
  d <- data_mrs %>% filter(les_med == q)
  m_e <- glm(mrs_dichot ~ age_decade + biol_sex + prestroke_mrs + yoe_c, 
             family = binomial, data = d)
  d$criq_work_10 <- d$criq_work / 10
  m_w <- glm(mrs_dichot ~ age_decade + biol_sex + prestroke_mrs + criq_work_10, 
             family = binomial, data = d)
  cat(q, "n:", nrow(d), "events:", sum(d$mrs_dichot), "\n")
  cat("  Edu:", round(exp(coef(m_e)["yoe_c"]), 3), "p:", round(summary(m_e)$coefficients["yoe_c", 4], 4), "\n")
  cat("  Work/10:", round(exp(coef(m_w)["criq_work_10"]), 3), "p:", round(summary(m_w)$coefficients["criq_work_10", 4], 4), "\n")
}

# Alternative mRS thresholds + ordinal
cat("\n--- Alternative thresholds ---\n")

data_mrs$mrs_12 <- ifelse(data_mrs$mrs_fu >= 2, 1, 0)
m_12 <- glm(mrs_12 ~ age_decade + log_lesion_c + biol_sex + prestroke_mrs + 
              yoe_c + yoe_c:log_lesion_c, family = binomial, data = data_mrs)
cat("0-1 vs 2-6: events:", sum(data_mrs$mrs_12),
    "p:", round(summary(m_12)$coefficients["yoe_c:log_lesion_c", 4], 4), "\n")

data_mrs$mrs_34 <- ifelse(data_mrs$mrs_fu >= 4, 1, 0)
m_34 <- glm(mrs_34 ~ age_decade + log_lesion_c + biol_sex + prestroke_mrs + 
              yoe_c + yoe_c:log_lesion_c, family = binomial, data = data_mrs)
cat("0-3 vs 4-6: events:", sum(data_mrs$mrs_34),
    "p:", round(summary(m_34)$coefficients["yoe_c:log_lesion_c", 4], 4), "\n")

m_ord <- polr(factor(mrs_fu) ~ age_decade + log_lesion_c + biol_sex + 
                prestroke_mrs + yoe_c + yoe_c:log_lesion_c, 
              data = data_mrs, Hess = TRUE)
ct <- coef(summary(m_ord))
p_ord <- pnorm(abs(ct["yoe_c:log_lesion_c", "t value"]), lower.tail = FALSE) * 2
cat("Ordinal: p:", round(p_ord, 4), "\n")

# Restricting to prestroke mRS = 0
d_pre0 <- data_mrs %>% filter(prestroke_mrs == 0)
m_pre0 <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
                yoe_c + yoe_c:log_lesion_c, 
              family = binomial, data = d_pre0)
cat("\n--- Prestroke mRS = 0 (n:", nrow(d_pre0), ") ---\n")
cat("Interaction p:", round(summary(m_pre0)$coefficients["yoe_c:log_lesion_c", 4], 4), "\n")
print(round(exp(cbind(OR = coef(m_pre0)["yoe_c:log_lesion_c"], 
                       confint(m_pre0, "yoe_c:log_lesion_c"))), 3))

# Additional covariates, each added separately to the primary model
cat("\n--- Covariate sensitivity (Table A.11) ---\n")

m_hemi <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
                prestroke_mrs + yoe_c + yoe_c:log_lesion_c + lesion_location, 
              family = binomial, data = data_mrs)
cat("+ Hemisphere:", round(summary(m_hemi)$coefficients["yoe_c:log_lesion_c", 4], 4), "\n")

m_cohort <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
                  prestroke_mrs + yoe_c + yoe_c:log_lesion_c + sample, 
                family = binomial, data = data_mrs)
cat("+ Cohort:", round(summary(m_cohort)$coefficients["yoe_c:log_lesion_c", 4], 4), "\n")

m_therapy <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
                   prestroke_mrs + yoe_c + yoe_c:log_lesion_c + therapy, 
                 family = binomial, data = data_mrs)
cat("+ Therapy:", round(summary(m_therapy)$coefficients["yoe_c:log_lesion_c", 4], 4), "\n")


# NIHSS as a potential mediator sensitivity 
m_nihss_adj <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
                     prestroke_mrs + yoe_c + yoe_c:log_lesion_c + nihss_24h, 
                   family = binomial, data = data_mrs)
cat("+ NIHSS 24h:", round(summary(m_nihss_adj)$coefficients["yoe_c:log_lesion_c", 4], 4), "\n")

# SVD marker adjustment
data_faz <- data_mrs %>% filter(!is.na(fazekas))
data_faz$yoe_c <- data_faz$yoe - mean(data_faz$yoe)
data_faz$log_lesion_c <- data_faz$log_lesion - mean(data_faz$log_lesion)
m_faz <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
               prestroke_mrs + yoe_c + yoe_c:log_lesion_c + fazekas, 
             family = binomial, data = data_faz)
cat("+ Fazekas (n:", nrow(data_faz), "):", round(summary(m_faz)$coefficients["yoe_c:log_lesion_c", 4], 4), "\n")

data_lac <- data_mrs %>% filter(!is.na(lacunes))
data_lac$yoe_c <- data_lac$yoe - mean(data_lac$yoe)
data_lac$log_lesion_c <- data_lac$log_lesion - mean(data_lac$log_lesion)
m_lac <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
               prestroke_mrs + yoe_c + yoe_c:log_lesion_c + I(lacunes > 0), 
             family = binomial, data = data_lac)
cat("+ Lacunes (n:", nrow(data_lac), "):", round(summary(m_lac)$coefficients["yoe_c:log_lesion_c", 4], 4), "\n")

data_cmb <- data_mrs %>% filter(!is.na(cmb))
data_cmb$yoe_c <- data_cmb$yoe - mean(data_cmb$yoe)
data_cmb$log_lesion_c <- data_cmb$log_lesion - mean(data_cmb$log_lesion)
m_cmb <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
               prestroke_mrs + yoe_c + yoe_c:log_lesion_c + I(cmb > 0), 
             family = binomial, data = data_cmb)
cat("+ Microbleeds (n:", nrow(data_cmb), "):", round(summary(m_cmb)$coefficients["yoe_c:log_lesion_c", 4], 4), "\n")


# Non-linearity checks

cat("\n--- Non-linearity ---\n")


# quadratic education x lesion interaction
m_quad <- glm(mrs_dichot ~ age_decade + log_lesion_c + I(log_lesion_c^2) + 
                biol_sex + prestroke_mrs + yoe_c + 
                yoe_c:log_lesion_c + yoe_c:I(log_lesion_c^2), 
              family = binomial, data = data_mrs)
cat("Quadratic interaction:", round(summary(m_quad)$coefficients["yoe_c:I(log_lesion_c^2)", 4], 4), "\n")

# quadratic education main effect
m_quad_edu <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
                    prestroke_mrs + yoe_c + I(yoe_c^2) + yoe_c:log_lesion_c, 
                  family = binomial, data = data_mrs)
cat("Quadratic education:", round(summary(m_quad_edu)$coefficients["I(yoe_c^2)", 4], 4), "\n")




# Firth penalised regression in the small tertile

d_small <- data_mrs %>% filter(lesion_t3 == "Small")
m_firth <- logistf(mrs_dichot ~ age_decade + biol_sex + prestroke_mrs + yoe, 
                   data = d_small)
cat("\n--- Firth (small tertile, n:", nrow(d_small), ") ---\n")
cat("OR:", round(exp(m_firth$coefficients["yoe"]), 3),
    "CI:", round(exp(m_firth$ci.lower["yoe"]), 3), "-", round(exp(m_firth$ci.upper["yoe"]), 3),
    "p:", round(m_firth$prob["yoe"], 4), "\n")

# Adding continuous lesion volume back into the small tertile model

d_small$yoe_c <- d_small$yoe - mean(data_mrs$yoe)
d_small$log_lesion_c <- d_small$log_lesion - mean(data_mrs$log_lesion)
m_small_adj <- glm(mrs_dichot ~ age_decade + log_lesion_c + biol_sex + 
                     prestroke_mrs + yoe_c, 
                   family = binomial, data = d_small)
cat("\n--- Small tertile + continuous lesion ---\n")
cat("OR:", round(exp(coef(m_small_adj)["yoe_c"]), 3),
    round(exp(confint(m_small_adj, "yoe_c")), 3),
    "p:", round(summary(m_small_adj)$coefficients["yoe_c", 4], 4), "\n")

#Education per SD 
sd_yoe <- sd(data_mrs$yoe)
beta <- coef(glm(mrs_dichot ~ age_decade + biol_sex + prestroke_mrs + yoe, 
                 family = binomial, 
                 data = data_mrs %>% filter(lesion_t3 == "Small")))["yoe"]
cat("\n--- Per SD (SD:", round(sd_yoe, 2), "years) ---\n")
cat("OR per SD:", round(exp(beta * sd_yoe), 3), "\n")



#  6. Exploratory prediction (Table A.12) 

# does adding CR to standard predictors improve discrimination?

cat("\n--- Full sample ---\n")
m_base <- glm(mrs_dichot ~ age_decade + log_lesion + biol_sex + prestroke_mrs, 
              family = binomial, data = data_mrs)
m_edu <- glm(mrs_dichot ~ age_decade + log_lesion + biol_sex + prestroke_mrs + yoe, 
             family = binomial, data = data_mrs)
cat("Base AUC:", round(auc(roc(data_mrs$mrs_dichot, fitted(m_base), quiet = TRUE)), 3), "\n")
cat("+ Edu AUC:", round(auc(roc(data_mrs$mrs_dichot, fitted(m_edu), quiet = TRUE)), 3),
    "LRT p:", round(anova(m_base, m_edu, test = "Chisq")$`Pr(>Chi)`[2], 4), "\n")

cat("\n--- Small tertile ---\n")
m_base_s <- glm(mrs_dichot ~ age_decade + biol_sex + prestroke_mrs, 
                family = binomial, data = d_small)
m_edu_s <- glm(mrs_dichot ~ age_decade + biol_sex + prestroke_mrs + yoe, 
               family = binomial, data = d_small)
d_small$criq_work_10 <- d_small$criq_work / 10
m_work_s <- glm(mrs_dichot ~ age_decade + biol_sex + prestroke_mrs + criq_work_10, 
                family = binomial, data = d_small)
cat("Base AUC:", round(auc(roc(d_small$mrs_dichot, fitted(m_base_s), quiet = TRUE)), 3), "\n")
cat("+ Edu AUC:", round(auc(roc(d_small$mrs_dichot, fitted(m_edu_s), quiet = TRUE)), 3),
    "LRT:", round(anova(m_base_s, m_edu_s, test = "Chisq")$`Pr(>Chi)`[2], 4), "\n")
cat("+ Work AUC:", round(auc(roc(d_small$mrs_dichot, fitted(m_work_s), quiet = TRUE)), 3),
    "LRT:", round(anova(m_base_s, m_work_s, test = "Chisq")$`Pr(>Chi)`[2], 4), "\n")



# 7. Figure 1: marginal effect of education across lesion volumes 



b <- coef(m_edu_mrs)
V <- vcov(m_edu_mrs)

les_range <- seq(min(data_mrs$log_lesion_c), 
                 quantile(data_mrs$log_lesion_c, 0.95), 
                 length.out = 100)

or_edu <- exp(b["yoe_c"] + b["yoe_c:log_lesion_c"] * les_range)

se <- sqrt(V["yoe_c", "yoe_c"] + 
             les_range^2 * V["yoe_c:log_lesion_c", "yoe_c:log_lesion_c"] + 
             2 * les_range * V["yoe_c", "yoe_c:log_lesion_c"])
or_lo <- exp((b["yoe_c"] + b["yoe_c:log_lesion_c"] * les_range) - 1.96 * se)
or_hi <- exp((b["yoe_c"] + b["yoe_c:log_lesion_c"] * les_range) + 1.96 * se)

# back-transform x-axis to cL
les_cL <- exp(les_range + mean(data_mrs$log_lesion)) - 1

fig1_df <- data.frame(lesion_cL = les_cL, OR = or_edu, lo = or_lo, hi = or_hi)

fig1 <- ggplot(fig1_df, aes(x = lesion_cL, y = OR)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "steelblue") +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = c(0, 1, 2, 5, 10, 20)) +
  coord_cartesian(ylim = c(0.4, 2.5)) +
  labs(x = "Lesion Volume (cL)", y = "OR for Education (per year)") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

ggsave("Figure1_marginal_effects.tiff", fig1, 
       width = 8, height = 5, dpi = 300, compression = "lzw")



# 8. Figure 2: predicted probabilities by tertile 


m_fig2 <- glm(mrs_dichot ~ age_decade + lesion_t3 + biol_sex + 
                prestroke_mrs + yoe + yoe:lesion_t3, 
              family = binomial, data = data_mrs)

yoe_lo <- quantile(data_mrs$yoe, 0.05)
yoe_hi <- quantile(data_mrs$yoe, 0.95)
pred2 <- ggpredict(m_fig2, terms = c(paste0("yoe [", yoe_lo, ":", yoe_hi, " by=0.5]"), 
                                      "lesion_t3"))

fig2 <- ggplot(pred2, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.40)) +
  scale_color_manual(values = c("Small" = "#2A9D8F", "Medium" = "#457B9D", "Large" = "#E63946"),
                     labels = c("Smaller", "Intermediate", "Larger")) +
  scale_fill_manual(values = c("Small" = "#2A9D8F", "Medium" = "#457B9D", "Large" = "#E63946"),
                    labels = c("Smaller", "Intermediate", "Larger")) +
  labs(x = "Years of Education", y = "Predicted P(mRS 3-6)",
       color = "Lesion Volume Tertile", fill = "Lesion Volume Tertile") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("Figure2_education_by_tertile.tiff", fig2, 
       width = 8, height = 5.5, dpi = 300, compression = "lzw")


# 9. Model diagnostics 

cat("\n--- Logistic (education x mRS) ---\n")

hl <- hoslem.test(m_edu_mrs$y, fitted(m_edu_mrs), g = 10)
cat("Hosmer-Lemeshow: chi-sq", round(hl$statistic, 2), "p", round(hl$p.value, 4), "\n")

cat("VIF:\n")
print(round(vif(m_edu_mrs), 2))

cooks_d <- cooks.distance(m_edu_mrs)
cat("Cook's D max:", round(max(cooks_d, na.rm = TRUE), 4),
    "| above 4/n:", sum(cooks_d > 4/nobs(m_edu_mrs), na.rm = TRUE), "\n")

cat("AUC:", round(auc(roc(m_edu_mrs$y, fitted(m_edu_mrs), quiet = TRUE)), 3), "\n")

cat("\n--- Gamma GLM (education x NIHSS) ---\n")
cat("Dispersion:", round(summary(m_edu_nihss)$dispersion, 4), "\n")



# 10. Missing data comparison (Table A.3) 

data$included <- data$record_id %in% data_mrs$record_id
inc <- data %>% filter(included)
exc <- data %>% filter(!included)

cat("\nIncluded:", nrow(inc), "| Excluded:", nrow(exc), "\n")
for (v in c("age", "yoe", "criq_work", "nihss_24h", "Volume_roi_cL")) {
  p <- wilcox.test(inc[[v]], exc[[v]])$p.value
  cat(v, "— included:", round(median(inc[[v]], na.rm = TRUE), 2),
      "excluded:", round(median(exc[[v]], na.rm = TRUE), 2),
      "p:", round(p, 3), "\n")
}



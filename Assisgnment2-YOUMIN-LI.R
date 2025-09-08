############################################################
## AEB 6933 – Problem Set 2 (Kremer & Miguel 2004 “Economtrica”)
## Replication + Extensions
## Author(s): <Youmin>
############################################################
## ---- 0. Setup library
library(tidyverse)
library(haven)      
library(fixest)     
library(sandwich)   # vcovCL for clustered/robust vcovs
library(lmtest)     # coeftest for printing with custom vcov
library(margins)    # marginal effects
install.packages("remotes") 
remotes::install_github("grantmcdermott/ritest")
library(ritest)    # randomization inference (McDermott)
library(sandwich)
library(lmtest)
library(margins)
li
theme_set(theme_minimal())

## Load data 
dpath <- "ps2.dta"
df <- read_dta(dpath) %>%
  janitor::clean_names()   
names(df)
df <- as.data.frame(df)

## ---- 0.2 Standardize key variables-
# Dependent variable (any moderate/heavy infection in 1999)
y <- "any_ics99"
stopifnot(y %in% names(df))
# Weight
wvar <- "indiv_weight"
stopifnot(wvar %in% names(df))
treat <- "wgrp1"
# Cluster id (school)
clustvar <- "sch98v1"
stopifnot(clustvar %in% names(df))
# Group 1 (treated in 1998, treated before 1999)
pre_treat <- c("t1","treat98","treat_1998","group1","g1","t98","t1i")
treat  <- pre_treat[pre_treat %in% names(df)][1]
# Convention in prompt: pop1_3km_original (treated pupils within 0–3km),
pop1_03 <- c("pop1_0_3km_original","pop1_03km_original","pop1_3km_original")
pop1_36 <- c("pop1_3_6km_original","pop1_36km_original")
popT_03 <- c("popt_0_3km_original","popt_03km_original","popt_3km_original","pop_0_3km_original","pop_03km_original")
popT_36 <- c("popt_3_6km_original","popt_36km_original")

v_pop1_03 <- pop1_03[pop1_03 %in% names(df)][1]
v_pop1_36 <- pop1_36[pop1_36 %in% names(df)][1]
v_popT_03 <- popT_03[popT_03 %in% names(df)][1]
v_popT_36 <- popT_36[popT_36 %in% names(df)][1]
# pop1_0_1km_updated, popT_0_1km_updated, ... for 0–1, 1–2, …, 5–6 km
mk_name <- function(prefix, lo, hi) sprintf("%s_%d_%dkm_updated", prefix, lo, hi)
dist_bins <- list(c(0,1), c(1,2), c(2,3), c(3,4), c(4,5), c(5,6))
pop1_upd <- setNames(
  vapply(dist_bins, \(b) mk_name("pop1", b[1], b[2]), character(1)),
  vapply(dist_bins, \(b) paste0(b[1], "-", b[2], "km"), character(1))
)
popT_upd <- setNames(
  vapply(dist_bins, \(b) mk_name("popT", b[1], b[2]), character(1)),
  names(pop1_upd)
)
pop1_upd <- pop1_upd[pop1_upd %in% names(df)]
popT_upd <- popT_upd[popT_upd %in% names(df)]

# Covariates X: sap1–sap4, Istd4–Istd9, mk96_s
covars <- c(paste0("sap", 1:4), paste0("istd", 4:9), "mk96_s")
covars <- covars[covars %in% names(df)]
ctest <- function(model, cluster_formula) {
  vc <- sandwich::vcovCL(model, cluster = cluster_formula, type = "HC1")
  lmtest::coeftest(model, vcov. = vc)
}
vcov = ~ sch98v1
df$any_ics99 <- as.integer(df$any_ics99 > 0)

### Q1  - weight linear prob model (pre-1999 treatment for G1; G2 untreated in 1999)
f_q1 <- as.formula(paste0(y, " ~ ", treat))
m_q1 <- feols(f_q1, data=df, weights=~indiv_weight, cluster=~sch98v1)
summary(m_q1)
fixest::etable(m_q1)      

wmean <- function(x, w){
  ok <- is.finite(x) & is.finite(w)
  weighted.mean(x[ok], w[ok])
}
g0 <- df$wgrp1 == 0 & is.finite(df$any_ics99) & is.finite(df$indiv_weight)
g1 <- df$wgrp1 == 1 & is.finite(df$any_ics99) & is.finite(df$indiv_weight)
mean_ctrl <- weighted.mean(df$any_ics99[g0], df$indiv_weight[g0])
mean_t1   <- weighted.mean(df$any_ics99[g1], df$indiv_weight[g1])
c(mean_ctrl = mean_ctrl, mean_t1 = mean_t1, diff = mean_t1 - mean_ctrl)


### Q2 - Randomization inference (RI) for β(wgrp1) in Q1
set.seed(1234)
obs_beta <- coef(m_q1)["wgrp1"]
B <- 2000  
schools <- df |>
  dplyr::distinct(sch98v1, wgrp1) |>
  dplyr::arrange(sch98v1)
S <- nrow(schools)
S_treated <- sum(schools$wgrp1 == 1, na.rm = TRUE)
ri_stat <- numeric(B)
for (b in seq_len(B)) {
  tr_sch <- rep(0L, S)
  tr_sch[sample.int(S, S_treated)] <- 1L  # reassign treated schools
  perm_map <- data.frame(sch98v1 = schools$sch98v1, T_perm = tr_sch)
  
  df_perm <- df |>
    dplyr::left_join(perm_map, by = "sch98v1")
  
  mb <- feols(any_ics99 ~ T_perm, data = df_perm,
              weights = ~ indiv_weight, cluster = ~ sch98v1)
  ri_stat[b] <- coef(mb)["T_perm"]
}

# two-sided RI p-value
p_ri <- mean(abs(ri_stat) >= abs(obs_beta))
p_ri
# quick visual
hist(ri_stat, breaks = 40, main = "RI null: β(T) under re-randomization",
     xlab = "β̂ under shuffled treatment (school-level)")
abline(v = obs_beta, lwd = 2)

## Q5: Probit with externality counts and totals, clustered SEs 
rhs <- c("wgrp1", covars, v_pop1_03, v_pop1_36, v_popT_03, v_popT_36)
rhs <- rhs[!is.na(rhs) & rhs %in% names(df)]
f_q5 <- as.formula(paste("any_ics99 ~", paste(rhs, collapse = " + ")))
m_q5 <- glm(f_q5, data = df,
            family = binomial(link = "probit"),
            weights = df$indiv_weight)

# cluster school-level SEs
vc_q5 <- sandwich::vcovCL(m_q5, cluster = ~ sch98v1, type = "HC1")
q5_tab <- lmtest::coeftest(m_q5, vcov. = vc_q5)
q5_tab

# compute robust SE & CI directly
rob_se <- sqrt(diag(V_q5))
q5_out <- data.frame(
  term     = names(coef(m_q5)),
  estimate = unname(coef(m_q5)),
  std.error= unname(rob_se)
)
q5_out$z    <- q5_out$estimate / pmax(q5_out$std.error, .Machine$double.eps)
q5_out$p    <- 2*pnorm(-abs(q5_out$z))
q5_out$ci.lo<- q5_out$estimate - 1.96*q5_out$std.error
q5_out$ci.hi<- q5_out$estimate + 1.96*q5_out$std.error
print(q5_out)


## ---- Q6-7 margins() at means is the default.
# Cluster-robust variance-covariance matrix for glm/probit
vcCL <- function(model, cluster) {
  cluster_var <- model.frame(cluster, data = model$data)[[1]]
  sandwich::vcovCL(model, cluster = cluster_var)
}

V_q5 <- vcCL(m_q5, cluster = ~ sch98v1)
ct_q5 <- lmtest::coeftest(m_q5, vcov. = V_q5)
print(ct_q5)
marg_q7 <- margins::margins(
  m_q5,
  vcov = V_q5,  
  at = NULL       
)
sum_marg_q7 <- summary(marg_q7)
print(sum_marg_q7)

get_var_ame <- function(marg_obj, var, per = 1000) {
  out <- as.data.frame(summary(marg_obj))
  deriv_col <- if ("dydx" %in% names(out)) "dydx" else
    if ("AME"  %in% names(out)) "AME"  else
      stop("No derivative column named 'dydx' or 'AME' in summary(marg_obj).")
  se_col    <- if ("SE" %in% names(out)) "SE" else
    if ("Std. Error" %in% names(out)) "Std. Error" else
      stop("No SE column named 'SE' or 'Std. Error' in summary(marg_obj).")
  out <- dplyr::filter(out, .data$factor == var)
  ame_raw <- as.numeric(out[[deriv_col]])
  se_raw  <- as.numeric(out[[se_col]])
  
  out |>
    dplyr::mutate(
      AME_per = ame_raw * per,
      SE_per  = se_raw  * per,
      z       = AME_per / SE_per,
      p       = 2 * pnorm(-abs(z))
    ) |>
    dplyr::select(factor, AME_per, SE_per, z, p)
}

ame_wgrp1   <- get_var_ame(marg_q7, "wgrp1", per = 1)
ame_3km_1k  <- get_var_ame(marg_q7, "pop1_3km_original", per = 1000)
print(ame_wgrp1)
print(ame_3km_1k)

## Q8 quadratic treated-count terms: pop1_d + I(pop1_d^2)
rhs_q8 <- c(
  "wgrp1", covars,
  # treated counts and their squares (only add squares if the base var exists)
  if (!is.na(v_pop1_03)) c(v_pop1_03, sprintf("I(%s^2)", v_pop1_03)) else NULL,
  if (!is.na(v_pop1_36)) c(v_pop1_36, sprintf("I(%s^2)", v_pop1_36)) else NULL,
  # totals remain linear
  v_popT_03, v_popT_36
)
rhs_q8 <- rhs_q8[!is.na(rhs_q8)]
f_q8 <- as.formula(paste("any_ics99 ~", paste(rhs_q8, collapse = " + ")))

m_q8 <- glm(f_q8, data = df, family = binomial(link = "probit"), weights = indiv_weight)
V_q8 <- vcCL(m_q8, ~ sch98v1)
marg_q8 <- margins::margins(m_q8, vcov = V_q8)
print(marg_q8)

# -------------------------------
# Q9) Marginal effects for quadratic model and comparison to Q7
if (!is.na(v_pop1_03)) ame_q8_0_3km_1k <- get_var_ame(marg_q8, v_pop1_03, per = 1000)
if (!is.na(v_pop1_36)) ame_q8_3_6km_1k <- get_var_ame(marg_q8, v_pop1_36, per = 1000)
if (exists("ame_q8_0_3km_1k")) print(ame_q8_0_3km_1k)
if (exists("ame_q8_3_6km_1k")) print(ame_q8_3_6km_1k)
marg_q9 <- margins::margins(m_q8, vcov = V_q8)
sum_marg_q9 <- summary(marg_q9)
print(sum_marg_q9)

# -------------------------------
# Q10) Fine-distance bands 0–1,1–2,2–3,3–4,4–5,5–6 km (updated variables)
bands <- c("0_1","1_2","2_3","3_4","4_5","5_6")
guess_band_vars <- function(prefix=c("pop1","popT"), band, suffix="_km_updated") {
  c(
    paste0(prefix, "_", band, suffix),                 # pop1_0_1_km_updated (if someone typed underscore before km)
    paste0(prefix, "_", band, "_km_updated"),          # pop1_0_1_km_updated
    paste0(prefix, "_", gsub("_","",band), "km_updated"), # pop1_01km_updated
    paste0(prefix, "_", band, "km_updated"),           # pop1_0_1km_updated
    paste0(prefix, "_", band, "_km_updated"),          # pop1_0_1_km_updated (dup)
    paste0(prefix, "_", band, "_updated"),             # pop1_0_1km_updated (missing 'km')
    paste0(prefix, "_", band, "km")                    # fallback without "_updated"
  )
}
pick_one <- function(candidates, nms) {
  hit <- candidates[candidates %in% nms]  
  if (length(hit)) {
    hit[1]   # return the first valid match
  } else {
    NA_character_   
  }
}
pop1_bands <- unlist(lapply(bands, function(b) pick_one(guess_band_vars("pop1", b), names(df))))
popT_bands <- unlist(lapply(bands, function(b) pick_one(guess_band_vars("popT", b), names(df))))

idx1 <- which(!is.na(pop1_bands))
idxT <- which(!is.na(popT_bands))
pop1_bands <- pop1_bands[idx1]
popT_bands <- popT_bands[idxT]

rhs_q10 <- c("wgrp1", covars, pop1_bands, popT_bands)
f_q10 <- as.formula(paste("any_ics99 ~", paste(rhs_q10, collapse = " + ")))
m_q10 <- glm(f_q10, data = df, family = binomial(link = "probit"), weights = indiv_weight)
V_q10 <- vcCL(m_q10, ~ sch98v1)
ct_q10 <- lmtest::coeftest(m_q10, vcov. = V_q10)
print(ct_q10)
q10_table <- lmtest::coeftest(m_q10, vcov. = V_q10)
print(q10_table)
# -------------------------------
## Q11) Marginal externality effect for each distance band
marg_q11 <- margins::margins(m_q10, vcov = V_q10)
sum_marg_q11 <- summary(marg_q11)
deriv_col <- if ("dydx" %in% names(sum_marg_q11)) "dydx" else if ("AME" %in% names(sum_marg_q11)) "AME" else
  stop("Q11: Couldn't find derivative column 'dydx' or 'AME' in margins summary.")
se_col    <- if ("SE" %in% names(sum_marg_q11)) "SE" else if ("Std. Error" %in% names(sum_marg_q11)) "Std. Error" else
  stop("Q11: Couldn't find SE column 'SE' or 'Std. Error' in margins summary.")

ame_bands <- sum_marg_q11 |>
  dplyr::filter(.data$factor %in% pop1_bands) |>
  dplyr::transmute(
    band  = .data$factor,
    AME   = suppressWarnings(as.numeric(.data[[deriv_col]])) * 1000,
    SE    = suppressWarnings(as.numeric(.data[[se_col]]))    * 1000
  ) |>
  dplyr::mutate(
    z     = AME / pmax(SE, .Machine$double.eps),
    p     = 2 * pnorm(-abs(z)),
    ci_lo = AME - 1.96 * SE,
    ci_hi = AME + 1.96 * SE
  )

print(ame_bands)
## Q12) Plot the marginal externality effects across distance
band_order <- c("pop1_0_1km_updated","pop1_1_2km_updated","pop1_2_3km_updated",
                "pop1_3_4km_updated","pop1_4_5km_updated","pop1_5_6km_updated")
band_labels <- c("0–1 km","1–2 km","2–3 km","3–4 km","4–5 km","5–6 km")

plot_df <- ame_bands %>%
  dplyr::mutate(band = factor(band, levels = band_order, labels = band_labels))

p_q12 <- ggplot(plot_df, aes(x = band, y = AME)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = AME - 1.96*SE, ymax = AME + 1.96*SE), width = 0.15) +
  labs(
    title = "Marginal externality effect by distance band",
    subtitle = "AME of +1,000 treated students in band (probit, clustered SE by school)",
    x = "Distance from school",
    y = "Δ Pr(any_ics99=1) per +1,000 treated students"
  )
print(p_q12)


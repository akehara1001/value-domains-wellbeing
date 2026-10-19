# =============================================================================
# 07_reference_one_domain.R -- every breadth model re-expressed against
# ONE ENDORSED DOMAIN as the reference category (co-author comments #11, #13,
# #14).
#
# Why. The commitment-adjusted column already used one endorsed domain as its
# reference, because average commitment intensity is undefined for respondents
# who endorsed none. The other columns used no endorsed domain, so the two were
# not on a common footing and the comparison the manuscript drew between them
# did not hold. The 215 respondents who endorsed no domain are also a mixed
# group: people who value something outside the nine offered domains cannot be
# told apart from people who value nothing. They stay in the analysis as a
# level; they no longer serve as the baseline.
#
# This is a REPARAMETERISATION, not a different model. Fit, the differences
# between non-reference levels, and the saturation test must come out
# numerically identical to the ref = 0 fits. The script asserts all three and
# stops if any of them moves.
#
#   input   values_analysis.rds   (00_derive.R)
#           lca_models.rds        (01_lca.R)
#   output  outputs/table22_breadth_ref1.csv        coefficients, all models
#           outputs/table23_ref1_wald.csv           joint and saturation tests
#           outputs/table24_ref_crosscheck.csv      ref = 0 against ref = 1
#           outputs/table25_firstorder_ref1_std.csv standardised, ref = 1
#           outputs/table26_breadth_by_class.csv    Table 1 breadth rows
#           outputs/table27_endorsement_rates.csv   observed rate per domain
# =============================================================================
suppressPackageStartupMessages({ library(lavaan); library(dplyr) })
options(stringsAsFactors = FALSE)
BASE <- "."  # the repository root; run every script from there
setwd(BASE)

d  <- readRDS("values_analysis.rds")
LC <- readRDS("lca_models.rds")
K  <- 4L
m  <- LC$models[[K]]
CLASS_LAB <- c("Minimal", "Relational-occupational", "Family-health",
               "Expansive")

data <- as.data.frame(d)
data$vclass   <- factor(m$predclass, levels = seq_len(K))
data$age_z    <- as.numeric(scale(data$age))
data$income_z <- as.numeric(scale(data$income))
data$k6_z     <- as.numeric(scale(data$k6))
data$mcom_z   <- as.numeric(scale(data$m_commitment))
FCTS <- c("vclass", "disaster", "fcancer2", "female", "edu2", "occupation2",
          "marriage", "comobidity", "fukushima", "health2")
data <- data |> mutate(across(all_of(FCTS), factor))
data <- cbind(data, model.matrix(~ vclass + disaster + fcancer2 + female + edu2 +
                                   occupation2 + marriage + comobidity +
                                   fukushima + health2, data)[, -1])
stopifnot(sum(data$disaster1) == 95, sum(data$disaster2) == 372)

COV   <- paste("disaster1 + disaster2 + fcancer21 + female1 + age_z + income_z +",
               "edu22 + edu23 + occupation22 + occupation23 + marriage2 +",
               "marriage3 + comobidity1 + fukushima1 + health22")
COV_B <- sub("fcancer21 + ", "", COV, fixed = TRUE)
MEAS_BASE <- paste("NFE =~ sadR + worryR",
                   "ResAll =~ ls + gh + es + vital + comp",
                   "WB =~ NFE + ResAll", "ls ~~ gh", "es ~~ vital", sep = "\n")
MEAS      <- paste0(MEAS_BASE, "\nResAll ~~ 0*ResAll")
MEAS_FREE <- MEAS_BASE                      # the K6 model keeps it free

LVL <- c("0", "1", "2", "3", "4", "5", "6+")
data$nsupcat <- cut(pmin(data$nsup, 6), breaks = -1:6, labels = LVL)
for (lv in 0:6) data[[paste0("nsc", lv)]] <- as.numeric(pmin(data$nsup, 6) == lv)
nscat_n     <- table(data$nsupcat)
nscat_n_grp <- table(data$nsupcat, data$fcancer2)

# reference = 1 endorsed domain: every level except nsc1 enters the model
REF1 <- paste0("nsc", c(0, 2:6))
REF0 <- paste0("nsc", 1:6)

terms_of <- function(prefix, v)
  paste(sprintf("%s%s*%s", prefix, sub("nsc", "", v), v), collapse = " + ")

# ---------------------------------------------------------------- ref = 1 fits
M1    <- paste0(MEAS,      "\nWB ~ ", terms_of("a", REF1), " + ", COV, "\n")
M1_k6 <- paste0(MEAS_FREE, "\nWB ~ ", terms_of("k", REF1), " + ", COV,
                " + k6_z\n")
M1_g  <- paste0(MEAS, "\nWB ~ ",
                paste(sprintf("c(g0_%s, g1_%s)*%s", sub("nsc", "", REF1),
                              sub("nsc", "", REF1), REF1), collapse = " + "),
                " + ", COV_B, "\n")
fit1    <- sem(M1,    data, estimator = "MLR", missing = "ML")
fit1_k6 <- sem(M1_k6, data, estimator = "MLR", missing = "ML")
fit1_g  <- sem(M1_g,  data, group = "fcancer2", estimator = "MLR",
               missing = "ML")

# the commitment column is unchanged: it already used one endorsed domain, and
# it is fitted to the 1964 respondents for whom commitment intensity is defined
dsub <- data[!is.na(data$m_commitment), ]
stopifnot(nrow(dsub) == 1964, all(dsub$nsup >= 1))
M1_com <- paste0(MEAS, "\nWB ~ ", terms_of("c", paste0("nsc", 2:6)), " + ",
                 COV, " + mcom_z\n")
fit1_com <- sem(M1_com, dsub, estimator = "MLR", missing = "ML")

# ------------------------------------------------- the ref = 0 fits, for check
M0    <- paste0(MEAS,      "\nWB ~ ", terms_of("b", REF0), " + ", COV, "\n")
M0_k6 <- paste0(MEAS_FREE, "\nWB ~ ", terms_of("m", REF0), " + ", COV,
                " + k6_z\n")
M0_g  <- paste0(MEAS, "\nWB ~ ",
                paste(sprintf("c(h0_%s, h1_%s)*%s", sub("nsc", "", REF0),
                              sub("nsc", "", REF0), REF0), collapse = " + "),
                " + ", COV_B, "\n")
fit0    <- sem(M0,    data, estimator = "MLR", missing = "ML")
fit0_k6 <- sem(M0_k6, data, estimator = "MLR", missing = "ML")
fit0_g  <- sem(M0_g,  data, group = "fcancer2", estimator = "MLR",
               missing = "ML")

# -------------------------------------------------------------- extraction
lev_rows <- function(fit, group = NULL, keep) {
  pe <- parameterEstimates(fit, standardized = TRUE)
  pe <- pe[pe$op == "~" & pe$lhs == "WB" & pe$rhs %in% keep, ]
  if (!is.null(group) && "group" %in% names(pe)) pe <- pe[pe$group == group, ]
  pe$level <- LVL[as.integer(sub("nsc", "", pe$rhs)) + 1]
  pe[order(as.integer(sub("nsc", "", pe$rhs))), ]
}

pack <- function(fit, label, keep, group = NULL, n_src = nscat_n) {
  r <- lev_rows(fit, group = group, keep = keep)
  data.frame(model = label, level = r$level,
             n = as.integer(n_src[r$level]),
             b = round(r$est, 4), se = round(r$se, 4),
             ci_low = round(r$ci.lower, 4), ci_high = round(r$ci.upper, 4),
             beta = round(r$std.all, 3), p = signif(r$pvalue, 4))
}

pack_raw <- function(fit, label, keep, group = NULL) {
  r <- lev_rows(fit, group = group, keep = keep)
  data.frame(model = label, level = r$level, b = r$est,
             ci_low = r$ci.lower, ci_high = r$ci.upper, p = r$pvalue)
}

out <- rbind(
  pack(fit1,     "Primary model (ref = 1 domain)",            REF1),
  pack(fit1_com, "Adjusted for commitment intensity (ref = 1 domain)",
       paste0("nsc", 2:6), n_src = table(dsub$nsupcat)),
  pack(fit1_k6,  "Adjusted for psychological distress (ref = 1 domain)", REF1),
  pack(fit1_g,   "No family cancer history (ref = 1 domain)", REF1, group = 1,
       n_src = nscat_n_grp[, 1]),
  pack(fit1_g,   "Family cancer history (ref = 1 domain)",    REF1, group = 2,
       n_src = nscat_n_grp[, 2]))
write.csv(out, "outputs/table22_breadth_ref1.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

# The same numbers at full working precision. Rounding a four-decimal value a
# second time to three moves anything sitting exactly on a boundary, so the
# manuscript and its checker read the printed digits from here.
outp <- rbind(
  pack_raw(fit1,     "Primary model (ref = 1 domain)",            REF1),
  pack_raw(fit1_com, "Adjusted for commitment intensity (ref = 1 domain)",
           paste0("nsc", 2:6)),
  pack_raw(fit1_k6,  "Adjusted for psychological distress (ref = 1 domain)",
           REF1),
  pack_raw(fit1_g,   "No family cancer history (ref = 1 domain)", REF1,
           group = 1),
  pack_raw(fit1_g,   "Family cancer history (ref = 1 domain)",    REF1,
           group = 2))
write.csv(outp, "outputs/table22b_breadth_ref1_precise.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

# ------------------------------------------------------------------- tests
wald <- function(fit, cons) {
  w <- lavTestWald(fit, constraints = paste(cons, collapse = "\n"))
  c(w$stat, w$df, w$p.value)
}
joint1  <- wald(fit1,    sprintf("a%s == 0", c(0, 2:6)))
joint1k <- wald(fit1_k6, sprintf("k%s == 0", c(0, 2:6)))
sat1    <- wald(fit1,    c("a4 == a5", "a5 == a6"))
sat1k   <- wald(fit1_k6, c("k4 == k5", "k5 == k6"))
sat1c   <- wald(fit1_com, c("c4 == c5", "c5 == c6"))
joint0  <- wald(fit0,    sprintf("b%s == 0", 1:6))
sat0    <- wald(fit0,    c("b4 == b5", "b5 == b6"))
grp1    <- wald(fit1_g,  sprintf("g0_%s == g1_%s", c(0, 2:6), c(0, 2:6)))
grp0    <- wald(fit0_g,  sprintf("h0_%s == h1_%s", 1:6, 1:6))

tests <- data.frame(
  test = c("Primary model, all levels = 0",
           "Primary model, levels 4 = 5 = 6+ (saturation)",
           "Distress-adjusted model, all levels = 0",
           "Distress-adjusted model, levels 4 = 5 = 6+",
           "Commitment-adjusted model, levels 4 = 5 = 6+",
           "Multiple-group model, breadth coefficients equal across groups",
           "ref = 0 multiple-group equality test (for comparison)",
           "ref = 0 primary model, all levels = 0 (for comparison)",
           "ref = 0 primary model, levels 4 = 5 = 6+ (for comparison)"),
  statistic = round(c(joint1[1], sat1[1], joint1k[1], sat1k[1], sat1c[1],
                      grp1[1], grp0[1], joint0[1], sat0[1]), 3),
  df = c(joint1[2], sat1[2], joint1k[2], sat1k[2], sat1c[2], grp1[2],
         grp0[2], joint0[2], sat0[2]),
  p = signif(c(joint1[3], sat1[3], joint1k[3], sat1k[3], sat1c[3], grp1[3],
               grp0[3], joint0[3], sat0[3]), 4))
write.csv(tests, "outputs/table23_ref1_wald.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

# ------------------------------------------------------------- cross-checks
# A reparameterisation leaves the fit, the differences between non-reference
# levels, and the saturation test untouched. Assert all three.
# SRMR is left out of the comparison on purpose. It averages residuals over the
# observed variables in the model, and the two parameterisations put different
# dummy columns there: nsc1 under ref = 0, nsc0 under ref = 1. The matrix it
# summarises therefore is not the same matrix, while the likelihood, the
# chi-square, CFI, RMSEA and every estimate are. Both values are written to the
# cross-check file so the difference is on the record rather than hidden.
fitm <- function(f) fitMeasures(f, c("chisq", "df", "cfi.robust",
                                     "rmsea.robust"))
fitm_all <- function(f) fitMeasures(f, c("chisq", "df", "cfi.robust",
                                         "rmsea.robust", "srmr"))
loglik <- function(f) as.numeric(logLik(f))
# Two runs of the same model under different parameterisations stop at slightly
# different points of the same optimum, so the comparison is relative and its
# tolerance is far tighter than any reported digit: 1e-5 relative on a
# chi-square near 620 is a difference in the seventh significant figure.
CHK <- list()
chk <- function(a, b, what, tol) {
  diff <- max(abs(a - b))
  ok <- diff < tol
  CHK[[what]] <<- c(diff, tol, as.numeric(ok))
  cat(sprintf("%-54s %s (max abs %.2e, tolerance %.0e)\n", what,
              if (ok) "MATCH" else "DIFFERS", diff, tol))
  ok
}
# Fit statistics are compared to 1e-3 absolute, which on a chi-square near 590
# is the sixth significant figure and on CFI or RMSEA the fourth decimal.
# Coefficients are reported to four decimals and compared to 1e-5.
ok <- c(
  chk(loglik(fit0),  loglik(fit1),  "primary model log-likelihood", 1e-4),
  chk(fitm(fit0),    fitm(fit1),    "primary model fit", 1e-3),
  chk(fitm(fit0_k6), fitm(fit1_k6), "distress-adjusted model fit", 1e-3),
  chk(fitm(fit0_g),  fitm(fit1_g),  "multiple-group model fit", 1e-3),
  chk(sat0[c(1, 3)], sat1[c(1, 3)], "saturation test, statistic and p", 1e-3),
  chk(joint0[2],     joint1[2],     "joint test degrees of freedom", 1e-9),
  chk(joint0[c(1, 3)], joint1[c(1, 3)], "joint test, statistic and p", 1e-3),
  chk(grp0[c(1, 3)], grp1[c(1, 3)],
      "group equality test, statistic and p", 1e-3))

# differences between non-reference levels, computed both ways
est0 <- setNames(lev_rows(fit0, keep = REF0)$est, 1:6)
est1 <- setNames(lev_rows(fit1, keep = REF1)$est, c(0, 2:6))
pairs_ <- t(combn(2:6, 2))
d0 <- est0[as.character(pairs_[, 2])] - est0[as.character(pairs_[, 1])]
d1 <- est1[as.character(pairs_[, 2])] - est1[as.character(pairs_[, 1])]
ok <- c(ok, chk(d0, d1, "differences between non-reference levels", 1e-5))
# the shift itself: every ref = 1 coefficient is the ref = 0 one minus b1
ok <- c(ok, chk(unname(est1[as.character(2:6)]),
                unname(est0[as.character(2:6)] - est0["1"]),
                "ref = 1 coefficients equal ref = 0 minus level 1", 1e-5))
ok <- c(ok, chk(unname(est1["0"]), unname(-est0["1"]),
                "level 0 coefficient equals minus level 1 under ref = 0",
                1e-5))
stopifnot(all(ok))

cross <- data.frame(
  quantity = c("chi-square", "degrees of freedom", "robust CFI",
               "robust RMSEA", "SRMR", "joint test statistic",
               "joint test df", "saturation statistic", "saturation df",
               "saturation p"),
  ref_0_domains = signif(c(fitm_all(fit0), joint0[1], joint0[2], sat0[1],
                           sat0[2], sat0[3]), 6),
  ref_1_domain  = signif(c(fitm_all(fit1), joint1[1], joint1[2], sat1[1],
                           sat1[2], sat1[3]), 6))
cross$abs_difference <- signif(abs(cross$ref_0_domains - cross$ref_1_domain), 3)
cross$agrees_to_reported_precision <-
  ifelse(abs(cross$ref_0_domains - cross$ref_1_domain) < 1e-3, "yes", "no")
write.csv(cross, "outputs/table24_ref_crosscheck.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

chk_out <- data.frame(quantity = names(CHK),
                      max_abs_difference = signif(sapply(CHK, `[`, 1), 3),
                      tolerance = sapply(CHK, `[`, 2),
                      matches = ifelse(sapply(CHK, `[`, 3) == 1, "yes", "no"))
write.csv(chk_out, "outputs/table24b_ref_reparameterisation_checks.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------------------- first-order factors, reference = 1
MFO1 <- paste("NFE =~ sadR + worryR",
              "ResAll =~ ls + gh + es + vital + comp",
              "ls ~~ gh", "es ~~ vital",
              paste0("NFE ~ ", terms_of("fn", REF1), " + ", COV),
              paste0("ResAll ~ ", terms_of("fr", REF1), " + ", COV),
              "NFE ~~ ResAll", sep = "\n")
fitFO1 <- sem(MFO1, data, estimator = "MLR", missing = "ML")
ss <- standardizedSolution(fitFO1)
ss <- ss[ss$op == "~" & ss$lhs %in% c("NFE", "ResAll") & ss$rhs %in% REF1, ]
pe <- parameterEstimates(fitFO1)
pe <- pe[pe$op == "~" & pe$lhs %in% c("NFE", "ResAll") & pe$rhs %in% REF1, ]
stopifnot(identical(paste(ss$lhs, ss$rhs), paste(pe$lhs, pe$rhs)))
fo1 <- data.frame(
  factor = ifelse(ss$lhs == "NFE", "NFE (negative-free emotion)",
                  "ResAll (resource-oriented well-being)"),
  level = LVL[as.integer(sub("nsc", "", ss$rhs)) + 1],
  b = round(pe$est, 4), b_ci_low = round(pe$ci.lower, 4),
  b_ci_high = round(pe$ci.upper, 4),
  beta = round(ss$est.std, 4), ci_low = round(ss$ci.lower, 4),
  ci_high = round(ss$ci.upper, 4), p = signif(ss$pvalue, 4))
write.csv(fo1, "outputs/table25_firstorder_ref1_std.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
fofit <- fitMeasures(fitFO1, c("chisq", "df", "cfi.robust", "rmsea.robust"))
write.csv(data.frame(quantity = names(fofit),
                     value = signif(as.numeric(fofit), 6)),
          "outputs/table25b_firstorder_ref1_fit.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
fo_n <- wald(fitFO1, sprintf("fn%s == 0", c(0, 2:6)))
fo_r <- wald(fitFO1, sprintf("fr%s == 0", c(0, 2:6)))

# The tests the manuscript reports that the 2026-08-30 inventory does not hold,
# because that inventory was written while the classes were still the exposure.
# Same columns as outputs/table9_test_inventory.csv, so the supplement can
# concatenate the two and count them together.
extra <- data.frame(
  family = c("first-order, breadth as exposure",
             "first-order, breadth as exposure",
             "breadth as exposure", "group"),
  test = c("breadth on NFE, all levels = 0 (ref = 1 domain)",
           "breadth on ResAll, all levels = 0 (ref = 1 domain)",
           "distress-adjusted model, all levels = 0 (ref = 1 domain)",
           "breadth coefficients equal across family cancer groups"),
  statistic = round(c(fo_n[1], fo_r[1], joint1k[1], grp1[1]), 3),
  df = c(fo_n[2], fo_r[2], joint1k[2], grp1[2]),
  p = signif(c(fo_n[3], fo_r[3], joint1k[3], grp1[3]), 4),
  status = "added 2026-09-06")
write.csv(extra, "outputs/table21_firstorder_breadth_tests.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------------------------- Table 1 rows, and Fig. 1 rates
bt <- table(data$nsupcat, data$vclass)
brd <- data.frame(level = LVL,
                  overall_n = as.integer(nscat_n),
                  overall_pct = round(100 * as.integer(nscat_n) / nrow(data), 1))
for (k in seq_len(K)) {
  brd[[paste0("class", k, "_n")]] <- as.integer(bt[, k])
  brd[[paste0("class", k, "_pct")]] <-
    round(100 * as.integer(bt[, k]) / sum(bt[, k]), 1)
}
write.csv(brd, "outputs/table26_breadth_by_class.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

DOM <- c("family", "friends", "romance", "work", "education", "leisure",
         "spiritual", "community", "health")
LAB <- c("Family", "Friends", "Romance", "Work", "Education", "Leisure",
         "Spiritual", "Community", "Health")
rate <- data.frame(
  domain = DOM, label = LAB,
  n_endorsed = sapply(DOM, function(x) sum(data[[paste0("v_", x)]])),
  n_total = nrow(data))
rate$pct_endorsed <- round(100 * rate$n_endorsed / rate$n_total, 1)
rate$p_endorsed <- round(rate$n_endorsed / rate$n_total, 4)
write.csv(rate, "outputs/table27_endorsement_rates.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

# the observed rates must reproduce the codebook written by 00_derive.R
cb <- read.csv("outputs/domain_codebook.csv", fileEncoding = "UTF-8")
cb <- cb[order(cb$domain_index), ]
stopifnot(identical(cb$name, DOM), all(abs(cb$pct_endorsed - rate$pct_endorsed) < 0.05))
cat("observed endorsement rates reproduce the codebook\n")

print(out, row.names = FALSE)
print(tests, row.names = FALSE)
print(cross, row.names = FALSE)

# ---------------------------------------------------------------------------
# The step from class to breadth, as a result table rather than a ledger entry
# (WRITING_RULES section 8: a result the Results section leans on must sit in a
# result table). Model A covariates, the three class dummies, and breadth
# entered as a continuous term. Reported as the second panel of Table S9.
#
#   output  outputs/table28_class_given_breadth.csv
# ---------------------------------------------------------------------------
CLSD <- paste0("vclass", 2:K)
data$nsup_z <- as.numeric(scale(data$nsup))
cls_terms <- function(prefix)
  paste(sprintf("%s%d*%s", prefix, 2:K, CLSD), collapse = " + ")
MA_ns <- paste0(MEAS, "\nWB ~ ", cls_terms("cN"), " + ", COV, " + nsup_z\n")
fitA_ns <- sem(MA_ns, data, estimator = "MLR", missing = "ML")

pe <- parameterEstimates(fitA_ns, standardized = TRUE)
pe <- pe[pe$op == "~" & pe$lhs == "WB" & pe$rhs %in% c(CLSD, "nsup_z"), ]
pe$term <- ifelse(pe$rhs == "nsup_z", "Breadth, per SD",
                  CLASS_LAB[as.integer(sub("vclass", "", pe$rhs))])
cgb <- data.frame(term = pe$term, b = round(pe$est, 4),
                  ci_low = round(pe$ci.lower, 4),
                  ci_high = round(pe$ci.upper, 4),
                  beta = round(pe$std.all, 3), p = signif(pe$pvalue, 4))
write.csv(cgb, "outputs/table28_class_given_breadth.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

w_joint <- wald(fitA_ns, sprintf("cN%d == 0", 2:K))
w_equal <- wald(fitA_ns, c("cN2 == cN3", "cN3 == cN4"))
cgb_tests <- data.frame(
  test = c("Class coefficients all zero, breadth adjusted",
           "Class coefficients equal, breadth adjusted"),
  statistic = round(c(w_joint[1], w_equal[1]), 3),
  df = c(w_joint[2], w_equal[2]),
  p = signif(c(w_joint[3], w_equal[3]), 4))
write.csv(cgb_tests, "outputs/table28b_class_given_breadth_tests.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
print(cgb, row.names = FALSE)
print(cgb_tests, row.names = FALSE)

# ---------------------------------------------------------------------------
# The second-order disturbance, estimated rather than fixed. The Limitations
# and the footnote of Table S4 both quote this, so it needs cells of its own
# instead of living only in a script comment.
#
#   output  outputs/table30_disturbance.csv
# ---------------------------------------------------------------------------
M1_free <- paste0(MEAS_BASE, "\nWB ~ ", terms_of("a", REF1), " + ", COV, "\n")
fit1_free <- sem(M1_free, data, estimator = "MLR", missing = "ML")
pf <- parameterEstimates(fit1_free)
pf <- pf[pf$op == "~~" & pf$lhs == "ResAll" & pf$rhs == "ResAll", ][1, ]

sf <- standardizedSolution(fit1)
sf <- sf[sf$op == "~" & sf$lhs == "WB" & sf$rhs %in% REF1, ]
sr <- standardizedSolution(fit1_free)
sr <- sr[sr$op == "~" & sr$lhs == "WB" & sr$rhs %in% REF1, ]
stopifnot(identical(sf$rhs, sr$rhs))
dist_out <- data.frame(
  quantity = c("Disturbance of the resource-oriented factor, estimated freely",
               "Standard error", "p",
               "Largest change in a standardised breadth coefficient when the disturbance is fixed at zero"),
  value = c(sprintf("%.3f", pf$est), sprintf("%.3f", pf$se),
            sprintf("%.3f", pf$pvalue),
            sprintf("%.4f", max(abs(sf$est.std - sr$est.std)))))
write.csv(dist_out, "outputs/table30_disturbance.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
print(dist_out, row.names = FALSE)

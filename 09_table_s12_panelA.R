# =============================================================================
# 09_table_s12_panelA.R -- Table S12 Panel A rebuilt on the DOMAIN-COUNT
# exposure.
#
# Why. Panel A of Table S12 carried the four-class exposure (class2, class3,
# class4) while the manuscript cites it, in step 8 and in the Results paragraph
# on psychological distress, for "full coefficients for models with and without
# distress" -- that is, for the seven-level domain-count model of Table 3.
# Panel A was left over from the version in which class membership was still
# the exposure. This script exports every predictor of the two models that
# Table 3 columns 1 and 3 already report, so the panel matches the text.
#
# NO NEW MODEL IS FITTED HERE. The model syntax is not retyped: lines 1-97 of
# 07_reference_one_domain.R are evaluated verbatim, which rebuilds the data and
# re-fits fit1 (primary) and fit1_k6 (primary + K6) exactly as delivered. The
# script then asserts that the six domain-count coefficients of each fit
# reproduce table22b_breadth_ref1_precise.csv to the last digit before writing
# anything.
#
#   input   07_reference_one_domain.R  (lines 1-97, evaluated as-is)
#           values_analysis.rds, lca_models.rds
#           outputs/table22b_breadth_ref1_precise.csv   (agreement check)
#   output  outputs/table36_s12_panelA_breadth.csv
# =============================================================================
BASE <- "."  # the repository root; run every script from there
setwd(BASE)
source("round_half_up.R")  # round() and signif() round half up from here on (see that file)

src <- readLines("07_reference_one_domain.R", encoding = "UTF-8")
eval(parse(text = paste(src[1:97], collapse = "\n")), envir = globalenv())

# ---------------------------------------------------------------- agreement
# fit1 and fit1_k6 must reproduce the delivered Table 3 columns exactly.
prec <- read.csv("outputs/table22b_breadth_ref1_precise.csv",
                 stringsAsFactors = FALSE)
chk <- function(fit, label) {
  pe <- parameterEstimates(fit)
  pe <- pe[pe$op == "~" & pe$lhs == "WB" & pe$rhs %in% REF1, ]
  pe$level <- c("0", "2", "3", "4", "5", "6+")[
    match(pe$rhs, paste0("nsc", c(0, 2:6)))]
  ref <- prec[prec$model == label, ]
  m <- merge(pe, ref, by = "level")
  stopifnot(nrow(m) == 6)
  d <- max(abs(m$est - m$b), abs(m$ci.lower - m$ci_low),
           abs(m$ci.upper - m$ci_high), abs(m$pvalue - m$p))
  cat(sprintf("  %-52s max abs deviation %.3e\n", label, d))
  stopifnot(d < 1e-10)
  invisible(d)
}
cat("Agreement with the delivered Table 3 columns:\n")
chk(fit1,    "Primary model (ref = 1 domain)")
chk(fit1_k6, "Adjusted for psychological distress (ref = 1 domain)")

# ------------------------------------------------------------------- export
# Row order: the six domain-count levels (reference = one endorsed domain),
# then the covariates in the order they enter COV, then K6.
ORD <- c(paste0("nsc", c(0, 2:6)),
         "disaster1", "disaster2", "fcancer21", "female1", "age_z", "income_z",
         "edu22", "edu23", "occupation22", "occupation23", "marriage2",
         "marriage3", "comobidity1", "fukushima1", "health22", "k6_z")

# Readable row labels. Reference categories, from 00_derive.R:
#   edu2 1 = junior/high school, occupation2 1 = self-employed,
#   marriage 1 = never married, health2 1 = low, disaster 0 = no damage,
#   comobidity 0 = none. Domain-count reference = one endorsed domain.
LAB <- c(
  nsc0 = "No domains endorsed",           nsc2 = "Two domains endorsed",
  nsc3 = "Three domains endorsed",        nsc4 = "Four domains endorsed",
  nsc5 = "Five domains endorsed",         nsc6 = "Six or more domains endorsed",
  disaster1 = "Disaster damage, not recovered",
  disaster2 = "Disaster damage, recovered",
  fcancer21 = "Family cancer history",    female1 = "Female",
  age_z = "Age, per SD",                  income_z = "Household income, per SD",
  edu22 = "Education beyond high school", edu23 = "Education, declined to answer",
  occupation22 = "Employee",              occupation23 = "Other occupation",
  marriage2 = "Married",                  marriage3 = "Divorced or widowed",
  comobidity1 = "Comorbidity",            fukushima1 = "Resident in Fukushima Prefecture",
  health22 = "Self-rated health, high",
  k6_z = "Psychological distress (K6), per SD")

grab <- function(fit) {
  pe <- parameterEstimates(fit, standardized = TRUE)
  pe <- pe[pe$op == "~" & pe$lhs == "WB", ]
  r <- pe[match(ORD, pe$rhs), c("rhs", "est", "ci.lower", "ci.upper",
                                "std.all", "pvalue")]
  # k6_z is absent from the model without distress; match() returns NA there,
  # which is the intended "n/a" row. Restore the label so the order is checkable.
  r$rhs <- ORD
  r
}
a <- grab(fit1)
b <- grab(fit1_k6)
stopifnot(identical(a$rhs, ORD), identical(b$rhs, ORD))
# every predictor except k6_z must be present in BOTH fits
stopifnot(!any(is.na(a$est[ORD != "k6_z"])), !any(is.na(b$est)))
stopifnot(is.na(a$est[ORD == "k6_z"]))

fmt_b  <- function(e, l, u) ifelse(is.na(e), "n/a",
  sprintf("%.3f [%.3f, %.3f]", round(e, 3), round(l, 3), round(u, 3)))
fmt_sd <- function(x) ifelse(is.na(x), "n/a", sprintf("%.3f", round(x, 3)))
# p is reported as the manuscript reports it: three decimals, or < 0.001.
fmt_p  <- function(p) ifelse(is.na(p), "n/a",
  ifelse(p < 0.0005, "<0.001", sprintf("%.3f", round(p, 3))))

out <- data.frame(
  predictor          = unname(LAB[ORD]),
  b_without_distress = fmt_b(a$est, a$ci.lower, a$ci.upper),
  std_without        = fmt_sd(a$std.all),
  p_without          = fmt_p(a$pvalue),
  b_with_distress    = fmt_b(b$est, b$ci.lower, b$ci.upper),
  std_with           = fmt_sd(b$std.all),
  p_with             = fmt_p(b$pvalue),
  change_in_std      = ifelse(is.na(a$std.all) | is.na(b$std.all), "n/a",
                              sprintf("%.3f", round(b$std.all - a$std.all, 3))),
  lavaan_rhs_raw     = ORD,
  stringsAsFactors   = FALSE)

write.csv(out, "outputs/table36_s12_panelA_breadth.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
cat("\nwrote outputs/table36_s12_panelA_breadth.csv,", nrow(out), "rows\n")
cat("p cells reading 0.000:",
    sum(c(out$p_without, out$p_with) == "0.000"), "\n")
print(out[, 1:8], row.names = FALSE)

# =============================================================================
# 10_sd_units_and_contrasts.R -- coefficients in standard deviation units, and
# contrasts between levels and between groups, from the models already fitted.
#
# NO NEW MODEL IS FITTED HERE. Every model is the delivered one: its syntax is
# not retyped but evaluated verbatim from the script that defines it, and each
# re-estimate is asserted to reproduce the delivered output before anything is
# written. The quantities below are read from those fits (lavInspect, vcov) or
# are linear combinations of their coefficients.
#
# 1. The divisor, SD_WB. The model-implied variance of overall well-being in
#    the primary model (Eq. (4), N = 2179, reference = one endorsed domain),
#    lavInspect(fit, "cov.lv")["WB", "WB"], and its square root. One divisor
#    serves every table, every column and every group: under strict
#    measurement invariance the latent scale is common, and the comparisons
#    between groups are made on that scale. The SD of overall well-being in
#    every other model and group is recorded beside it, with its relative
#    difference from the primary value, so that the choice is on the record.
#    Check: the std.lv coefficients of the primary model times SD_WB return the
#    raw coefficients to 1e-6.
#
# 2. SD units. b / SD_WB and CI / SD_WB, for Table 3 (three columns), Table 4
#    (two groups), Table 5 (two groups), Table S5 Panel B, Table S6 and Table
#    S12 Panel A (both columns, and the change between them). For the primary
#    model the delta-method interval of the std.lv coefficient is written
#    beside the simple quotient, to show the two do not differ at the printed
#    precision.
#
# 3. Contrasts, raw and in SD units, with SE and 95% CI, from vcov():
#    step 2 model (Eq. (3)): delta3 - delta2, delta3 - delta4, delta2 - delta4
#    primary model (Eq. (4)): beta5 - beta4, beta6+ - beta4, beta6+ - beta5
#    Table 4 multiple-group model: beta_k(history) - beta_k(no history)
#    Table 5 two separate fits: beta_k(damage) - beta_k(no damage),
#      SE = sqrt(se1^2 + se0^2) because the fits are independent
#
# 4. Table S7 reproduced: std.lv of the two-component model (step 5), which is
#    what the "Difference in SD units" columns of Table S7 carry.
#
#   input   07_reference_one_domain.R  (lines 1-97, then the step 2 and
#                                       two-component blocks, evaluated as-is)
#           02_structural.R            (setup, step 2 model, matched block)
#           08_disaster_modification.R (everything up to the step 7 fits)
#           values_analysis.rds, lca_models.rds
#           outputs/table22b_breadth_ref1_precise.csv, table28_class_given_
#           breadth.csv, table10_breadth_matched.csv, table31b_disaster_
#           precise.csv, table36_s12_panelA_breadth.csv  (agreement checks)
#   output  outputs/table37_sd_units.csv          every SD-unit cell
#           outputs/table38_sd_wb.csv             SD of overall well-being,
#                                                 every model and group
#           outputs/table39_contrasts.csv         contrasts, raw and SD units
#           outputs/table40_stdlv_delta.csv       primary model, quotient
#                                                 against delta-method std.lv
#           outputs/table41_s7_stdlv.csv          Table S7 SD-unit columns
# =============================================================================
suppressPackageStartupMessages({ library(lavaan) })
options(stringsAsFactors = FALSE)
BASE <- "."  # the repository root; run every script from there
setwd(BASE)

# ------------------------------------------------------------ evaluation
# Evaluate lines of a delivered script in an environment of their own, so the
# three scripts' objects (each has its own `data`, `grab`, `MEAS`) cannot
# collide. Blocks are located by their first and last line, each of which must
# occur exactly once, so a later edit that moves lines cannot shift the block.
src_of <- function(file) readLines(file, encoding = "UTF-8")
line_of <- function(src, pattern) {
  i <- grep(pattern, src)
  if (length(i) != 1) stop("pattern matched ", length(i), " lines: ", pattern)
  i
}
run_lines <- function(src, idx, env) {
  eval(parse(text = paste(src[idx], collapse = "\n")), envir = env)
  invisible(env)
}

# ---- 07: data, the primary model and its two adjusted columns, the multiple-
# group model by family cancer history, the step 2 model, the two-component
# model. Lines 1-97 are what 09_table_s12_panelA.R evaluates.
s07 <- src_of("07_reference_one_domain.R")
e07 <- new.env(parent = globalenv())
run_lines(s07, 1:97, e07)
run_lines(s07, line_of(s07, "^CLSD <- paste0\\(\"vclass\", 2:K\\)"):
               line_of(s07, "^fitA_ns <- sem\\(MA_ns"), e07)
run_lines(s07, line_of(s07, "^MFO1 <- paste\\("):
               line_of(s07, "^fitFO1 <- sem\\(MFO1"), e07)

# ---- 02: the same step 2 model as fitted there, and the three matched fits of
# Table S6. The command-line handling is skipped and K given directly.
s02 <- src_of("02_structural.R")
e02 <- new.env(parent = globalenv())
run_lines(s02, line_of(s02, "^suppressPackageStartupMessages"):
               (line_of(s02, "^args <- commandArgs") - 1), e02)
assign("K", 4L, envir = e02)
assign("CLASS_LAB", c("Minimal", "Relational-occupational", "Family-health",
                      "Expansive"), envir = e02)
run_lines(s02, line_of(s02, "^t_start <- Sys.time"):
               (line_of(s02, "^# Model A -- main") - 2), e02)
run_lines(s02, line_of(s02, "^mnsub <- data\\[data\\$nsup %in% c\\(3, 4\\)"):
               line_of(s02, "^row_m123 <- struct_rows\\(fit_m123\\)"), e02)
run_lines(s02, c(line_of(s02, "^data\\$nsup_z <- as.numeric"),
                 line_of(s02, "^MA_ns   <- paste0"),
                 line_of(s02, "^fitA_ns   <- sem\\(MA_ns")), e02)

# ---- 08: the step 7 fits. `grab` is wrapped so that the fitted objects are
# kept; the lines that call it are evaluated unchanged.
s08 <- src_of("08_disaster_modification.R")
e08 <- new.env(parent = globalenv())
i_chk <- line_of(s08, "reproduction check, step 6")
run_lines(s08, 1:(i_chk - 1), e08)
local({
  FITS <- list()
  grab_orig <- grab
  grab <- function(fit) { FITS[[length(FITS) + 1]] <<- fit; grab_orig(fit) }
}, envir = e08)
run_lines(s08, i_chk:line_of(s08, "^de1 <- grab\\("), e08)
FITS08 <- get("FITS", envir = e08)
stopifnot(length(FITS08) == 4)              # fc0, fc1, de0, de1, in that order
fit_de0 <- FITS08[[3]]
fit_de1 <- FITS08[[4]]

setwd(BASE)                                 # the evaluated scripts set it too

fit1     <- e07$fit1
fit1_com <- e07$fit1_com
fit1_k6  <- e07$fit1_k6
fit1_g   <- e07$fit1_g
fitFO1   <- e07$fitFO1
fitA_ns  <- e02$fitA_ns
REF1     <- e07$REF1
LVL      <- e07$LVL

# ------------------------------------------------------------- agreement
# Each re-estimate must reproduce what was delivered before anything is used.
TOL <- 1e-10
coef_rows <- function(fit, keep, group = NULL, lhs = "WB") {
  pe <- parameterEstimates(fit)
  pe <- pe[pe$op == "~" & pe$lhs == lhs & pe$rhs %in% keep, ]
  if (!is.null(group)) pe <- pe[pe$group == group, ]
  pe[match(keep[keep %in% pe$rhs], pe$rhs), ]
}
level_of <- function(rhs) LVL[as.integer(sub("nsc", "", rhs)) + 1]

prec <- read.csv("outputs/table22b_breadth_ref1_precise.csv")
agree <- function(fit, label, keep, group = NULL) {
  r <- coef_rows(fit, keep, group)
  p <- prec[prec$model == label, ]
  p <- p[match(level_of(r$rhs), p$level), ]
  d <- max(abs(r$est - p$b), abs(r$ci.lower - p$ci_low),
           abs(r$ci.upper - p$ci_high), abs(r$pvalue - p$p))
  cat(sprintf("  %-55s max abs deviation %.2e\n", label, d))
  stopifnot(nrow(r) == nrow(p), d < TOL)
}
cat("Agreement with delivered outputs:\n")
agree(fit1,     "Primary model (ref = 1 domain)", REF1)
agree(fit1_com, "Adjusted for commitment intensity (ref = 1 domain)",
      paste0("nsc", 2:6))
agree(fit1_k6,  "Adjusted for psychological distress (ref = 1 domain)", REF1)
agree(fit1_g,   "No family cancer history (ref = 1 domain)", REF1, group = 1)
agree(fit1_g,   "Family cancer history (ref = 1 domain)",    REF1, group = 2)

# group 1 of the multiple-group fit is the group without a family cancer history
glab <- lavInspect(fit1_g, "group.label")
gn   <- lavInspect(fit1_g, "nobs")
stopifnot(identical(as.character(glab), c("0", "1")), gn[1] == 1324,
          gn[2] == 855)

precd <- read.csv("outputs/table31b_disaster_precise.csv")
for (k in 1:2) {
  f   <- list(fit_de0, fit_de1)[[k]]
  lab <- c("No disaster damage (ref = 1 domain)",
           "Disaster damage (ref = 1 domain)")[k]
  r <- coef_rows(f, REF1)
  p <- precd[precd$model == lab, ]
  p <- p[match(level_of(r$rhs), p$level), ]
  d <- max(abs(r$est - p$b), abs(r$ci.lower - p$ci_low),
           abs(r$ci.upper - p$ci_high), abs(r$pvalue - p$p))
  cat(sprintf("  %-55s max abs deviation %.2e\n", lab, d))
  stopifnot(d < TOL)
}
stopifnot(lavInspect(fit_de0, "nobs") == 1712, lavInspect(fit_de1, "nobs") == 467)

# step 2: 02 and 07 fit the same syntax; both must agree with the delivered
# Table S5 Panel B to its four printed decimals, and with each other exactly
CLSD <- paste0("vclass", 2:4)
ra <- coef_rows(fitA_ns, c(CLSD, "nsup_z"))
rb <- coef_rows(e07$fitA_ns, c(CLSD, "nsup_z"))
d_ab <- max(abs(ra$est - rb$est), abs(ra$se - rb$se))
cat(sprintf("  %-55s max abs deviation %.2e\n",
            "step 2 model, 02 against 07", d_ab))
stopifnot(d_ab < TOL)
t28 <- read.csv("outputs/table28_class_given_breadth.csv")
stopifnot(max(abs(round(ra$est, 4) - t28$b)) == 0,
          max(abs(round(ra$ci.lower, 4) - t28$ci_low)) == 0)

# Table S6: the three matched fits against the delivered full-precision file
t10 <- read.csv("outputs/table10_breadth_matched.csv")
mfits <- list(e02$fit_mA, e02$fit_mA2b, e02$fit_m123)
for (k in 1:3) {
  r <- coef_rows(mfits[[k]], "c3")
  d <- max(abs(r$est - t10$b[k]), abs(r$ci.lower - t10$ci_low[k]),
           abs(r$ci.upper - t10$ci_high[k]))
  cat(sprintf("  %-55s max abs deviation %.2e\n", t10$model[k], d))
  stopifnot(d < TOL, lavInspect(mfits[[k]], "nobs") == t10$n[k])
}

# ------------------------------------------------------------- 1. SD_WB
var_wb <- function(fit) {
  cl <- lavInspect(fit, "cov.lv")
  if (is.list(cl)) sapply(cl, function(z) z["WB", "WB"]) else cl["WB", "WB"]
}
VAR_WB <- var_wb(fit1)
SD_WB  <- sqrt(VAR_WB)
cat(sprintf("\nSD_WB (primary model) = %.10f  (variance %.10f)\n", SD_WB,
            VAR_WB))

# std.lv check: standardised on the latent variance only, so times SD_WB it
# must give back the raw coefficient
slv <- standardizedSolution(fit1, type = "std.lv")
slv <- slv[slv$op == "~" & slv$lhs == "WB" & slv$rhs %in% REF1, ]
r1  <- coef_rows(fit1, REF1)
slv <- slv[match(r1$rhs, slv$rhs), ]
d_lv <- max(abs(slv$est.std * SD_WB - r1$est))
cat(sprintf("std.lv x SD_WB against raw b, max abs difference %.2e\n", d_lv))
stopifnot(d_lv < 1e-6)

sdtab <- rbind(
  data.frame(model = "Primary model (Table 3, step 3)", group = "all",
             n = lavInspect(fit1, "nobs"), var_wb = VAR_WB),
  data.frame(model = "Adjusted for commitment intensity (Table 3, step 4)",
             group = "all", n = lavInspect(fit1_com, "nobs"),
             var_wb = var_wb(fit1_com)),
  data.frame(model = "Adjusted for psychological distress (Table 3, step 8)",
             group = "all", n = lavInspect(fit1_k6, "nobs"),
             var_wb = var_wb(fit1_k6)),
  data.frame(model = "Family cancer history, multiple-group (Table 4, step 6)",
             group = c("No family cancer history", "Family cancer history"),
             n = as.integer(gn), var_wb = unname(var_wb(fit1_g))),
  data.frame(model = "Disaster exposure, two separate fits (Table 5, step 7)",
             group = c("No disaster damage", "Disaster damage"),
             n = c(lavInspect(fit_de0, "nobs"), lavInspect(fit_de1, "nobs")),
             var_wb = c(var_wb(fit_de0), var_wb(fit_de1))),
  data.frame(model = "Class and standardised count (Table S5 Panel B, step 2)",
             group = "all", n = lavInspect(fitA_ns, "nobs"),
             var_wb = var_wb(fitA_ns)),
  data.frame(model = c("Matched, primary adjustment set (Table S6 row 1)",
                       "Matched, adding commitment intensity (Table S6 row 2)",
                       "Matched, retaining the Minimal class (Table S6 row 3)"),
             group = "three or four domains",
             n = sapply(mfits, lavInspect, "nobs"),
             var_wb = sapply(mfits, var_wb)))
sdtab$sd_wb <- sqrt(sdtab$var_wb)
sdtab$rel_diff_vs_primary <- sdtab$sd_wb / SD_WB - 1
write.csv(sdtab, "outputs/table38_sd_wb.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
cat("\nSD of overall well-being, every model and group:\n")
print(transform(sdtab, sd_wb = round(sd_wb, 5),
                rel_diff_vs_primary = round(rel_diff_vs_primary, 4))[
  , c("model", "group", "n", "sd_wb", "rel_diff_vs_primary")], row.names = FALSE)
cat(sprintf("largest relative difference from the primary SD: %.4f\n",
            max(abs(sdtab$rel_diff_vs_primary))))

# -------------------------------------------------------- 2. SD units
f3 <- function(x) {
  s <- sprintf("%.3f", x)
  ifelse(s == "-0.000", "0.000", s)
}
cell3 <- function(e, l, u) sprintf("%s [%s, %s]", f3(e), f3(l), f3(u))

sd_rows <- function(tab, column, fit, keep, group = NULL, labels = NULL) {
  r <- coef_rows(fit, keep, group)
  data.frame(table = tab, column = column,
             row = if (is.null(labels)) level_of(r$rhs) else labels[r$rhs],
             lavaan_rhs = r$rhs, b = r$est, ci_low = r$ci.lower,
             ci_high = r$ci.upper, p = r$pvalue, divisor = SD_WB,
             sd = r$est / SD_WB, sd_low = r$ci.lower / SD_WB,
             sd_high = r$ci.upper / SD_WB,
             b_cell = cell3(r$est, r$ci.lower, r$ci.upper),
             sd_cell = cell3(r$est / SD_WB, r$ci.lower / SD_WB,
                             r$ci.upper / SD_WB),
             sd_point = f3(r$est / SD_WB))
}
# Table S5 Panel B row labels, as printed
S5_LAB <- c(vclass2 = "Relational-occupational", vclass3 = "Family-health",
            vclass4 = "Expansive",
            nsup_z = "Number of endorsed domains, per SD")
units <- rbind(
  sd_rows("Table 3", "Primary model", fit1, REF1),
  sd_rows("Table 3", "Adjusted for commitment intensity", fit1_com,
          paste0("nsc", 2:6)),
  sd_rows("Table 3", "Adjusted for psychological distress", fit1_k6, REF1),
  sd_rows("Table 4", "No family cancer history", fit1_g, REF1, group = 1),
  sd_rows("Table 4", "Family cancer history", fit1_g, REF1, group = 2),
  sd_rows("Table 5", "No disaster damage", fit_de0, REF1),
  sd_rows("Table 5", "Disaster damage", fit_de1, REF1),
  sd_rows("Table S5 Panel B", "Difference in SD units", fitA_ns,
          c(CLSD, "nsup_z"), labels = S5_LAB))
S6_LAB <- c("Relational-occupational versus Family-health, primary adjustment set",
            "Relational-occupational versus Family-health, adding commitment intensity",
            "As above, retaining the Minimal class in the model")
for (k in 1:3)
  units <- rbind(units, sd_rows("Table S6", "Difference in SD units",
                                mfits[[k]], "c3",
                                labels = setNames(S6_LAB[k], "c3")))

# Table S12 Panel A: every predictor of the two Table 3 models, both columns,
# and the change between them on the same divisor
ORD12 <- c(paste0("nsc", c(0, 2:6)),
           "disaster1", "disaster2", "fcancer21", "female1", "age_z",
           "income_z", "edu22", "edu23", "occupation22", "occupation23",
           "marriage2", "marriage3", "comobidity1", "fukushima1", "health22",
           "k6_z")
t36 <- read.csv("outputs/table36_s12_panelA_breadth.csv")
stopifnot(identical(t36$lavaan_rhs_raw, ORD12))
LAB12 <- setNames(t36$predictor, ORD12)
u_wo <- sd_rows("Table S12 Panel A", "Without distress", fit1,
                ORD12[ORD12 != "k6_z"], labels = LAB12)
u_wi <- sd_rows("Table S12 Panel A", "With distress", fit1_k6, ORD12,
                labels = LAB12)
units <- rbind(units, u_wo, u_wi)
m12 <- merge(u_wo[, c("lavaan_rhs", "b")], u_wi[, c("lavaan_rhs", "b")],
             by = "lavaan_rhs", suffixes = c("_without", "_with"))
m12 <- m12[match(ORD12[ORD12 != "k6_z"], m12$lavaan_rhs), ]
chg <- data.frame(table = "Table S12 Panel A", column = "Change in SD units",
                  row = unname(LAB12[m12$lavaan_rhs]),
                  lavaan_rhs = m12$lavaan_rhs, b = m12$b_with - m12$b_without,
                  ci_low = NA, ci_high = NA, p = NA, divisor = SD_WB,
                  sd = (m12$b_with - m12$b_without) / SD_WB, sd_low = NA,
                  sd_high = NA, b_cell = NA, sd_cell = NA,
                  sd_point = f3((m12$b_with - m12$b_without) / SD_WB))
units <- rbind(units, chg)
# the S12 b column must reproduce the delivered panel to its printed digits
# (the one hand-widened cell, 0.0005, is compared on its three-digit form)
chk12 <- function(u, col) {
  got <- u$b_cell[match(ORD12, u$lavaan_rhs)]
  want <- t36[[col]]
  ok <- is.na(got) & want == "n/a" | (!is.na(got) & got == want)
  stopifnot(all(ok))
}
chk12(u_wo, "b_without_distress")
chk12(u_wi, "b_with_distress")
write.csv(units, "outputs/table37_sd_units.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
cat(sprintf("\nwrote outputs/table37_sd_units.csv, %d rows\n", nrow(units)))

# --------------------------------------- primary model: quotient vs std.lv
slv_ci <- data.frame(
  level = level_of(r1$rhs), b = r1$est, b_low = r1$ci.lower,
  b_high = r1$ci.upper,
  quotient = r1$est / SD_WB, quotient_low = r1$ci.lower / SD_WB,
  quotient_high = r1$ci.upper / SD_WB,
  stdlv = slv$est.std, stdlv_se = slv$se, stdlv_low = slv$ci.lower,
  stdlv_high = slv$ci.upper)
slv_ci$diff_low  <- slv_ci$stdlv_low - slv_ci$quotient_low
slv_ci$diff_high <- slv_ci$stdlv_high - slv_ci$quotient_high
slv_ci$quotient_cell <- cell3(slv_ci$quotient, slv_ci$quotient_low,
                              slv_ci$quotient_high)
slv_ci$stdlv_cell <- cell3(slv_ci$stdlv, slv_ci$stdlv_low, slv_ci$stdlv_high)
slv_ci$same_at_3dp <- slv_ci$quotient_cell == slv_ci$stdlv_cell
write.csv(slv_ci, "outputs/table40_stdlv_delta.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
cat("\nPrimary model, simple quotient against delta-method std.lv:\n")
print(slv_ci[, c("level", "quotient_cell", "stdlv_cell", "diff_low",
                 "diff_high")], row.names = FALSE)

# -------------------------------------------------------- 3. contrasts
Z <- qnorm(0.975)
lin <- function(fit, w) {
  cf <- coef(fit)
  V  <- vcov(fit)
  stopifnot(all(names(w) %in% names(cf)), all(names(w) %in% colnames(V)))
  # vcov() must be the matrix behind the reported robust standard errors
  pe <- parameterEstimates(fit)
  se_pe <- pe$se[match(names(w), pe$label)]
  stopifnot(max(abs(sqrt(diag(V)[names(w)]) - se_pe)) < 1e-10)
  est <- sum(w * cf[names(w)])
  se  <- sqrt(as.numeric(t(w) %*% V[names(w), names(w)] %*% w))
  c(est = est, se = se)
}
con_row <- function(section, model, contrast, est, se) {
  data.frame(section = section, model = model, contrast = contrast,
             est = est, se = se, ci_low = est - Z * se, ci_high = est + Z * se,
             p = 2 * pnorm(-abs(est / se)),
             est_sd = est / SD_WB, se_sd = se / SD_WB,
             ci_low_sd = (est - Z * se) / SD_WB,
             ci_high_sd = (est + Z * se) / SD_WB)
}
cons <- list()
# step 2, labels cN2 (Relational-occupational), cN3 (Family-health),
# cN4 (Expansive); the reference is Minimal
for (cc in list(c("cN3", "cN2", "delta3 - delta2 (Family-health minus Relational-occupational)"),
                c("cN3", "cN4", "delta3 - delta4 (Family-health minus Expansive)"),
                c("cN2", "cN4", "delta2 - delta4 (Relational-occupational minus Expansive)"))) {
  v <- lin(fitA_ns, setNames(c(1, -1), cc[1:2]))
  cons[[length(cons) + 1]] <- con_row("step 2 model, Eq. (3)",
                                      "Class and standardised count",
                                      cc[3], v["est"], v["se"])
}
# primary model, labels a4, a5, a6 (a6 = six or more)
for (cc in list(c("a5", "a4", "beta5 - beta4"), c("a6", "a4", "beta6+ - beta4"),
                c("a6", "a5", "beta6+ - beta5"))) {
  v <- lin(fit1, setNames(c(1, -1), cc[1:2]))
  cons[[length(cons) + 1]] <- con_row("primary model, Eq. (4)",
                                      "Primary model", cc[3], v["est"], v["se"])
}
# Table 4: family cancer history minus none, at each level; group 1 = none
for (k in c(0, 2:6)) {
  v <- lin(fit1_g, setNames(c(1, -1), c(sprintf("g1_%d", k), sprintf("g0_%d", k))))
  cons[[length(cons) + 1]] <- con_row(
    "Table 4 multiple-group model", "Family cancer history minus none",
    sprintf("level %s", LVL[k + 1]), v["est"], v["se"])
}
# Table 5: damage minus no damage, independent fits
r0 <- coef_rows(fit_de0, REF1)
r1d <- coef_rows(fit_de1, REF1)
stopifnot(identical(r0$rhs, r1d$rhs))
for (i in seq_len(nrow(r0))) {
  cons[[length(cons) + 1]] <- con_row(
    "Table 5 two separate fits", "Disaster damage minus no damage",
    sprintf("level %s", level_of(r0$rhs[i])), r1d$est[i] - r0$est[i],
    sqrt(r1d$se[i]^2 + r0$se[i]^2))
}
cons <- do.call(rbind, cons)
rownames(cons) <- NULL
write.csv(cons, "outputs/table39_contrasts.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
cat("\nContrasts (raw, then SD units):\n")
print(transform(cons, est = round(est, 4), se = round(se, 4),
                ci_low = round(ci_low, 4), ci_high = round(ci_high, 4),
                p = signif(p, 3), est_sd = round(est_sd, 3),
                ci_low_sd = round(ci_low_sd, 3),
                ci_high_sd = round(ci_high_sd, 3))[
  , c("model", "contrast", "est", "se", "ci_low", "ci_high", "p", "est_sd",
      "ci_low_sd", "ci_high_sd")], row.names = FALSE)

# ------------------------------------------------- 4. Table S7, std.lv
s7 <- standardizedSolution(fitFO1, type = "std.lv")
s7 <- s7[s7$op == "~" & s7$lhs %in% c("NFE", "ResAll") & s7$rhs %in% REF1, ]
s7out <- data.frame(component = s7$lhs, level = level_of(s7$rhs),
                    stdlv = s7$est.std, low = s7$ci.lower, high = s7$ci.upper,
                    p = s7$pvalue,
                    cell = cell3(s7$est.std, s7$ci.lower, s7$ci.upper))
write.csv(s7out, "outputs/table41_s7_stdlv.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
cat("\nTable S7, SD units (std.lv):\n")
print(s7out[, c("component", "level", "cell")], row.names = FALSE)

cat("\ndone -- outputs/table37 to table41\n")

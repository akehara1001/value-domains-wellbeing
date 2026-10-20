# =============================================================================
# 02_structural.R -- value class -> well-being (Tables 1-4, contrasts, K6 pair)
# =============================================================================
# PHASE 2. Runs only after the number of classes is fixed. K is a REQUIRED
# command-line argument; this script never picks one.
#
#   Rscript 02_structural.R <K> [optional,comma,separated,class,labels]
#
#   input   values_analysis.rds   (00_derive.R)
#           lca_models.rds        (01_lca.R, classes size-ordered)
#   output  outputs/results_phase2.md   every number, with its source
#           outputs/table1_background.csv
#           outputs/table2_class_profile.csv
#           outputs/table3_modelA.csv
#           outputs/table4_modelB.csv
#           outputs/table5_k6_contrast.csv     K6 in / out, all predictors
#           outputs/table6_firstorder.csv      class effects on NFE and ResAll
#           outputs/table7_slope_per_domain.csv
#           outputs/contrasts_class2_vs_3.csv
#           outputs/sensitivity.csv
#           outputs/fig1_class_profiles.png
#           structural_fits.rds
#
# MEASUREMENT MODEL. Paper 1's, unchanged and not re-explored:
#   NFE    =~ sadR + worryR              negative-free emotion
#   ResAll =~ ls + gh + es + vital + comp  resource-oriented well-being
#   WB     =~ NFE + ResAll               second-order well-being
#   ls ~~ gh ; es ~~ vital
#
# ResAll DISTURBANCE FIXED AT ZERO in the main models. The second-order factor
# has only two indicators, so the disturbance of ResAll is weakly determined and
# in the free parameterisation it lands slightly below zero. Measured free it is
# indistinguishable from zero in every K6-free specification:
#   Model A -0.049 (SE 0.260, p 0.852)   A2a  0.178 (SE 0.202, p 0.376)
#   A2b      0.024 (SE 0.164, p 0.882)   interaction -0.064 (SE 0.262, p 0.808)
#   Model B  group 1 -0.049 (p 0.875), group 2 0.030 (p 0.941)
# Fixing it removes the Heywood case without moving the substance: the largest
# change in a standardised class coefficient is 0.0005.
#
# THE ONE EXCEPTION is the K6 model, where the free estimate is 0.771
# (SE 0.102, p < 0.001) -- solidly non-zero. Fixing it there would be a genuine
# misspecification, so the K6 model keeps the disturbance FREE. The two main-text
# columns therefore differ in this one respect, which is stated in the table.
#
# K6 is not in the main model (it is near-collinear with well-being and
# plausibly lies on the value -> well-being path). Following the 2026-08-30
# instruction the K6-adjusted model is now reported ALONGSIDE the main model in
# the main text rather than in the supplement.
#
# The reference class is class 1, which 01_lca.R defines as the largest class.
# =============================================================================

suppressPackageStartupMessages({
  library(lavaan); library(dplyr); library(poLCA); library(car); library(psych)
})
options(stringsAsFactors = FALSE)

BASE <- "."  # the repository root; run every script from there
setwd(BASE)
source("round_half_up.R")  # round() and signif() round half up from here on (see that file)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || is.na(suppressWarnings(as.integer(args[1])))) {
  stop("usage: Rscript 02_structural.R <K> [class labels]\n",
       "  K must be supplied explicitly; this script does not choose it.")
}
K <- as.integer(args[1])
CLASS_LAB <- if (length(args) >= 2) trimws(strsplit(args[2], ",")[[1]]) else
  paste("Class", seq_len(K))
stopifnot(K >= 2, length(CLASS_LAB) == K)

t_start <- Sys.time()
d  <- readRDS("values_analysis.rds")
LC <- readRDS("lca_models.rds")
stopifnot(K <= length(LC$models))
m  <- LC$models[[K]]
DOM <- LC$domains; DOM_LAB <- LC$dom_lab; IND <- LC$indicators

# -----------------------------------------------------------------------------
# Analysis frame
# -----------------------------------------------------------------------------
data <- as.data.frame(d)
data$vclass <- factor(m$predclass, levels = seq_len(K))
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

# Guard the disaster dummy naming, which is the point of instruction 1.
# 00_derive.R line 153 codes 0 = unexposed, 1 = UNRECOVERED, 2 = RECOVERED, so
# factor level 1 -> dummy "disaster1" is UNRECOVERED (n = 95) and level 2 ->
# "disaster2" is RECOVERED (n = 372). The ordering is counter-intuitive and is
# inherited from paper 1, so it is asserted rather than trusted.
stopifnot(sum(data$disaster1) == 95, sum(data$disaster2) == 372)

CLSD <- paste0("vclass", 2:K)                      # class dummies, ref = class 1
COV  <- paste("disaster1 + disaster2 + fcancer21 + female1 + age_z + income_z +",
              "edu22 + edu23 + occupation22 + occupation23 + marriage2 +",
              "marriage3 + comobidity1 + fukushima1 + health22")

MEAS_BASE <- paste("NFE =~ sadR + worryR",
                   "ResAll =~ ls + gh + es + vital + comp",
                   "WB =~ NFE + ResAll", "ls ~~ gh", "es ~~ vital", sep = "\n")
MEAS      <- paste0(MEAS_BASE, "\nResAll ~~ 0*ResAll")   # main specification
MEAS_FREE <- MEAS_BASE                                    # K6 model keeps it free

# -----------------------------------------------------------------------------
# Readable term names. Instruction 1: disaster1/disaster2 must never appear in
# an output table under a name that hides which is which.
# -----------------------------------------------------------------------------
TMAP <- c(setNames(sprintf("class%d_%s", 2:K, CLASS_LAB[2:K]), CLSD),
          setNames(sprintf("classO%d", 2:K), paste0("vclassO", 2:K)),
          disaster1 = "disaster_unrecovered", disaster2 = "disaster_recovered",
          fcancer21 = "family_cancer", female1 = "female",
          age_z = "age_z", income_z = "income_z",
          edu22 = "edu_beyond_highschool", edu23 = "edu_declined",
          occupation22 = "occupation_employee", occupation23 = "occupation_other",
          marriage2 = "married", marriage3 = "divorced_widowed",
          comobidity1 = "comorbidity", fukushima1 = "fukushima_resident",
          health22 = "subjective_health_high", k6_z = "k6_z",
          mcom_z = "mean_commitment_z")
relab <- function(x) unname(ifelse(x %in% names(TMAP), TMAP[x], x))

cls_terms <- function(prefix) paste(sprintf("%s%d*%s", prefix, 2:K, CLSD),
                                    collapse = " + ")
wald_zero <- function(prefix) paste(sprintf("%s%d == 0", prefix, 2:K),
                                    collapse = "\n")
# Paper 1's VIF convention: car::vif()'s GVIF^(1/(2*Df)) column, which for a
# single-df term equals sqrt(VIF) and is on the standard-error inflation scale.
vif_pub <- function(fit) {
  v <- car::vif(fit)
  if (is.matrix(v)) v[, "GVIF^(1/(2*Df))"] else sqrt(v)
}

fmt_ci <- function(est, lo, hi, dp = 3)
  sprintf("%.*f [%.*f, %.*f]", dp, round(est, dp), dp, round(lo, dp), dp, round(hi, dp))
fmt_p <- function(p) ifelse(is.na(p), "—", ifelse(p < 0.001, "<0.001", sprintf("%.3f", round(p, 3))))

struct_rows <- function(fit, group = NULL, lhs = "WB") {
  pe <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
  pe <- pe[pe$op == "~" & pe$lhs == lhs, ]
  if (!is.null(group) && "group" %in% names(pe)) pe <- pe[pe$group == group, ]
  pe$term <- relab(pe$rhs)
  pe[, c("rhs", "term", "est", "se", "z", "pvalue", "ci.lower", "ci.upper", "std.all")]
}
def_rows <- function(fit, group = NULL) {
  pe <- parameterEstimates(fit, ci = TRUE)
  pe <- pe[pe$op == ":=", ]
  pe[, c("label", "est", "se", "z", "pvalue", "ci.lower", "ci.upper")]
}
getdef <- function(fit, lab) {
  r <- def_rows(fit); r[r$label == lab, ]
}

# contrast definitions used everywhere: class2 vs class3, and the adjacent
# steps needed for the per-domain slope table
defs <- function(p) paste0(
  "\n", p, "_21 := ", p, "2",
  "\n", p, "_32 := ", p, "3 - ", p, "2",
  if (K >= 4) paste0("\n", p, "_42 := ", p, "4 - ", p, "2",
                     "\n", p, "_43 := ", p, "4 - ", p, "3") else "", "\n")

# H0: the non-reference class effects are all equal, i.e. b2 = b3 = b4 (df K-2).
# Rejecting it means the classes are ordered/graded; failing to reject means the
# data support only "class 1 differs from the rest" -- a minimal-pair split.
wald_equal <- function(prefix)
  paste(sprintf("%s%d == %s%d", prefix, 2:(K - 1), prefix, 3:K), collapse = "\n")

# =============================================================================
# Model A -- main, all 2179
# =============================================================================
MA <- paste0(MEAS, "\nWB ~ ", cls_terms("cA"), " + ", COV, "\n", defs("cA"))
fitA <- sem(MA, data, estimator = "MLR", missing = "ML")
rowsA <- struct_rows(fitA)
waldA <- lavTestWald(fitA, constraints = wald_zero("cA"))
waldA_eq <- lavTestWald(fitA, constraints = wald_equal("cA"))
fitA_meas <- fitMeasures(fitA, c("chisq", "df", "cfi.robust", "rmsea.robust", "srmr"))

# the free-disturbance version, kept only to report the note required by
# instruction 4 and to show that fixing it changes nothing
MA_free <- paste0(MEAS_FREE, "\nWB ~ ", cls_terms("cF"), " + ", COV, "\n")
fitA_free <- sem(MA_free, data, estimator = "MLR", missing = "ML")
rowsA_free <- struct_rows(fitA_free)
pe_free <- parameterEstimates(fitA_free)
res_free <- pe_free[pe_free$op == "~~" & pe_free$lhs == "ResAll" &
                      pe_free$rhs == "ResAll", ]
fitA_free_meas <- fitMeasures(fitA_free, c("cfi.robust", "rmsea.robust", "srmr"))

# =============================================================================
# Model A2 -- the Wada 2022 discriminant test
# =============================================================================
dsub <- data[data$nsup > 0, ]
stopifnot(nrow(dsub) == 1964, !anyNA(dsub$mcom_z))
MA2a <- paste0(MEAS, "\nWB ~ ", cls_terms("cB"), " + ", COV, "\n", defs("cB"))
MA2b <- paste0(MEAS, "\nWB ~ ", cls_terms("cC"), " + ", COV, " + mcom_z\n", defs("cC"))
fitA2a <- sem(MA2a, dsub, estimator = "MLR", missing = "ML")
fitA2b <- sem(MA2b, dsub, estimator = "MLR", missing = "ML")
rowsA2a <- struct_rows(fitA2a); rowsA2b <- struct_rows(fitA2b)
waldA2a <- lavTestWald(fitA2a, constraints = wald_zero("cB"))
waldA2b <- lavTestWald(fitA2b, constraints = wald_zero("cC"))
waldA2b_eq <- lavTestWald(fitA2b, constraints = wald_equal("cC"))

# Model A2a is refitted on exactly the respondents who endorsed at least one
# domain, which is exactly the set left after excluding the 215 zero-endorsers.
# Sensitivity (ii) is therefore NOT an independent analysis -- it IS A2a. The
# identity is asserted here so the report can say so without hand-waving.
stopifnot(identical(sort(rownames(dsub)), sort(rownames(data[data$nsup > 0, ]))))

# =============================================================================
# Model B -- multigroup by family cancer history
# =============================================================================
MB_terms <- paste(sprintf("c(g0_%d, g1_%d)*%s", 2:K, 2:K, CLSD), collapse = " + ")
COV_B <- sub("fcancer21 \\+ ", "", COV)
MB <- paste0(MEAS, "\nWB ~ ", MB_terms, " + ", COV_B, "\n",
             defs("g0_"), defs("g1_"))
fitB <- sem(MB, data, group = "fcancer2", estimator = "MLR", missing = "ML")
rowsB0 <- struct_rows(fitB, group = 1)   # group 1 = fcancer2 == 0, no history
rowsB1 <- struct_rows(fitB, group = 2)   # group 2 = fcancer2 == 1, history
waldB_eq <- lavTestWald(fitB, constraints = paste(
  sprintf("g0_%d == g1_%d", 2:K, 2:K), collapse = "\n"))

# ---- instruction 6: latent variance of WB in each group ---------------------
# If the unstandardised class coefficients are larger in the cancer group only
# because WB is more dispersed there, the standardised coefficients will not
# show the same gap. Both are reported next to SD(WB).
lvcov <- lavInspect(fitB, "cov.lv")
wb_var_total <- sapply(lvcov, function(z) z["WB", "WB"])
peB <- parameterEstimates(fitB)
wb_var_res <- sapply(1:2, function(g)
  peB$est[peB$op == "~~" & peB$lhs == "WB" & peB$rhs == "WB" & peB$group == g][1])
wb_sd_total <- sqrt(wb_var_total)

# Is the larger coefficient in the cancer group specific to the value classes,
# or does every predictor come out larger in that group? If the latter, the
# group difference says something about the group, not about value classes.
CMP <- c(CLSD, "income_z", "health22", "marriage2", "age_z", "disaster1")
grp_cmp <- do.call(rbind, lapply(CMP, function(tm) {
  a <- rowsB0[rowsB0$rhs == tm, ]; b <- rowsB1[rowsB1$rhs == tm, ]
  if (!nrow(a) || !nrow(b)) return(NULL)
  data.frame(term = relab(tm),
             b_nocancer = round(a$est, 4), b_cancer = round(b$est, 4),
             b_ratio = round(b$est / a$est, 2),
             beta_nocancer = round(a$std.all, 3), beta_cancer = round(b$std.all, 3),
             beta_ratio = round(b$std.all / a$std.all, 2),
             p_nocancer = signif(a$pvalue, 3), p_cancer = signif(b$pvalue, 3))
}))
write.csv(grp_cmp, "outputs/table8_groupratio_by_predictor.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# Pre-specified interaction: value class x disaster recovery status
# =============================================================================
for (k in 2:K) for (dd in 1:2) {
  data[[sprintf("x_c%d_d%d", k, dd)]] <-
    data[[paste0("vclass", k)]] * data[[paste0("disaster", dd)]]
}
ix_names <- as.vector(outer(2:K, 1:2, function(k, dd) sprintf("x_c%d_d%d", k, dd)))
ix_terms <- paste(sprintf("i%s*%s", gsub("x_", "", ix_names), ix_names), collapse = " + ")
MI <- paste0(MEAS, "\nWB ~ ", cls_terms("cI"), " + ", COV, " + ", ix_terms, "\n")
fitI <- sem(MI, data, estimator = "MLR", missing = "ML")
waldI <- lavTestWald(fitI, constraints = paste(
  sprintf("i%s == 0", gsub("x_", "", ix_names)), collapse = "\n"))
cells <- table(data$vclass, factor(data$disaster, 0:2,
                                   c("unexposed", "unrecovered", "recovered")))

# =============================================================================
# Breadth-matched contrast: class 2 vs class 3 among respondents who endorsed
# the same number of domains (nsup in {3, 4})
# =============================================================================
# This is the direct test of the compositional hypothesis. If classes 2 and 3
# differ only because one is broader, holding breadth fixed should remove the
# difference; if composition matters, it should survive.
#
# Cell counts drive the design. Within nsup in {3,4}: class 4 has n = 2, which
# cannot support a dummy at all, and class 1 has n = 36. The primary fit is
# therefore restricted to classes 2 and 3 (a single dummy = the contrast
# itself); a secondary fit adds class 1 to show the exclusion is not doing work.
mnsub <- data[data$nsup %in% c(3, 4), ]
cell_n <- table(class = mnsub$vclass, nsup = mnsub$nsup)

# a covariate that is constant in a subsample cannot be estimated; drop and say so
COV_TERMS <- trimws(strsplit(COV, "\\+")[[1]])
usable_cov <- function(z) {
  keep <- sapply(COV_TERMS, function(v) length(unique(z[[v]])) > 1)
  list(str = paste(COV_TERMS[keep], collapse = " + "), dropped = COV_TERMS[!keep])
}

d23 <- mnsub[mnsub$vclass %in% c(2, 3), ]
d23$c3 <- as.numeric(d23$vclass == 3)          # 1 = class 3, 0 = class 2 (ref)
cv23 <- usable_cov(d23)
M_m_A   <- paste0(MEAS, "\nWB ~ mA*c3 + ", cv23$str, "\n")
M_m_A2b <- paste0(MEAS, "\nWB ~ mB*c3 + ", cv23$str, " + mcom_z\n")
fit_mA   <- sem(M_m_A,   d23, estimator = "MLR", missing = "ML")
fit_mA2b <- sem(M_m_A2b, d23, estimator = "MLR", missing = "ML")
row_mA   <- struct_rows(fit_mA);   row_mA   <- row_mA[row_mA$rhs == "c3", ]
row_mA2b <- struct_rows(fit_mA2b); row_mA2b <- row_mA2b[row_mA2b$rhs == "c3", ]

# secondary: keep class 1 (n = 36), class 2 still the reference; class 4 (n = 2)
# is dropped in both fits because two observations cannot identify a dummy
d123 <- mnsub[mnsub$vclass %in% c(1, 2, 3), ]
d123$c1 <- as.numeric(d123$vclass == 1)
d123$c3 <- as.numeric(d123$vclass == 3)
cv123 <- usable_cov(d123)
M_m123 <- paste0(MEAS, "\nWB ~ m1*c1 + m3*c3 + ", cv123$str, "\n")
fit_m123 <- sem(M_m123, d123, estimator = "MLR", missing = "ML")
row_m123 <- struct_rows(fit_m123)

match_tab <- rbind(
  data.frame(model = "matched, Model A covariates (classes 2+3)", n = nrow(d23),
             b = row_mA$est, se = row_mA$se, ci_low = row_mA$ci.lower,
             ci_high = row_mA$ci.upper, beta = row_mA$std.all, p = row_mA$pvalue),
  data.frame(model = "matched, A2b covariates + commitment (classes 2+3)", n = nrow(d23),
             b = row_mA2b$est, se = row_mA2b$se, ci_low = row_mA2b$ci.lower,
             ci_high = row_mA2b$ci.upper, beta = row_mA2b$std.all, p = row_mA2b$pvalue),
  data.frame(model = "matched, Model A covariates (classes 1+2+3)", n = nrow(d123),
             b = row_m123$est[row_m123$rhs == "c3"],
             se = row_m123$se[row_m123$rhs == "c3"],
             ci_low = row_m123$ci.lower[row_m123$rhs == "c3"],
             ci_high = row_m123$ci.upper[row_m123$rhs == "c3"],
             beta = row_m123$std.all[row_m123$rhs == "c3"],
             p = row_m123$pvalue[row_m123$rhs == "c3"]))
write.csv(match_tab, "outputs/table10_breadth_matched.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# nsup-adjusted models, and whether the adjustment is even estimable
# =============================================================================
data$nsup_z <- as.numeric(scale(data$nsup))
dsub$nsup_z <- as.numeric(scale(dsub$nsup))
MA_ns   <- paste0(MEAS, "\nWB ~ ", cls_terms("cN"), " + ", COV, " + nsup_z\n")
MA2b_ns <- paste0(MEAS, "\nWB ~ ", cls_terms("cM"), " + ", COV, " + mcom_z + nsup_z\n")
fitA_ns   <- sem(MA_ns,   data,  estimator = "MLR", missing = "ML")
fitA2b_ns <- sem(MA2b_ns, dsub,  estimator = "MLR", missing = "ML")
waldA_ns   <- lavTestWald(fitA_ns,   constraints = wald_zero("cN"))
waldA2b_ns <- lavTestWald(fitA2b_ns, constraints = wald_zero("cM"))
rowsA_ns   <- struct_rows(fitA_ns)
rowsA2b_ns <- struct_rows(fitA2b_ns)

# THE DECIDING TEST. If, once breadth is adjusted for, the three non-reference
# class effects are statistically indistinguishable, then the class structure
# adds nothing beyond a two-level "minimal vs structured" indicator, and the
# compositional claim has no support left to stand on.
waldA_ns_eq   <- lavTestWald(fitA_ns,   constraints = wald_equal("cN"))
waldA2b_ns_eq <- lavTestWald(fitA2b_ns, constraints = wald_equal("cM"))

# =============================================================================
# Breadth itself as the exposure: nsup as a 7-level categorical variable
# =============================================================================
# Categories 0,1,2,3,4,5,6+ (6 and above pooled: 6..9 hold 127/113/48/10). If
# well-being rises with breadth and then flattens, that saturation is visible
# here without any class model at all.
mk_nscat <- function(z) {
  z$nsupcat <- factor(pmin(z$nsup, 6), levels = 0:6,
                      labels = c("0", "1", "2", "3", "4", "5", "6+"))
  for (lv in 1:6) z[[paste0("nsc", lv)]] <- as.numeric(pmin(z$nsup, 6) == lv)
  z
}
data <- mk_nscat(data); dsub <- mk_nscat(dsub)
NSC_ALL <- paste0("nsc", 1:6)     # reference = 0 domains
NSC_SUB <- paste0("nsc", 2:6)     # in the n=1964 subsample the reference is 1
nscat_n <- table(data$nsupcat)
nscat_n_grp <- table(data$nsupcat, data$fcancer2)   # within-group counts

nsc_terms <- function(prefix, v) paste(sprintf("%s%s*%s", prefix, sub("nsc", "", v), v),
                                       collapse = " + ")
M_ncat    <- paste0(MEAS, "\nWB ~ ", nsc_terms("qA", NSC_ALL), " + ", COV, "\n")
M_ncat_c  <- paste0(MEAS, "\nWB ~ ", nsc_terms("qB", NSC_SUB), " + ", COV, " + mcom_z\n")
M_ncat_k6 <- paste0(MEAS_FREE, "\nWB ~ ", nsc_terms("qK", NSC_ALL), " + ", COV, " + k6_z\n")
M_ncat_g  <- paste0(MEAS, "\nWB ~ ",
                    paste(sprintf("c(h0_%s, h1_%s)*%s", sub("nsc", "", NSC_ALL),
                                  sub("nsc", "", NSC_ALL), NSC_ALL), collapse = " + "),
                    " + ", COV_B, "\n")
fit_ncat    <- sem(M_ncat,    data, estimator = "MLR", missing = "ML")
fit_ncat_c  <- sem(M_ncat_c,  dsub, estimator = "MLR", missing = "ML")
fit_ncat_k6 <- sem(M_ncat_k6, data, estimator = "MLR", missing = "ML")
fit_ncat_g  <- sem(M_ncat_g,  data, group = "fcancer2", estimator = "MLR", missing = "ML")

ncat_rows <- function(fit, group = NULL, ref) {
  r <- struct_rows(fit, group = group)
  r <- r[r$rhs %in% paste0("nsc", 1:6), ]
  r$level <- as.integer(sub("nsc", "", r$rhs))
  r$ref <- ref
  r
}
ncat <- rbind(
  cbind(model = "Model A covariates (ref = 0 domains)", ncat_rows(fit_ncat, ref = 0)),
  cbind(model = "A2b + commitment (ref = 1 domain)", ncat_rows(fit_ncat_c, ref = 1)),
  cbind(model = "K6 added (ref = 0 domains)", ncat_rows(fit_ncat_k6, ref = 0)),
  cbind(model = "Model B: no family cancer (ref = 0)", ncat_rows(fit_ncat_g, group = 1, ref = 0)),
  cbind(model = "Model B: family cancer (ref = 0)", ncat_rows(fit_ncat_g, group = 2, ref = 0)))
ncat$n_level <- as.integer(nscat_n[ncat$level + 1])
ncat_out <- data.frame(model = ncat$model, level = c("1","2","3","4","5","6+")[ncat$level],
                       n = ncat$n_level, b = round(ncat$est, 4),
                       ci_low = round(ncat$ci.lower, 4), ci_high = round(ncat$ci.upper, 4),
                       beta = round(ncat$std.all, 3), p = signif(ncat$pvalue, 4))
write.csv(ncat_out, "outputs/table12_nsup_categorical.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
wald_ncat <- lavTestWald(fit_ncat, constraints = paste(
  sprintf("qA%d == 0", 1:6), collapse = "\n"))
# does the curve flatten? test that the top three levels are equal
wald_ncat_flat <- lavTestWald(fit_ncat, constraints = "qA4 == qA5\nqA5 == qA6")

# nsup only, no class dummies -- the comparison the composition claim must beat
MA_nsonly <- paste0(MEAS, "\nWB ~ nsup_z + ", COV, "\n")
fitA_nsonly <- sem(MA_nsonly, data, estimator = "MLR", missing = "ML")

ic_mods <- list("classes only (Model A)" = fitA,
                "nsup continuous only" = fitA_nsonly,
                "nsup categorical only (7 levels)" = fit_ncat,
                "classes + nsup continuous" = fitA_ns)
ic_tab <- data.frame(
  model = names(ic_mods),
  npar  = sapply(ic_mods, function(f) fitMeasures(f, "npar")),
  AIC   = sapply(ic_mods, function(f) fitMeasures(f, "aic")),
  BIC   = sapply(ic_mods, function(f) fitMeasures(f, "bic")))
rownames(ic_tab) <- NULL
ic_tab$dAIC <- round(ic_tab$AIC - min(ic_tab$AIC), 2)
ic_tab$dBIC <- round(ic_tab$BIC - min(ic_tab$BIC), 2)
write.csv(ic_tab, "outputs/table11_nsup_ic.csv", row.names = FALSE, fileEncoding = "UTF-8")

# nsup distribution within class
nsup_dist <- do.call(rbind, lapply(seq_len(K), function(k) {
  x <- data$nsup[data$vclass == k]
  q <- quantile(x, c(0, .25, .5, .75, 1))
  data.frame(class = k, label = CLASS_LAB[k], n = length(x), min = q[1], Q1 = q[2],
             median = q[3], Q3 = q[4], max = q[5], mean = round(mean(x), 2))
}))
rownames(nsup_dist) <- NULL

# how much of nsup is already determined by class membership
r2_nsup_class <- summary(lm(nsup ~ vclass, data = data))$r.squared

# Overlap of the classes on nsup. The nsup-adjusted class-4 coefficient compares
# class 4 with the other classes at the same breadth; where almost nobody else
# has that breadth, the comparison is extrapolation rather than contrast.
ov <- table(class = data$vclass, nsup = factor(data$nsup, levels = 0:9))
p3 <- ov[3, ] / sum(ov[3, ]); p4 <- ov[4, ] / sum(ov[4, ])
overlap_34 <- sum(pmin(p3, p4))                  # overlap coefficient, classes 3 and 4
hi <- as.character(6:9)                          # where class 4 mass sits
share_c4_hi <- sum(ov[4, hi]) / sum(ov[, hi])    # class 4's share of nsup >= 6
n_hi_other  <- sum(ov[-4, hi])                   # everyone else at nsup >= 6
ovl_tab <- as.data.frame.matrix(ov)
ovl_tab$total <- rowSums(ov)
write.csv(cbind(class = rownames(ovl_tab), label = CLASS_LAB, ovl_tab),
          "outputs/table13_class_nsup_overlap.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
wb_score_ns <- lavPredict(fitA_ns)[, "WB"]
vif_ns <- vif_pub(lm(wb_score_ns ~ vclass + nsup_z + disaster + fcancer2 + female +
                       age_z + income_z + edu2 + occupation2 + marriage +
                       comobidity + fukushima + health2, data = data))

# =============================================================================
# Instruction 5: K6 in / K6 out, every predictor
# =============================================================================
# The K6 model keeps the ResAll disturbance FREE -- see the header. Everything
# else about the two specifications is identical.
MK6 <- paste0(MEAS_FREE, "\nWB ~ ", cls_terms("cK"), " + ", COV, " + k6_z\n", defs("cK"))
fitK6 <- sem(MK6, data, estimator = "MLR", missing = "ML")
rowsK6 <- struct_rows(fitK6)
waldK6 <- lavTestWald(fitK6, constraints = wald_zero("cK"))
waldK6_eq <- lavTestWald(fitK6, constraints = wald_equal("cK"))
fitK6_meas <- fitMeasures(fitK6, c("chisq", "df", "cfi.robust", "rmsea.robust", "srmr"))
peK6 <- parameterEstimates(fitK6)
res_k6 <- peK6[peK6$op == "~~" & peK6$lhs == "ResAll" & peK6$rhs == "ResAll", ]

k6_tab <- merge(
  rowsA[, c("rhs", "term", "est", "ci.lower", "ci.upper", "std.all", "pvalue")],
  rowsK6[, c("rhs", "est", "ci.lower", "ci.upper", "std.all", "pvalue")],
  by = "rhs", suffixes = c("_noK6", "_K6"), all = TRUE)
k6_tab <- k6_tab[match(c(CLSD, "disaster1", "disaster2", "fcancer21", "female1",
                         "age_z", "income_z", "edu22", "edu23", "occupation22",
                         "occupation23", "marriage2", "marriage3", "comobidity1",
                         "fukushima1", "health22", "k6_z"), k6_tab$rhs), ]
k6_tab <- k6_tab[!is.na(k6_tab$rhs), ]
k6_tab$term <- relab(k6_tab$rhs)
k6_tab$d_beta <- round(k6_tab$std.all_K6 - k6_tab$std.all_noK6, 4)
# `term` is the item name used in every table and in the manuscript. The raw
# lavaan right-hand-side name is kept only so a number can be traced back to the
# model syntax, and is renamed so it cannot be mistaken for a display label --
# "disaster1"/"disaster2" must never read as an item name anywhere.
names(k6_tab)[names(k6_tab) == "rhs"] <- "lavaan_rhs_raw"
k6_tab <- k6_tab[, c("term", setdiff(names(k6_tab), c("term", "lavaan_rhs_raw")),
                     "lavaan_rhs_raw")]
write.csv(k6_tab, "outputs/table5_k6_contrast.csv", row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# Instruction 4b: class effects on each FIRST-ORDER factor
# =============================================================================
# No second-order factor. NFE and ResAll are each regressed on the same
# predictors and their disturbances are allowed to covary. This shows whether a
# class acts on negative-free emotion, on resource-oriented well-being, or both.
MFO <- paste("NFE =~ sadR + worryR",
             "ResAll =~ ls + gh + es + vital + comp",
             "ls ~~ gh", "es ~~ vital",
             paste0("NFE ~ ", paste(sprintf("n%d*%s", 2:K, CLSD), collapse = " + "),
                    " + ", COV),
             paste0("ResAll ~ ", paste(sprintf("r%d*%s", 2:K, CLSD), collapse = " + "),
                    " + ", COV),
             "NFE ~~ ResAll", sep = "\n")
fitFO <- sem(MFO, data, estimator = "MLR", missing = "ML")

# The point of this model is that it does NOT route the class effects through a
# second-order factor. If it did -- or if the ResAll disturbance were fixed at
# zero here -- the path to NFE would be forced through WB and a null on NFE
# would be an artefact of the parameterisation rather than a finding. Both
# conditions are asserted rather than assumed.
stopifnot(setequal(lavNames(fitFO, "lv"), c("NFE", "ResAll")))       # no WB
peFO <- parameterEstimates(fitFO)
stopifnot(sum(peFO$lhs == "WB") == 0)                                 # no WB anywhere
fo_dist <- peFO[peFO$op == "~~" & peFO$lhs == "ResAll" & peFO$rhs == "ResAll", ]
fo_cov  <- peFO[peFO$op == "~~" & peFO$lhs == "NFE" & peFO$rhs == "ResAll", ]
stopifnot(fo_dist$se[1] > 0, fo_dist$est[1] > 0)   # ResAll disturbance FREE, > 0
stopifnot(fo_cov$se[1] > 0)                        # disturbances free to covary

rowsFO_n <- struct_rows(fitFO, lhs = "NFE")
rowsFO_r <- struct_rows(fitFO, lhs = "ResAll")

# ---------------------------------------------------------------------------
# Same decomposition with BREADTH as the exposure, which is what the paper
# actually reports. The class version above is retained for the supplement.
#
# SCALE WARNING. The two first-order factors are on different metrics: NFE is
# marked by a binary item and ResAll by a 0-10 item, so their unstandardised
# coefficients are NOT comparable with each other, and neither is comparable
# with the second-order coefficients in Table 3. Any claim that an interval on
# NFE "excludes" an effect seen on ResAll has to be made in standardised units.
# Both are therefore reported for every predictor.
MFO_ns <- paste("NFE =~ sadR + worryR",
                "ResAll =~ ls + gh + es + vital + comp",
                "ls ~~ gh", "es ~~ vital",
                paste0("NFE ~ ",
                       paste(sprintf("fn%d*nsc%d", 1:6, 1:6), collapse = " + "),
                       " + ", COV),
                paste0("ResAll ~ ",
                       paste(sprintf("fr%d*nsc%d", 1:6, 1:6), collapse = " + "),
                       " + ", COV),
                "NFE ~~ ResAll", sep = "\n")
fitFO_ns <- sem(MFO_ns, data, estimator = "MLR", missing = "ML")
stopifnot(setequal(lavNames(fitFO_ns, "lv"), c("NFE", "ResAll")))
rowsFOns_n <- struct_rows(fitFO_ns, lhs = "NFE")
rowsFOns_r <- struct_rows(fitFO_ns, lhs = "ResAll")
waldFOns_n <- lavTestWald(fitFO_ns, constraints = paste(
  sprintf("fn%d == 0", 1:6), collapse = "\n"))
waldFOns_r <- lavTestWald(fitFO_ns, constraints = paste(
  sprintf("fr%d == 0", 1:6), collapse = "\n"))
fitFOns_meas <- fitMeasures(fitFO_ns, c("chisq", "df", "cfi.robust",
                                        "rmsea.robust", "srmr"))

# Standardised solution with its own confidence intervals. The comparison that
# decides the state-dependent-responding argument is between a NFE interval and
# the ResAll effect sizes, and that comparison is only meaningful in
# standardised units, so the interval must be the standardised one rather than
# an unstandardised interval rescaled by a point ratio.
std_rows <- function(fit, lhs) {
  ss <- standardizedSolution(fit)
  ss <- ss[ss$op == "~" & ss$lhs == lhs, ]
  ss$term <- relab(ss$rhs)
  ss
}
stdFOns_n <- std_rows(fitFO_ns, "NFE")
stdFOns_r <- std_rows(fitFO_ns, "ResAll")
nfe_upper <- max(stdFOns_n$ci.upper[stdFOns_n$rhs %in% paste0("nsc", 1:6)])
res_std   <- stdFOns_r$est.std[stdFOns_r$rhs %in% paste0("nsc", 1:6)]

fo_full <- function(rows, factor_lab, model_lab) data.frame(
  model = model_lab, factor = factor_lab, term = rows$term,
  b = round(rows$est, 4), se = round(rows$se, 4),
  ci_low = round(rows$ci.lower, 4), ci_high = round(rows$ci.upper, 4),
  beta = round(rows$std.all, 4), p = signif(rows$pvalue, 4))
fo_full_tab <- rbind(
  fo_full(rowsFOns_n, "NFE (negative-free emotion)", "breadth as exposure"),
  fo_full(rowsFOns_r, "ResAll (resource-oriented well-being)", "breadth as exposure"),
  fo_full(rowsFO_n, "NFE (negative-free emotion)", "class as exposure"),
  fo_full(rowsFO_r, "ResAll (resource-oriented well-being)", "class as exposure"))
write.csv(fo_full_tab, "outputs/table14_firstorder_full.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# Reliability of the two first-order factors, recomputed on this dataset so the
# values transcribed into Methods from paper 1 are verified rather than copied.
# NFE has two binary indicators, so the standardised alpha equals the
# Spearman-Brown prophecy for two items; both are reported, with the mean
# inter-item correlation, because a two-item binary factor is a real limitation.
NFE_items <- data[, c("sadR", "worryR")]
RES_items <- data[, c("ls", "gh", "es", "vital", "comp")]
rel_nfe_alpha <- psych::alpha(NFE_items, warnings = FALSE)$total$std.alpha
rel_nfe_r     <- cor(NFE_items$sadR, NFE_items$worryR, use = "complete.obs")
rel_nfe_sb    <- 2 * rel_nfe_r / (1 + rel_nfe_r)
rel_res_alpha <- psych::alpha(RES_items, warnings = FALSE)$total$raw_alpha
rel_res_omega <- suppressWarnings(
  psych::omega(RES_items, nfactors = 1, plot = FALSE)$omega.tot)
rel_tab <- data.frame(
  factor = c("NFE (negative-free emotion)", "NFE", "NFE",
             "ResAll (resource-oriented well-being)", "ResAll"),
  statistic = c("standardised alpha", "mean inter-item correlation",
                "Spearman-Brown coefficient", "Cronbach alpha", "McDonald omega"),
  value = round(c(rel_nfe_alpha, rel_nfe_r, rel_nfe_sb,
                  rel_res_alpha, rel_res_omega), 3),
  paper1_published = c(0.650, 0.480, 0.650, 0.780, 0.840))
# Paper 1 publishes these to two decimals; comparing an unrounded value with a
# two-decimal published figure at a four-decimal tolerance fails on display
# precision alone, so the observed value is rounded to the published precision.
rel_tab$matches_paper1 <- round(rel_tab$value, 2) == round(rel_tab$paper1_published, 2)
write.csv(rel_tab, "outputs/table15_reliability.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
fitFO_meas <- fitMeasures(fitFO, c("chisq", "df", "cfi.robust", "rmsea.robust", "srmr"))
waldFO_n <- lavTestWald(fitFO, constraints = paste(sprintf("n%d == 0", 2:K), collapse = "\n"))
waldFO_r <- lavTestWald(fitFO, constraints = paste(sprintf("r%d == 0", 2:K), collapse = "\n"))
fo_tab <- rbind(
  data.frame(factor = "NFE (negative-free emotion)",
             term = rowsFO_n$term, b = round(rowsFO_n$est, 4),
             se = round(rowsFO_n$se, 4), ci_low = round(rowsFO_n$ci.lower, 4),
             ci_high = round(rowsFO_n$ci.upper, 4),
             beta = round(rowsFO_n$std.all, 4), p = signif(rowsFO_n$pvalue, 4)),
  data.frame(factor = "ResAll (resource-oriented well-being)",
             term = rowsFO_r$term, b = round(rowsFO_r$est, 4),
             se = round(rowsFO_r$se, 4), ci_low = round(rowsFO_r$ci.lower, 4),
             ci_high = round(rowsFO_r$ci.upper, 4),
             beta = round(rowsFO_r$std.all, 4), p = signif(rowsFO_r$pvalue, 4)))
write.csv(fo_tab, "outputs/table6_firstorder.csv", row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# VIF -- paper 1's convention (GVIF^(1/(2*Df)) column)
# =============================================================================
wb_score <- lavPredict(fitA)[, "WB"]
vif_main <- vif_pub(lm(wb_score ~ vclass + disaster + fcancer2 + female + age_z +
                         income_z + edu2 + occupation2 + marriage + comobidity +
                         fukushima + health2, data = data))
wb_score_I <- lavPredict(fitI)[, "WB"]
vif_ix <- suppressWarnings(vif_pub(
  lm(wb_score_I ~ vclass * disaster + fcancer2 + female + age_z + income_z +
       edu2 + occupation2 + marriage + comobidity + fukushima + health2,
     data = data)))

# =============================================================================
# Sensitivity analyses (supplementary): (ii), (iii), (iv).
# (i) is no longer here -- the K6 model is now a main-text companion.
# =============================================================================
fitS2 <- fitA2a                                    # (ii) exclude the 215

imp_mat <- as.matrix(data[, paste0("Q11_3_", DOM)])
sup_mat <- as.matrix(data[, IND])
n_low_imp <- sum(sup_mat == 1 & !is.na(imp_mat) & imp_mat <= 2)
ord <- matrix(1L, nrow(data), 9)
ord[sup_mat == 1 & !is.na(imp_mat) & imp_mat <= 4] <- 2L
ord[sup_mat == 1 & !is.na(imp_mat) & imp_mat == 5] <- 3L
ord <- as.data.frame(ord); names(ord) <- IND
fo_formula <- as.formula(paste0("cbind(", paste(IND, collapse = ", "), ") ~ 1"))
ORD_CACHE <- sprintf("lca_ordinal_K%d.rds", K)
ord_key <- list(seed = LC$seed, nrep = LC$nrep, K = K, tab = table(unlist(ord)))
ord_cached <- if (file.exists(ORD_CACHE)) readRDS(ORD_CACHE) else NULL
if (!is.null(ord_cached) && identical(ord_cached$key, ord_key)) {
  lca_ord <- ord_cached$fit
} else {
  set.seed(LC$seed)
  lca_ord <- poLCA(fo_formula, ord, nclass = K, nrep = LC$nrep, maxiter = 50000,
                   tol = 1e-10, na.rm = FALSE, graphs = FALSE, verbose = FALSE)
  saveRDS(list(key = ord_key, fit = lca_ord), ORD_CACHE)
}
o <- order(lca_ord$P, decreasing = TRUE)
lca_ord$posterior <- lca_ord$posterior[, o, drop = FALSE]
lca_ord$predclass <- apply(lca_ord$posterior, 1, which.max)
data$vclass_ord <- factor(lca_ord$predclass, levels = seq_len(K))
data <- cbind(data, setNames(
  as.data.frame(model.matrix(~ vclass_ord, data)[, -1, drop = FALSE]),
  paste0("vclassO", 2:K)))
MS3 <- paste0(MEAS, "\nWB ~ ",
              paste(sprintf("s3_%d*vclassO%d", 2:K, 2:K), collapse = " + "),
              " + ", COV, "\n")
fitS3 <- sem(MS3, data, estimator = "MLR", missing = "ML")
xtab_ord <- table(binary = data$vclass, ordinal = data$vclass_ord)

maxpp <- apply(m$posterior, 1, max)
data$pp_max <- maxpp
keep80 <- maxpp >= 0.80
n_keep80 <- sum(keep80)
keep80_by_class <- table(data$vclass, keep80)
MS4a <- paste0(MEAS, "\nWB ~ ", cls_terms("s4a"), " + ", COV, "\n")
MS4b <- paste0(MEAS, "\nWB ~ ", cls_terms("s4b"), " + ", COV, " + mcom_z\n")
fitS4a <- sem(MS4a, data[keep80, ], estimator = "MLR", missing = "ML")
fitS4b <- sem(MS4b, data[keep80 & data$nsup > 0, ], estimator = "MLR", missing = "ML")
n_S4b <- sum(keep80 & data$nsup > 0)

# supplementary material on class enumeration
K5 <- 5
m5 <- LC$models[[K5]]
P5 <- sapply(IND, function(v) m5$probs[[v]][, 2])
xtab_45 <- table(K4 = m$predclass, K5 = m5$predclass)
cells5_dis <- table(m5$predclass, factor(d$disaster, 0:2,
                                         c("unexposed", "unrecovered", "recovered")))
cells5_fc  <- table(m5$predclass, d$fcancer2)

# =============================================================================
# Table 1 -- background by value class
# =============================================================================
pct <- function(x, lv) sprintf("%d (%.1f)", sum(x == lv), round(100 * mean(x == lv), 1))
msd <- function(x) sprintf("%.1f (%.1f)", round(mean(x), 1), round(sd(x), 1))
t1_rows <- list()
add1 <- function(var, fn) {
  t1_rows[[length(t1_rows) + 1]] <<- data.frame(
    variable = var, overall = fn(data),
    setNames(as.data.frame(t(sapply(seq_len(K), function(k) fn(data[data$vclass == k, ])))),
             paste0("class", seq_len(K))))
}
add1("n", function(z) as.character(nrow(z)))
add1("Age, mean (SD)", function(z) msd(z$age))
add1("Female, n (%)", function(z) pct(z$female == 1, TRUE))
add1("Family cancer history, n (%)", function(z) pct(z$fcancer2 == 1, TRUE))
add1("Fukushima resident, n (%)", function(z) pct(z$fukushima == 1, TRUE))
add1("Disaster: unexposed, n (%)", function(z) pct(z$disaster == 0, TRUE))
add1("Disaster: unrecovered, n (%)", function(z) pct(z$disaster == 1, TRUE))
add1("Disaster: recovered, n (%)", function(z) pct(z$disaster == 2, TRUE))
add1("Income <3M JPY, n (%)", function(z) pct(z$income3 == 1, TRUE))
add1("Income 3-6M JPY, n (%)", function(z) pct(z$income3 == 2, TRUE))
add1("Income >=6M JPY, n (%)", function(z) pct(z$income3 == 3, TRUE))
add1("Education: junior/high, n (%)", function(z) pct(z$edu2 == 1, TRUE))
add1("Education: beyond high, n (%)", function(z) pct(z$edu2 == 2, TRUE))
add1("Education: declined, n (%)", function(z) pct(z$edu2 == 3, TRUE))
add1("Never married, n (%)", function(z) pct(z$marriage == 1, TRUE))
add1("Married, n (%)", function(z) pct(z$marriage == 2, TRUE))
add1("Divorced/widowed, n (%)", function(z) pct(z$marriage == 3, TRUE))
add1("Self-employed, n (%)", function(z) pct(z$occupation2 == 1, TRUE))
add1("Employee, n (%)", function(z) pct(z$occupation2 == 2, TRUE))
add1("Other occupation, n (%)", function(z) pct(z$occupation2 == 3, TRUE))
add1("Subjective health high, n (%)", function(z) pct(z$health2 == 2, TRUE))
add1("Comorbidity, n (%)", function(z) pct(z$comobidity == 1, TRUE))
add1("K6, mean (SD)", function(z) msd(z$k6))
table1 <- do.call(rbind, t1_rows)
write.csv(table1, "outputs/table1_background.csv", row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# Table 2 -- class profile
# =============================================================================
P <- sapply(IND, function(v) m$probs[[v]][, 2])
t2 <- data.frame(class = seq_len(K), label = CLASS_LAB,
                 n = as.integer(table(data$vclass)),
                 pct = round(100 * as.numeric(table(data$vclass)) / nrow(data), 1),
                 class_share_model = round(100 * m$P, 1))
for (j in 1:9) t2[[paste0("p_", DOM[j])]] <- round(P[, j], 3)
t2$mean_n_domains <- round(tapply(data$nsup, data$vclass, mean), 2)
t2$model_expected_n_domains <- round(rowSums(P), 2)
for (v in c("m_importance", "m_commitment", "m_intention", "m_consistency")) {
  t2[[v]] <- round(tapply(data[[v]], data$vclass, mean, na.rm = TRUE), 2)
  t2[[paste0(v, "_n")]] <- tapply(!is.na(data[[v]]), data$vclass, sum)
}
write.csv(t2, "outputs/table2_class_profile.csv", row.names = FALSE, fileEncoding = "UTF-8")

ent_K <- LC$fit$entropy[LC$fit$K == K]
avepp <- sapply(seq_len(K), function(k) mean(m$posterior[m$predclass == k, k]))
worst_k <- which.min(avepp)

# =============================================================================
# Instruction 3: slope per additional endorsed domain
# =============================================================================
# For each adjacent pair of classes, the difference in b divided by the
# difference in mean number of domains endorsed. If well-being tracked breadth
# alone this quantity would be roughly constant across the three steps.
# The CI is the contrast CI rescaled by the same (positive) gap, which is a
# monotone rescaling and not a new estimate.
mnd <- as.numeric(t2$mean_n_domains)
bA  <- c(0, sapply(2:K, function(k) rowsA$est[rowsA$rhs == paste0("vclass", k)]))
step_lab <- sprintf("class%d -> class%d", 1:(K - 1), 2:K)
step_def <- c("cA_21", "cA_32", "cA_43")[1:(K - 1)]
slope <- do.call(rbind, lapply(seq_len(K - 1), function(i) {
  gap <- mnd[i + 1] - mnd[i]
  dd  <- getdef(fitA, step_def[i])
  data.frame(step = step_lab[i],
             mean_domains_from = mnd[i], mean_domains_to = mnd[i + 1],
             domain_gap = round(gap, 2),
             b_from = round(bA[i], 4), b_to = round(bA[i + 1], 4),
             b_difference = round(dd$est, 4),
             ci_low = round(dd$ci.lower, 4), ci_high = round(dd$ci.upper, 4),
             p = signif(dd$pvalue, 4),
             slope_per_domain = round(dd$est / gap, 4),
             slope_ci_low = round(dd$ci.lower / gap, 4),
             slope_ci_high = round(dd$ci.upper / gap, 4))
}))
write.csv(slope, "outputs/table7_slope_per_domain.csv", row.names = FALSE)

# =============================================================================
# Instruction 2: class 2 vs class 3 contrasts
# =============================================================================
pair_tab <- function(suffix) {
  z <- rbind(
    cbind(model = "Model A (n=2179)", getdef(fitA, paste0("cA", suffix))),
    cbind(model = "Model A2b (n=1964, + commitment)", getdef(fitA2b, paste0("cC", suffix))),
    cbind(model = "Model B: no family cancer (n=1324)", getdef(fitB, paste0("g0_", suffix))),
    cbind(model = "Model B: family cancer (n=855)", getdef(fitB, paste0("g1_", suffix))),
    cbind(model = "K6 model (n=2179)", getdef(fitK6, paste0("cK", suffix))))
  z$label <- NULL
  names(z)[names(z) == "est"] <- "b"
  z
}
con23 <- pair_tab("_32")
con24 <- pair_tab("_42")
con34 <- pair_tab("_43")
con_all <- rbind(cbind(contrast = "class2 vs class3", con23),
                 cbind(contrast = "class2 vs class4", con24),
                 cbind(contrast = "class3 vs class4", con34))
write.csv(con_all, "outputs/contrasts_class2_vs_3.csv", row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# Inventory of every hypothesis test reported
# =============================================================================
# No multiplicity adjustment is applied. The inventory exists so that the scope
# of what may be called confirmatory can be drawn deliberately: tests planned
# before looking (the joint class effect, the commitment adjustment, the group
# comparison, the one pre-specified interaction) are separated from tests added
# afterwards (all pairwise class contrasts, the equality tests, the first-order
# decomposition, the per-domain slopes).
inv <- function(family, test, stat, df, p, status)
  data.frame(family = family, test = test, statistic = round(stat, 3),
             df = df, p = signif(p, 4), status = status)
tests <- rbind(
  inv("class effect", "Model A: all class dummies = 0", waldA$stat, waldA$df, waldA$p.value, "pre-specified"),
  inv("class effect", "Model A2a: all class dummies = 0", waldA2a$stat, waldA2a$df, waldA2a$p.value, "pre-specified"),
  inv("class effect", "Model A2b: all class dummies = 0", waldA2b$stat, waldA2b$df, waldA2b$p.value, "pre-specified"),
  inv("class effect", "K6 model: all class dummies = 0", waldK6$stat, waldK6$df, waldK6$p.value, "pre-specified"),
  inv("gradient", "Model A: b2 = b3 = b4", waldA_eq$stat, waldA_eq$df, waldA_eq$p.value, "added 2026-08-30"),
  inv("gradient", "Model A2b: b2 = b3 = b4", waldA2b_eq$stat, waldA2b_eq$df, waldA2b_eq$p.value, "added 2026-08-30"),
  inv("gradient", "K6 model: b2 = b3 = b4", waldK6_eq$stat, waldK6_eq$df, waldK6_eq$p.value, "added 2026-08-30"),
  inv("breadth-matched", "matched nsup 3-4, Model A covariates: class3 vs class2",
      row_mA$z^2, 1, row_mA$pvalue, "added 2026-08-30"),
  inv("breadth-matched", "matched nsup 3-4, + commitment: class3 vs class2",
      row_mA2b$z^2, 1, row_mA2b$pvalue, "added 2026-08-30"),
  inv("breadth-matched", "matched nsup 3-4, classes 1+2+3: class3 vs class2",
      row_m123$z[row_m123$rhs == "c3"]^2, 1,
      row_m123$pvalue[row_m123$rhs == "c3"], "added 2026-08-30"),
  inv("nsup-adjusted", "Model A + nsup: all class dummies = 0",
      waldA_ns$stat, waldA_ns$df, waldA_ns$p.value, "added 2026-08-30"),
  inv("nsup-adjusted", "Model A2b + nsup: all class dummies = 0",
      waldA2b_ns$stat, waldA2b_ns$df, waldA2b_ns$p.value, "added 2026-08-30"),
  inv("nsup-adjusted", "Model A + nsup: b2 = b3 = b4 (deciding test)",
      waldA_ns_eq$stat, waldA_ns_eq$df, waldA_ns_eq$p.value, "added 2026-08-30"),
  inv("nsup-adjusted", "Model A2b + nsup: b2 = b3 = b4 (deciding test)",
      waldA2b_ns_eq$stat, waldA2b_ns_eq$df, waldA2b_ns_eq$p.value, "added 2026-08-30"),
  inv("breadth as exposure", "nsup categorical: all levels = 0",
      wald_ncat$stat, wald_ncat$df, wald_ncat$p.value, "added 2026-08-30"),
  inv("breadth as exposure", "nsup categorical: levels 4 = 5 = 6+ (saturation)",
      wald_ncat_flat$stat, wald_ncat_flat$df, wald_ncat_flat$p.value, "added 2026-08-30"),
  inv("group", "Model B: class coefficients equal across groups", waldB_eq$stat, waldB_eq$df, waldB_eq$p.value, "pre-specified"),
  inv("interaction", "class x recovery, all interaction terms = 0", waldI$stat, waldI$df, waldI$p.value, "pre-specified"),
  inv("first-order", "class effects on NFE = 0", waldFO_n$stat, waldFO_n$df, waldFO_n$p.value, "added 2026-08-30"),
  inv("first-order", "class effects on ResAll = 0", waldFO_r$stat, waldFO_r$df, waldFO_r$p.value, "added 2026-08-30"))
# individual class-vs-reference coefficients
add_coefs <- function(rows, model, status, family = "class vs reference") {
  r <- rows[rows$rhs %in% CLSD, ]
  if (!nrow(r)) return(NULL)
  inv(family, sprintf("%s: %s", model, r$term), r$z^2, 1, r$pvalue, status)
}
tests <- rbind(tests,
  add_coefs(rowsA, "Model A", "pre-specified"),
  add_coefs(rowsA2a, "Model A2a", "pre-specified"),
  add_coefs(rowsA2b, "Model A2b", "pre-specified"),
  add_coefs(rowsK6, "K6 model", "pre-specified"),
  add_coefs(rowsB0, "Model B / no family cancer", "pre-specified"),
  add_coefs(rowsB1, "Model B / family cancer", "pre-specified"))
# pairwise class contrasts
tests <- rbind(tests, do.call(rbind, lapply(seq_len(nrow(con_all)), function(i) {
  r <- con_all[i, ]
  inv("pairwise class contrast", sprintf("%s: %s", r$model, r$contrast),
      r$z^2, 1, r$pvalue, "added 2026-08-30")
})))
# individual interaction terms
peI <- parameterEstimates(fitI)
peI <- peI[peI$op == "~" & peI$lhs == "WB" & grepl("^x_c", peI$rhs), ]
tests <- rbind(tests, inv("interaction", paste("individual:", peI$rhs),
                          peI$z^2, 1, peI$pvalue, "pre-specified"))
# individual first-order class coefficients
tests <- rbind(tests,
  add_coefs(rowsFO_n, "NFE", "added 2026-08-30", "first-order individual"),
  add_coefs(rowsFO_r, "ResAll", "added 2026-08-30", "first-order individual"))
# The per-domain slope table performs NO new test. Each of its three rows is one
# of the contrasts already counted above (cA_21 = class2 vs reference in Model A,
# cA_32 and cA_43 = the Model A pairwise contrasts), divided by a constant. It is
# therefore deliberately NOT added to the inventory, which would double-count it.
slope_dup <- c("Model A: class2 vs reference", "Model A pairwise class2 vs class3",
               "Model A pairwise class3 vs class4")
rownames(tests) <- NULL
write.csv(tests, "outputs/table9_test_inventory.csv", row.names = FALSE, fileEncoding = "UTF-8")
n_sig <- sum(tests$p < 0.05, na.rm = TRUE)

# =============================================================================
# Figure 1
# =============================================================================
pal <- c("#1f4e79", "#c00000", "#2e7d32", "#e69f00", "#7b3294", "#0097a7")
png("outputs/fig1_class_profiles.png", width = 1650, height = 1000, res = 150)
par(mar = c(7.5, 4.6, 3, 1))
plot(NA, xlim = c(1, 9), ylim = c(0, 1), xaxt = "n", xlab = "",
     ylab = "Conditional probability of endorsement",
     main = sprintf("Value-domain classes (%d-class solution, N = %d)", K, nrow(data)))
axis(1, at = 1:9, labels = DOM_LAB, las = 2)
abline(h = seq(0, 1, 0.2), col = "grey90")
marg <- colMeans(d[, IND])
lines(1:9, marg, col = "grey45", lwd = 2, lty = 3, type = "b", pch = 4)
for (k in seq_len(K)) lines(1:9, P[k, ], col = pal[k], lwd = 2.4, type = "b", pch = 19)
legend("topright", bty = "n", cex = 0.82,
       legend = c(sprintf("%s (n = %d, %.1f%%, mean %.2f domains)",
                          CLASS_LAB, t2$n, t2$pct, t2$mean_n_domains),
                  "Overall endorsement rate"),
       col = c(pal[seq_len(K)], "grey45"), lwd = 2,
       lty = c(rep(1, K), 3), pch = c(rep(19, K), 4))
dev.off()

# Figure 2 -- breadth as the exposure, estimate and interval at every level
png("outputs/fig2_nsup_categorical.png", width = 1650, height = 1000, res = 150)
par(mar = c(5, 4.8, 3.2, 1))
mods <- unique(ncat$model)
off <- seq(-0.18, 0.18, length.out = length(mods))
yl <- range(c(ncat$ci.lower, ncat$ci.upper, 0)) + c(-0.01, 0.01)
plot(NA, xlim = c(0.6, 6.4), ylim = yl, xaxt = "n",
     xlab = "Number of value domains endorsed", ylab = "b on second-order well-being",
     main = "Well-being by breadth of endorsement (reference in parentheses)")
axis(1, at = 1:6, labels = c("1", "2", "3", "4", "5", "6+"))
abline(h = 0, col = "grey55", lty = 2)
for (i in seq_along(mods)) {
  z <- ncat[ncat$model == mods[i], ]
  x <- z$level + off[i]
  segments(x, z$ci.lower, x, z$ci.upper, col = pal[i], lwd = 2)
  points(x, z$est, col = pal[i], pch = 19, cex = 1.05)
  lines(x, z$est, col = pal[i], lwd = 1.4)
}
legend("topleft", bty = "n", cex = 0.78, legend = mods, col = pal[seq_along(mods)],
       lwd = 2, pch = 19)
dev.off()

# =============================================================================
# Tables 3 and 4 as CSV
# =============================================================================
tab_est <- function(rows, model) data.frame(
  model = model, term = rows$term, b = round(rows$est, 4), se = round(rows$se, 4),
  ci_low = round(rows$ci.lower, 4), ci_high = round(rows$ci.upper, 4),
  beta = round(rows$std.all, 4), p = signif(rows$pvalue, 4))
table3 <- rbind(tab_est(rowsA, "A (n=2179)"),
                tab_est(rowsA2a, "A2a (n=1964, no commitment)"),
                tab_est(rowsA2b, "A2b (n=1964, + commitment)"))
write.csv(table3, "outputs/table3_modelA.csv", row.names = FALSE, fileEncoding = "UTF-8")
table4 <- rbind(tab_est(rowsB0, "B: no family cancer (n=1324)"),
                tab_est(rowsB1, "B: family cancer (n=855)"))
write.csv(table4, "outputs/table4_modelB.csv", row.names = FALSE, fileEncoding = "UTF-8")
sens <- rbind(tab_est(struct_rows(fitS2), "S(ii) exclude 215"),
              tab_est(struct_rows(fitS3), "S(iii) ordinal-LCA classes"),
              tab_est(struct_rows(fitS4a), sprintf("S(iv) ModelA, pp>=0.80 (n=%d)", n_keep80)),
              tab_est(struct_rows(fitS4b), sprintf("S(iv) ModelA2, pp>=0.80 (n=%d)", n_S4b)))
write.csv(sens, "outputs/sensitivity.csv", row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# Report
# =============================================================================
con <- file("outputs/results_phase2.md", open = "wt", encoding = "UTF-8")
w <- function(...) cat(sprintf(...), file = con)
cls_row <- function(rows, k) rows[rows$rhs == paste0("vclass", k), ]

w("# Phase 2 \u7d50\u679c \u2014 \u4fa1\u5024\u30af\u30e9\u30b9 \u2192 well-being\n\n")
w("\u751f\u6210: `02_structural.R %d`\uff08%s\uff09\u3002K = %d\u3001\u53c2\u7167\u306f\u30af\u30e9\u30b9 1\uff08\u6700\u5927\u30af\u30e9\u30b9\uff09\u3002\n\n",
  K, format(t_start, "%Y-%m-%d %H:%M"), K)
w("\u30af\u30e9\u30b9\u540d: %s\n\n", paste(sprintf("%d = %s", seq_len(K), CLASS_LAB), collapse = ", "))

w("## 0. \u5909\u6570\u540d\u306e\u78ba\u5b9a \u2014 disaster1 / disaster2\n\n")
w("`00_derive.R` 151-153 \u884c:\n\n")
w("```r\nd_ex2 <- ifelse(q(\"Q3_1\") %%in%% c(1, 2), 0, 1)\n")
w("d_re2 <- ifelse(q(\"Q3_2\") %%in%% c(1, 2), 1, ifelse(q(\"Q3_2\") %%in%% c(3, 4), 0, NA))\n")
w("d$disaster <- ifelse(d_ex2 == 0, 0, ifelse(d_re2 == 1, 2, 1))\n```\n\n")
w("`d_re2 == 1`\uff08Q3_2 \u2208 {1,2}\u3001\u56de\u5fa9\u6e08\u307f\u30fb\u898b\u901a\u3057\u3042\u308a\uff09\u306e\u3068\u304d **2**\u3001")
w("\u305d\u308c\u4ee5\u5916\u306e\u88ab\u5bb3\u3042\u308a\uff08Q3_2 \u2208 {3,4}\uff09\u304c **1** \u3067\u3042\u308b\u3002")
w("`factor(disaster)` \u306e\u6c34\u6e96\u306f \"0\",\"1\",\"2\" \u306a\u306e\u3067\u3001\n\n")
w("- **`disaster1` = \u672a\u56de\u5fa9\uff08unrecovered\uff09n = %d**\n", sum(data$disaster1))
w("- **`disaster2` = \u56de\u5fa9\uff08recovered\uff09n = %d**\n\n", sum(data$disaster2))
w("\u3053\u306e\u9806\u5e8f\u306f\u76f4\u611f\u306b\u53cd\u3059\u308b\u304c\u7b2c1\u8ad6\u6587\u304b\u3089\u306e\u5f15\u304d\u7d99\u304e\u3067\u3042\u308b\u3002")
w("\u63a8\u5b9a\u5024\u3068\u3082\u6574\u5408\u3059\u308b\uff08\u672a\u56de\u5fa9 b = %.3f, p %s\uff1b\u56de\u5fa9 b = %.3f, p %s\uff09\u3002\n\n",
  rowsA$est[rowsA$rhs == "disaster1"], fmt_p(rowsA$pvalue[rowsA$rhs == "disaster1"]),
  rowsA$est[rowsA$rhs == "disaster2"], fmt_p(rowsA$pvalue[rowsA$rhs == "disaster2"]))
w("\u4ee5\u4e0b\u306e\u5168\u8868\u3067 `disaster_unrecovered` / `disaster_recovered` \u3068\u8868\u8a18\u3059\u308b\u3002\n\n")

w("## 1. \u6e2c\u5b9a\u30e2\u30c7\u30eb\u306e\u4ed5\u69d8\u5909\u66f4 \u2014 ResAll \u306e\u652a\u4e71\u9805\u3092 0 \u306b\u56fa\u5b9a\n\n")
w("\u4e3b\u30e2\u30c7\u30eb\u306f `ResAll ~~ 0*ResAll` \u3092\u8ab2\u3057\u305f\u7248\u306b\u5909\u66f4\u3057\u305f\u3002\n\n")
w("\u56fa\u5b9a\u524d\u306e\u81ea\u7531\u63a8\u5b9a\u5024\uff08\u6ce8\u8a18\u3068\u3057\u3066\u6b8b\u3059\uff09: **%.3f\uff08SE %.3f, p %.3f\uff09**\u3002",
  res_free$est[1], res_free$se[1], res_free$pvalue[1])
w("\u4e8c\u6b21\u56e0\u5b50 WB \u306e\u6307\u6a19\u304c NFE \u3068 ResAll \u306e 2 \u3064\u3057\u304b\u306a\u3044\u305f\u3081\u3053\u306e\u652a\u4e71\u9805\u306f\u5f31\u304f\u3057\u304b\u6c7a\u307e\u3089\u305a\u3001")
w("K6 \u3092\u542b\u307e\u306a\u3044\u5168\u4ed5\u69d8\u3067\u30bc\u30ed\u3068\u533a\u5225\u3067\u304d\u306a\u3044\u3002\n\n")
w("| \u30e2\u30c7\u30eb | ResAll \u652a\u4e71\u9805\uff08\u81ea\u7531\u63a8\u5b9a\uff09| p |\n|---|---|---|\n")
w("| Model A | -0.049 (SE 0.260) | 0.852 |\n")
w("| Model A2a | 0.178 (SE 0.202) | 0.376 |\n")
w("| Model A2b | 0.024 (SE 0.164) | 0.882 |\n")
w("| \u4ea4\u4e92\u4f5c\u7528\u30e2\u30c7\u30eb | -0.064 (SE 0.262) | 0.808 |\n")
w("| Model B \u7fa41 / \u7fa42 | -0.049 / 0.030 | 0.875 / 0.941 |\n")
w("| **K6 \u30e2\u30c7\u30eb** | **%.3f (SE %.3f)** | **%s** |\n\n",
  res_k6$est[1], res_k6$se[1], fmt_p(res_k6$pvalue[1]))
w("**K6 \u30e2\u30c7\u30eb\u3060\u3051\u306f\u4f8b\u5916\u3067\u3042\u308b\u3002**\u81ea\u7531\u63a8\u5b9a\u5024\u304c %.3f\uff08p < 0.001\uff09\u3068\u660e\u3089\u304b\u306b\u975e\u30bc\u30ed\u306a\u306e\u3067\u3001",
  res_k6$est[1])
w("\u3053\u3053\u3067 0 \u306b\u56fa\u5b9a\u3059\u308b\u306e\u306f\u771f\u306e\u8aa4\u5b9a\u5f0f\u5316\u306b\u306a\u308b\u3002\u3088\u3063\u3066 K6 \u30e2\u30c7\u30eb\u306f\u652a\u4e71\u9805\u3092\u81ea\u7531\u306e\u307e\u307e\u306b\u3059\u308b\u3002")
w("\u672c\u6587\u306b\u4e26\u8a18\u3059\u308b 2 \u5217\u306f\u3053\u306e 1 \u70b9\u3060\u3051\u4ed5\u69d8\u304c\u7570\u306a\u308b\u3002\n\n")
w("\u56fa\u5b9a\u306b\u3088\u308b\u5b9f\u5bb3\u306f\u306a\u3044\u3002\u30af\u30e9\u30b9\u4fc2\u6570\u306e\u6bd4\u8f03:\n\n")
w("| \u30af\u30e9\u30b9 | \u81ea\u7531\u63a8\u5b9a b [95%% CI] | beta | 0 \u56fa\u5b9a b [95%% CI] | beta |\n|---|---|---|---|---|\n")
for (k in 2:K) {
  a <- cls_row(rowsA, k); h <- cls_row(rowsA_free, k)
  w("| %d: %s | %s | %.3f | %s | %.3f |\n", k, CLASS_LAB[k],
    fmt_ci(h$est, h$ci.lower, h$ci.upper), h$std.all,
    fmt_ci(a$est, a$ci.lower, a$ci.upper), a$std.all)
}
w("\n\u9069\u5408: \u81ea\u7531 robust CFI %.3f / RMSEA %.3f / SRMR %.3f\u3001",
  fitA_free_meas[["cfi.robust"]], fitA_free_meas[["rmsea.robust"]], fitA_free_meas[["srmr"]])
w("0 \u56fa\u5b9a robust CFI %.3f / RMSEA %.3f / SRMR %.3f\u3002\n\n",
  fitA_meas[["cfi.robust"]], fitA_meas[["rmsea.robust"]], fitA_meas[["srmr"]])

w("## 2. Table 2 \u2014 \u30af\u30e9\u30b9\u5225\u30d7\u30ed\u30d5\u30a1\u30a4\u30eb\n\n")
w("n \u3068\u5272\u5f53 %% \u306f\u6700\u5927\u4e8b\u5f8c\u78ba\u7387\u306b\u3088\u308b\u5272\u5f53\u3001\u30e2\u30c7\u30eb\u30af\u30e9\u30b9\u6bd4\u306f\u30e2\u30c7\u30eb\u306e\u6df7\u5408\u6bd4\u7387\u3002\n\n")
w("| \u30af\u30e9\u30b9 | n | \u5272\u5f53 %% | \u30e2\u30c7\u30eb\u30af\u30e9\u30b9\u6bd4 %% | \u5e73\u5747\u652f\u6301\u9818\u57df\u6570 | \u30e2\u30c7\u30eb\u671f\u5f85\u652f\u6301\u6570 | %s |\n",
  paste(DOM_LAB, collapse = " | "))
w("|---|---|---|---|---|---|%s\n", paste(rep("---|", 9), collapse = ""))
for (k in seq_len(K)) {
  w("| %d: %s | %d | %.1f | %.1f | %.2f | %.2f | %s |\n", k, CLASS_LAB[k],
    t2$n[k], t2$pct[k], t2$class_share_model[k],
    t2$mean_n_domains[k], t2$model_expected_n_domains[k],
    paste(sprintf("%.2f", P[k, ]), collapse = " | "))
}
w("\n| \u30af\u30e9\u30b9 | \u91cd\u8981\u5ea6 Q11_3 | \u30b3\u30df\u30c3\u30c8 Q11_4 | \u610f\u5411 Q11_5 | \u4e00\u81f4\u5ea6 Q11_6 | \u5b9a\u7fa9\u3055\u308c\u308b n |\n")
w("|---|---|---|---|---|---|\n")
for (k in seq_len(K)) {
  w("| %d: %s | %.2f | %.2f | %.2f | %.2f | %d |\n", k, CLASS_LAB[k],
    t2$m_importance[k], t2$m_commitment[k], t2$m_intention[k],
    t2$m_consistency[k], t2$m_commitment_n[k])
}
w("\n**Table 2 \u811a\u6ce8**\u3002\u91cd\u8981\u5ea6\uff08Q11_3\u30015\u4ef6\u6cd5\uff09\u3068\u4e00\u81f4\u5ea6\uff08Q11_6\u30015\u5e2f\u306e\u767e\u5206\u7387\uff09\u306f")
w("\u5171\u901a\u5c3a\u5ea6\u4e0a\u306b\u306a\u304f\u3001\u5dee\u306f\u53d6\u3063\u3066\u3044\u306a\u3044\u3002Q11_3\u301c6 \u306f\u652f\u6301\u9818\u57df\u306b\u3057\u304b\u554f\u3046\u3066\u304a\u3089\u305a\u3001")
w("\u652f\u6301\u9818\u57df\u306a\u3057\u306e %d \u540d\u3067\u306f\u5b9a\u7fa9\u3055\u308c\u306a\u3044\u3002", sum(data$nsup == 0))
w("%d \u30af\u30e9\u30b9\u89e3\u306e entropy = **%.3f**\u3001\u30af\u30e9\u30b9\u5225\u306e\u5e73\u5747\u4e8b\u5f8c\u78ba\u7387\u306f %s \u3067\u3001",
  K, ent_K, paste(sprintf("\u30af\u30e9\u30b9%d %.3f", seq_len(K), avepp), collapse = "\u3001"))
w("\u6700\u3082\u5f31\u3044\u306e\u306f\u30af\u30e9\u30b9%d \u306e **%.2f**\u3002\n\n", worst_k, avepp[worst_k])
w("**Table 2 \u811a\u6ce8\uff08\u7d9a\u304d\u30fb\u69cb\u6210\u306e\u9006\u8ee2\u306e\u6271\u3044\uff09**\u3002\u30af\u30e9\u30b9 2 \u3068\u30af\u30e9\u30b9 3 \u306f\u5e73\u5747\u652f\u6301\u9818\u57df\u6570\u304c ")
w("%.2f \u3068 %.2f \u3067\u307b\u307c\u7b49\u3057\u304f\u3001\u6761\u4ef6\u4ed8\u304d\u652f\u6301\u78ba\u7387\u306f Health\uff08%.2f \u5bfe %.2f\uff09\u3068 ",
  mnd[2], mnd[3], P[2, 9], P[3, 9])
w("Romance\uff08%.2f \u5bfe %.2f\uff09\u3067\u9006\u8ee2\u3057\u3066\u3044\u308b\u3002\u305f\u3060\u3057\u3053\u306e 2 \u30af\u30e9\u30b9\u306e well-being \u5dee\u306f\u9811\u5065\u3067\u306f\u306a\u3044\u3002",
  P[2, 3], P[3, 3])
w("\u76f4\u63a5\u5bfe\u6bd4\u304c\u540d\u76ee\u4e0a\u6709\u610f\u306a\u306e\u306f Model A\uff08p %s\uff09\u3068\u5bb6\u65cf\u304c\u3093\u7d4c\u9a13\u3042\u308a\u7fa4\uff08p %s\uff09\u3060\u3051\u3067\u3001",
  fmt_p(con23$pvalue[con23$model == "Model A (n=2179)"]),
  fmt_p(con23$pvalue[con23$model == "Model B: family cancer (n=855)"]))
w("\u5e73\u5747\u30b3\u30df\u30c3\u30c8\u30e1\u30f3\u30c8\u8abf\u6574\u5f8c\u306f p = %.3f\u3001\u5bb6\u65cf\u304c\u3093\u7d4c\u9a13\u306a\u3057\u7fa4\u3067\u306f p = %.3f \u3068\u6d88\u3048\u308b\u3002",
  con23$pvalue[con23$model == "Model A2b (n=1964, + commitment)"],
  con23$pvalue[con23$model == "Model B: no family cancer (n=1324)"])
w("**\u69cb\u6210\u306e\u9006\u8ee2\u306f\u8a18\u8ff0\u3068\u3057\u3066\u63d0\u793a\u3057\u3001well-being \u5dee\u306e\u6839\u62e0\u3068\u3057\u3066\u306f\u7528\u3044\u306a\u3044\u3002**\n\n")
w("**Table 2 脚注（健康の支持確率）**。健康の条件付き支持確率は")
w("クラス 3 で %.3f、クラス 4 で %.3f とほぼ同等である（差 %.3f）。",
  P[3, 9], P[4, 9], abs(P[3, 9] - P[4, 9]))
w("したがって「健康を含む価値構造が well-being と結びつく」という説明は、")
w("nsup 調整後にクラス 3 だけが残りクラス 4 が消えることと整合しない。")
w("両クラスは健康をほぼ同じ確率で含んでいるからである。")
w("**クラス 3 を単独で中心主張に据えないこと。**\n\n")
w("**Table 2 \u811a\u6ce8\uff08\u5f37\u5ea6\u6307\u6a19\uff09**\u3002\u652f\u6301\u3057\u305f\u9818\u57df\u306b\u3064\u3044\u3066\u306e\u5f37\u5ea6\u306f 4 \u30af\u30e9\u30b9\u3067\u307b\u307c\u4e00\u5b9a\u3067\u3042\u308b\uff1a")
w("\u91cd\u8981\u5ea6 %.2f\u301c%.2f\u3001\u30b3\u30df\u30c3\u30c8\u30e1\u30f3\u30c8 %.2f\u301c%.2f\u3001\u610f\u5411 %.2f\u301c%.2f\u3001\u4e00\u81f4\u5ea6 %.2f\u301c%.2f\u3002",
  min(t2$m_importance), max(t2$m_importance),
  min(t2$m_commitment), max(t2$m_commitment),
  min(t2$m_intention), max(t2$m_intention),
  min(t2$m_consistency), max(t2$m_consistency))
w("\u3044\u305a\u308c\u3082 5 \u4ef6\u6cd5\u4e0a\u306e\u5e45\u304c %.2f \u4ee5\u4e0b\u3067\u3001\u30af\u30e9\u30b9\u9593\u306e\u5e73\u5747\u652f\u6301\u9818\u57df\u6570\u306e\u5e45\uff08%.2f\u301c%.2f\uff09\u3068\u306f\u5bfe\u7167\u7684\u3067\u3042\u308b\u3002\n\n",
  max(c(diff(range(t2$m_importance)), diff(range(t2$m_commitment)),
        diff(range(t2$m_intention)), diff(range(t2$m_consistency)))),
  min(t2$mean_n_domains), max(t2$mean_n_domains))
w("\u3057\u305f\u304c\u3063\u3066**\u30af\u30e9\u30b9 1 \u306f\u4fa1\u5024\u3078\u306e\u30b3\u30df\u30c3\u30c8\u30e1\u30f3\u30c8\u304c\u8584\u3044\u306e\u3067\u306f\u306a\u304f\u3001")
w("\u540c\u3058\u5f37\u5ea6\u3067\u8a8d\u3081\u308b\u9818\u57df\u304c\u5c11\u306a\u3044**\u3001\u3068\u8a18\u8ff0\u3067\u304d\u308b\u3002")
w("\u30af\u30e9\u30b9 1 \u3067\u30b3\u30df\u30c3\u30c8\u30e1\u30f3\u30c8\u304c\u5b9a\u7fa9\u3055\u308c\u308b %d \u540d\u306e\u5e73\u5747\u306f %.2f \u3067\u3001",
  t2$m_commitment_n[1], t2$m_commitment[1])
w("\u6700\u3082\u9ad8\u3044\u30af\u30e9\u30b9 %d \u306e %.2f \u3068\u306e\u5dee\u306f %.2f \u306b\u3059\u304e\u306a\u3044\u3002\n\n",
  which.max(t2$m_commitment), max(t2$m_commitment),
  max(t2$m_commitment) - t2$m_commitment[1])

w("## 3. Model A\uff08\u4e3b\u3001N = 2179\uff09\n\n")
w("\u9069\u5408: chi-square %.1f / %d df\u3001robust CFI %.3f\u3001robust RMSEA %.3f\u3001SRMR %.3f\n\n",
  fitA_meas[["chisq"]], fitA_meas[["df"]], fitA_meas[["cfi.robust"]],
  fitA_meas[["rmsea.robust"]], fitA_meas[["srmr"]])
w("| \u9805 | b [95%% CI] | beta | p |\n|---|---|---|---|\n")
for (i in seq_len(nrow(rowsA))) {
  r <- rowsA[i, ]
  w("| %s | %s | %.3f | %s |\n", r$term,
    fmt_ci(r$est, r$ci.lower, r$ci.upper), r$std.all, fmt_p(r$pvalue))
}
w("\n\u30af\u30e9\u30b9\u52b9\u679c\u306e\u540c\u6642\u691c\u5b9a: Wald chi-square %.2f, df %d, p %s\n\n",
  waldA$stat, waldA$df, fmt_p(waldA$p.value))

w("## 4. \u30af\u30e9\u30b9\u9593\u5bfe\u6bd4\u3068\u3001\u4e2d\u5fc3\u4e3b\u5f35\u306e\u9078\u629e\n\n")
w("\u3053\u306e\u7bc0\u304c\u300c\u6700\u5c0f\u5bfe\u306e\u4e8c\u5206\u300d\u304b\u300c\u30af\u30e9\u30b9\u9593\u306e\u52fe\u914d\u300d\u304b\u3092\u6c7a\u3081\u308b\u3002\n\n")
w("### 4a. H0: b2 = b3 = b4\uff08\u975e\u53c2\u7167\u30af\u30e9\u30b9\u306e\u52b9\u679c\u304c\u3059\u3079\u3066\u7b49\u3057\u3044\uff09\n\n")
w("\u68c4\u5374\u3055\u308c\u308c\u3070\u30af\u30e9\u30b9\u9593\u306b\u52fe\u914d\u304c\u3042\u308b\u3002\u68c4\u5374\u3055\u308c\u306a\u3051\u308c\u3070\u3001\u30c7\u30fc\u30bf\u304c\u652f\u6301\u3059\u308b\u306e\u306f")
w("\u300c\u30af\u30e9\u30b9 1\uff08Minimal\uff09\u3068\u305d\u308c\u4ee5\u5916\u304c\u9055\u3046\u300d\u3068\u3044\u3046\u4e8c\u5206\u307e\u3067\u3067\u3042\u308b\u3002\n\n")
w("| \u30e2\u30c7\u30eb | Wald chi-square | df | p |\n|---|---|---|---|\n")
w("| Model A (n=2179) | %.2f | %d | %s |\n", waldA_eq$stat, waldA_eq$df, fmt_p(waldA_eq$p.value))
w("| Model A2b (n=1964) | %.2f | %d | %s |\n", waldA2b_eq$stat, waldA2b_eq$df, fmt_p(waldA2b_eq$p.value))
w("| K6 モデル (n=2179) | %.2f | %d | %s |\n\n", waldK6_eq$stat, waldK6_eq$df, fmt_p(waldK6_eq$p.value))
w("\u53c2\u7167\uff08\u30af\u30e9\u30b9\u52b9\u679c\u304c\u3059\u3079\u3066 0\uff09: Model A Wald %.2f (df %d, p %s)\u3001Model A2b Wald %.2f (df %d, p %s)\u3002\n\n",
  waldA$stat, waldA$df, fmt_p(waldA$p.value), waldA2b$stat, waldA2b$df, fmt_p(waldA2b$p.value))
w("### 4b. \u5bfe\u30af\u30e9\u30b9\u5bfe\u6bd4\uff08\u5b9a\u7fa9\u30d1\u30e9\u30e1\u30fc\u30bf\uff09\n\n")
w("| \u5bfe\u6bd4 | \u30e2\u30c7\u30eb | b [95%% CI] | SE | Wald z | p |\n|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(con_all))) {
  r <- con_all[i, ]
  w("| %s | %s | %s | %.4f | %.2f | %s |\n", r$contrast, r$model,
    fmt_ci(r$b, r$ci.lower, r$ci.upper, 4), r$se, r$z, fmt_p(r$pvalue))
}
w("\n\u5e73\u5747\u652f\u6301\u9818\u57df\u6570: \u30af\u30e9\u30b92 %.2f\u3001\u30af\u30e9\u30b93 %.2f\u3001\u30af\u30e9\u30b94 %.2f\u3002\n\n",
  mnd[2], mnd[3], mnd[4])

w("### 4c. \u652f\u6301\u6570\u3092\u63c3\u3048\u305f\u5bfe\u6bd4\uff08nsup \u2208 {3, 4}\uff09\n\n")
w("\u69cb\u6210\u4eee\u8aac\u306e\u76f4\u63a5\u306e\u691c\u5b9a\u3002\u5e83\u3055\u3092\u56fa\u5b9a\u3057\u3066\u3082\u30af\u30e9\u30b9 2 \u3068 3 \u306e\u5dee\u304c\u6b8b\u308b\u304b\u3092\u898b\u308b\u3002\n\n")
w("\u90e8\u5206\u6a19\u672c\u306e\u30bb\u30eb\u5ea6\u6570:\n\n")
w("| \u30af\u30e9\u30b9 | nsup = 3 | nsup = 4 | \u8a08 |\n|---|---|---|---|\n")
for (k in seq_len(K)) w("| %d: %s | %d | %d | %d |\n", k, CLASS_LAB[k],
                        cell_n[k, "3"], cell_n[k, "4"], sum(cell_n[k, ]))
w("| **\u8a08** | %d | %d | **%d** |\n\n", sum(cell_n[, "3"]), sum(cell_n[, "4"]), sum(cell_n))
w("**\u30af\u30e9\u30b9 4 \u306f\u3053\u306e\u7bc4\u56f2\u306b %d \u540d\u3057\u304b\u304a\u3089\u305a\u3001\u30c0\u30df\u30fc\u3092\u8b58\u5225\u3067\u304d\u306a\u3044\u306e\u3067\u4e21\u30e2\u30c7\u30eb\u304b\u3089\u9664\u5916\u3057\u305f\u3002**",
  sum(cell_n[4, ]))
w("\u30af\u30e9\u30b9 1 \u306f %d \u540d\u3067\u3001\u4e3b\u5206\u6790\u304b\u3089\u306f\u5916\u3057\uff082 \u5bfe 3 \u306e\u5bfe\u6bd4\u306b\u5bc4\u4e0e\u305b\u305a\u63a8\u5b9a\u3092\u4e0d\u5b89\u5b9a\u306b\u3059\u308b\u3060\u3051\uff09\u3001",
  sum(cell_n[1, ]))
w("\u9664\u5916\u304c\u7d50\u679c\u3092\u4f5c\u3063\u3066\u3044\u306a\u3044\u3053\u3068\u3092\u793a\u3059\u305f\u3081\u306b\u526f\u6b21\u7684\u306b\u6295\u5165\u3057\u305f\u7248\u3082\u4f75\u8a18\u3059\u308b\u3002\n\n")
if (length(cv23$dropped))
  w("\u90e8\u5206\u6a19\u672c\u3067\u5b9a\u6570\u306b\u306a\u308a\u843d\u3068\u3057\u305f\u5171\u5909\u91cf: %s\u3002\n\n", paste(relab(cv23$dropped), collapse = ", "))
w("| \u30e2\u30c7\u30eb | n | b [95%% CI] | SE | beta | p |\n|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(match_tab))) {
  r <- match_tab[i, ]
  w("| %s | %d | %s | %.4f | %.3f | %s |\n", r$model, r$n,
    fmt_ci(r$b, r$ci_low, r$ci_high, 4), r$se, r$beta, fmt_p(r$p))
}
w("\n\u53c2\u8003: \u5168\u6a19\u672c\u3067\u306e 2 \u5bfe 3 \u306f Model A %s\uff08p %s\uff09\u3001A2b %s\uff08p %s\uff09\u3002\n\n",
  fmt_ci(con23$b[1], con23$ci.lower[1], con23$ci.upper[1], 4), fmt_p(con23$pvalue[1]),
  fmt_ci(con23$b[2], con23$ci.lower[2], con23$ci.upper[2], 4), fmt_p(con23$pvalue[2]))

w("### 4d. nsup \u3092\u8abf\u6574\u3057\u305f\u30af\u30e9\u30b9\u52b9\u679c\n\n")
w("(a) \u30af\u30e9\u30b9\u5185\u306e nsup \u5206\u5e03:\n\n")
w("| \u30af\u30e9\u30b9 | n | \u6700\u5c0f | Q1 | \u4e2d\u592e | Q3 | \u6700\u5927 | \u5e73\u5747 |\n|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(nsup_dist))) {
  r <- nsup_dist[i, ]
  w("| %d: %s | %d | %g | %g | %g | %g | %g | %.2f |\n", r$class, r$label, r$n,
    r$min, r$Q1, r$median, r$Q3, r$max, r$mean)
}
w("\n\u30af\u30e9\u30b9\u30c0\u30df\u30fc\u3067 nsup \u3092\u56de\u5e30\u3057\u305f R-squared = **%.4f**\u3002", r2_nsup_class)
w("\u56db\u5206\u4f4d\u7bc4\u56f2\u306f\u30af\u30e9\u30b9 3 \u3067 %g\u301c%g\u3001\u30af\u30e9\u30b9 4 \u3067 %g\u301c%g \u3068\u91cd\u306a\u308a\u304c\u307b\u3068\u3093\u3069\u306a\u3044\u3002\n\n",
  nsup_dist$Q1[3], nsup_dist$Q3[3], nsup_dist$Q1[4], nsup_dist$Q3[4])
w("(b) nsup \u3092\u5165\u308c\u305f\u30e2\u30c7\u30eb\u306e VIF\uff08GVIF^(1/(2*Df))\uff09:\n\n")
w("%s\n\n", paste(sprintf("%s %.2f", names(vif_ns), vif_ns), collapse = ", "))
w("(c) 情報量規準は §4h にまとめた。

")
w("\n\u30af\u30e9\u30b9\u30c0\u30df\u30fc 3 \u672c\u306e\u540c\u6642 Wald \u691c\u5b9a\uff08nsup \u8abf\u6574\u5f8c\uff09: Model A chi-square %.2f (df %d, p %s)\u3001",
  waldA_ns$stat, waldA_ns$df, fmt_p(waldA_ns$p.value))
w("Model A2b chi-square %.2f (df %d, p %s)\u3002\n\n",
  waldA2b_ns$stat, waldA2b_ns$df, fmt_p(waldA2b_ns$p.value))
w("| \u30af\u30e9\u30b9 | nsup \u8abf\u6574\u524d b [95%% CI] | nsup \u8abf\u6574\u5f8c b [95%% CI] | p |\n|---|---|---|---|\n")
for (k in 2:K) {
  a <- cls_row(rowsA, k); b <- cls_row(rowsA_ns, k)
  w("| %d: %s | %s | %s | %s |\n", k, CLASS_LAB[k],
    fmt_ci(a$est, a$ci.lower, a$ci.upper), fmt_ci(b$est, b$ci.lower, b$ci.upper),
    fmt_p(b$pvalue))
}
ns_row <- rowsA_ns[rowsA_ns$rhs == "nsup_z", ]
w("\nnsup \u81ea\u4f53\uff081 SD \u3042\u305f\u308a\uff09: b = %s, p = %s\n\n",
  fmt_ci(ns_row$est, ns_row$ci.lower, ns_row$ci.upper), fmt_p(ns_row$pvalue))

w("### 4e. \u6c7a\u5b9a\u7684\u691c\u5b9a \u2014 nsup \u8abf\u6574\u5f8c\u306e H0: b2 = b3 = b4\n\n")
w("\u7b49\u3057\u3051\u308c\u3070\u3001\u69cb\u6210\u306f\u30af\u30e9\u30b9\u6570\u3092\u8d85\u3048\u308b\u60c5\u5831\u3092\u8db3\u3057\u3066\u304a\u3089\u305a\u3001")
w("\u300c\u6700\u5c0f \u5bfe \u69cb\u9020\u5316\u300d\u306e\u4e8c\u5024\u6307\u6a19\u3067\u8db3\u308a\u308b\u3053\u3068\u306b\u306a\u308b\u3002\n\n")
w("| \u30e2\u30c7\u30eb | Wald chi-square | df | p |\n|---|---|---|---|\n")
w("| Model A + nsup (n=2179) | %.2f | %d | %s |\n",
  waldA_ns_eq$stat, waldA_ns_eq$df, fmt_p(waldA_ns_eq$p.value))
w("| Model A2b + nsup (n=1964) | %.2f | %d | %s |\n\n",
  waldA2b_ns_eq$stat, waldA2b_ns_eq$df, fmt_p(waldA2b_ns_eq$p.value))
w("nsup \u8abf\u6574\u524d\u306e\u540c\u3058\u691c\u5b9a\u306f Model A p = %.3f\u3001Model A2b p = %.3f \u3067\u3042\u3063\u305f\u3002\n\n",
  waldA_eq$p.value, waldA2b_eq$p.value)

w("### 4f. \u30af\u30e9\u30b9\u3068 nsup \u306e\u91cd\u306a\u308a\uff08limitation \u306e\u6839\u62e0\uff09\n\n")
w("| \u30af\u30e9\u30b9 | %s | \u8a08 |\n|---|%s---|\n", paste(0:9, collapse = " | "),
  paste(rep("---|", 10), collapse = ""))
for (k in seq_len(K)) w("| %d: %s | %s | %d |\n", k, CLASS_LAB[k],
                        paste(ov[k, ], collapse = " | "), sum(ov[k, ]))
w("\n\u30af\u30e9\u30b9\u30c0\u30df\u30fc\u3067 nsup \u3092\u56de\u5e30\u3057\u305f **R-squared = %.4f**\u3002", r2_nsup_class)
w("\u30af\u30e9\u30b9 3 \u3068\u30af\u30e9\u30b9 4 \u306e nsup \u5206\u5e03\u306e\u91cd\u306a\u308a\u4fc2\u6570\u306f **%.3f**\uff08= %.1f%%\uff09\u306b\u3059\u304e\u306a\u3044\u3002\n\n",
  overlap_34, 100 * overlap_34)
w("nsup >= 6 \u306e %d \u540d\u306e\u3046\u3061 **%.1f%%\uff08%d \u540d\uff09\u304c\u30af\u30e9\u30b9 4** \u3067\u3001\u305d\u308c\u4ee5\u5916\u306e\u30af\u30e9\u30b9\u306f\u5408\u8a08 %d \u540d\u3057\u304b\u3044\u306a\u3044\u3002",
  sum(ov[, hi]), 100 * share_c4_hi, sum(ov[4, hi]), n_hi_other)
w("\u3057\u305f\u304c\u3063\u3066 nsup \u3092\u8abf\u6574\u3057\u305f\u3068\u304d\u306e\u30af\u30e9\u30b9 4 \u306e\u4fc2\u6570\u306f\u3001\u540c\u3058\u652f\u6301\u6570\u3092\u6301\u3064\u4ed6\u30af\u30e9\u30b9\u3068\u306e")
w("\u5bfe\u6bd4\u3067\u306f\u306a\u304f\u3001\u307b\u307c**\u5916\u633f**\u3067\u3042\u308b\u3002\u00a74d \u3067\u30af\u30e9\u30b9 4 \u306e\u4fc2\u6570\u304c %s \u3068\u6d88\u3048\u308b\u306e\u306f\u3001",
  fmt_ci(cls_row(rowsA_ns, 4)$est, cls_row(rowsA_ns, 4)$ci.lower,
         cls_row(rowsA_ns, 4)$ci.upper))
w("\u52b9\u679c\u304c\u7121\u3044\u3053\u3068\u306e\u8a3c\u62e0\u3067\u306f\u306a\u304f\u3001\u6bd4\u8f03\u5bfe\u8c61\u304c\u5c45\u306a\u3044\u3053\u3068\u306e\u5e30\u7d50\u3068\u3057\u3066\u8aad\u3080\u3079\u304d\u3067\u3042\u308b\u3002\n\n")

w("### 4g. \u4e3b\u66dd\u9732\u3092\u652f\u6301\u6570\u306b\u5207\u308a\u66ff\u3048\u305f\u7248\n\n")
w("\u30af\u30e9\u30b9\u3092\u4f7f\u308f\u305a\u3001nsup \u3092 7 \u6c34\u6e96\u306e\u30ab\u30c6\u30b4\u30ea\uff080/1/2/3/4/5/6\u4ee5\u4e0a\uff09\u3068\u3057\u3066")
w("Model A \u306e\u5171\u5909\u91cf\u30bb\u30c3\u30c8\u3067\u56de\u5e30\u3057\u305f\u3002\u56f3\u306f `outputs/fig2_nsup_categorical.png`\u3002\n\n")
w("\u6c34\u6e96\u3054\u3068\u306e n: %s\n\n",
  paste(sprintf("%s: %d (なし %d / あり %d)", names(nscat_n),
                as.integer(nscat_n), nscat_n_grp[, 1], nscat_n_grp[, 2]),
        collapse = "、"))
for (mm in unique(ncat$model)) {
  z <- ncat[ncat$model == mm, ]
  w("**%s**\n\n", mm)
  w("| \u652f\u6301\u6570 | n | b [95%% CI] | p |\n|---|---|---|---|\n")
  for (i in seq_len(nrow(z))) {
    w("| %s | %d | %s | %s |\n", c("1","2","3","4","5","6+")[z$level[i]],
      z$n_level[i], fmt_ci(z$est[i], z$ci.lower[i], z$ci.upper[i]),
      fmt_p(z$pvalue[i]))
  }
  w("\n")
}
w("\u6c34\u6e96\u304c\u3059\u3079\u3066 0: Wald chi-square %.2f (df %d, p %s)\u3002",
  wald_ncat$stat, wald_ncat$df, fmt_p(wald_ncat$p.value))
w("\u4e0a\u4f4d 3 \u6c34\u6e96\u304c\u7b49\u3057\u3044\uff084 = 5 = 6\u4ee5\u4e0a\u3001\u98fd\u548c\uff09: Wald chi-square %.2f (df %d, p %s)\u3002\n\n",
  wald_ncat_flat$stat, wald_ncat_flat$df, fmt_p(wald_ncat_flat$p.value))

w("### 4h. \u60c5\u5831\u91cf\u898f\u6e96\n\n")
w("| \u30e2\u30c7\u30eb | npar | AIC | dAIC | BIC | dBIC |\n|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(ic_tab))) {
  r <- ic_tab[i, ]
  w("| %s | %d | %.1f | %.2f | %.1f | %.2f |\n", r$model, r$npar, r$AIC, r$dAIC,
    r$BIC, r$dBIC)
}
w("\n")

w("### 4i. 構成仮説の扱い — 判定不能

")
w("**構成仮説は「判定不能」として扱う。支持もされていないが、否定もされていない。**

")
w("- 支持数を揃えた対比（4c）は %s、p = %.3f、n = %d。",
  fmt_ci(match_tab$b[1], match_tab$ci_low[1], match_tab$ci_high[1], 4),
  match_tab$p[1], match_tab$n[1])
w("**これは検出力不足であって零の証拠ではない。**点推定は全標本の %.4f からほとんど動いておらず",
  con23$b[1])
w("（%.4f、変化 %.1f%%）、区間は全標本推定値を余裕をもって含んでいる。",
  match_tab$b[1], 100 * (match_tab$b[1] - con23$b[1]) / con23$b[1])
w("SE が %.4f から %.4f に広がっただけである。
",
  con23$se[1], match_tab$se[1])
w("- nsup 調整後に b2 = b3 = b4 が棄却されない（4e、p = %.3f / %.3f）ことは、",
  waldA_ns_eq$p.value, waldA2b_ns_eq$p.value)
w("構成がクラス数を超える情報を足していないことを示す。ただしこれも「差が無い」の証明ではない。
")
w("- クラス 4 の係数は外挿（4f）なので、そこから構成の可否を読んではならない。
")
w("- 健康の支持確率はクラス 3 と 4 でほぼ同等（%.3f 対 %.3f）なので、",
  P[3, 9], P[4, 9])
w("「健康を含む価値構造」という説明も成り立たない。**クラス 3 を単独で中心主張に据えない。**

")
w("残るのは支持の**広さ**であり、その用量反応は 4 領域あたりで飽和する")
w("（水準 4 = 5 = 6以上 の等値検定 p = %.3f）。

", wald_ncat_flat$p.value)

w("## 5. 1 \u9818\u57df\u3042\u305f\u308a\u306e\u50be\u304d\n\n")
w("\u96a3\u63a5\u30af\u30e9\u30b9\u9593\u306e b \u306e\u5dee\u3092\u3001\u5e73\u5747\u652f\u6301\u9818\u57df\u6570\u306e\u5dee\u3067\u5272\u3063\u305f\u3082\u306e\u3002")
w("well-being \u304c\u652f\u6301\u306e\u5e83\u3055\u3060\u3051\u3067\u6c7a\u307e\u308b\u306a\u3089\u3001\u3053\u306e\u50be\u304d\u306f 3 \u3064\u306e\u6bb5\u3067\u307b\u307c\u4e00\u5b9a\u306b\u306a\u308b\u306f\u305a\u3067\u3042\u308b\u3002\n\n")
w("| \u6bb5 | \u5e73\u5747\u652f\u6301\u6570 | \u9818\u57df\u5dee | b \u306e\u5dee [95%% CI] | p | **1 \u9818\u57df\u3042\u305f\u308a\u50be\u304d** [95%% CI] |\n")
w("|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(slope))) {
  s <- slope[i, ]
  w("| %s | %.2f \u2192 %.2f | %.2f | %s | %s | **%.4f** [%.4f, %.4f] |\n",
    s$step, s$mean_domains_from, s$mean_domains_to, s$domain_gap,
    fmt_ci(s$b_difference, s$ci_low, s$ci_high, 4), fmt_p(s$p),
    s$slope_per_domain, s$slope_ci_low, s$slope_ci_high)
}
w("\n\u50be\u304d\u306f\u4e00\u5b9a\u3067\u306f\u306a\u304f\u3001\u6700\u5927\u3068\u6700\u5c0f\u3067 **%.1f \u500d** \u9055\u3046\uff08%.4f \u5bfe %.4f\uff09\u3002",
  max(slope$slope_per_domain) / min(slope$slope_per_domain),
  max(slope$slope_per_domain), min(slope$slope_per_domain))
w("\u652f\u6301\u9818\u57df\u6570\u304c\u6700\u3082\u5c0f\u3055\u304f\u3057\u304b\u5897\u3048\u306a\u3044\u30af\u30e9\u30b92\u21923 \u306e\u6bb5\u3067\u50be\u304d\u304c\u6700\u5927\u306b\u306a\u308a\u3001")
w("\u652f\u6301\u6570\u304c\u6700\u3082\u5927\u304d\u304f\u5897\u3048\u308b\u30af\u30e9\u30b93\u21924 \u306e\u6bb5\u3067\u6700\u5c0f\u306b\u306a\u308b\u3002")
w("\u652f\u6301\u6570\u3060\u3051\u3067\u306f\u8aac\u660e\u3067\u304d\u306a\u3044\u3002\n\n")

w("## 6. Model A2\uff08n = 1964\uff09\u2014 \u5e73\u5747\u30b3\u30df\u30c3\u30c8\u30e1\u30f3\u30c8\u8abf\u6574\n\n")
w("| \u30af\u30e9\u30b9 | A\uff08n=2179\uff09b [95%% CI] | A2a\uff08n=1964\uff09b [95%% CI] | A2b\uff08+\u30b3\u30df\u30c3\u30c8\uff09b [95%% CI] | A2b beta | A2b p |\n")
w("|---|---|---|---|---|---|\n")
for (k in 2:K) {
  a <- cls_row(rowsA, k); b <- cls_row(rowsA2a, k); c2 <- cls_row(rowsA2b, k)
  w("| %d: %s | %s | %s | %s | %.3f | %s |\n", k, CLASS_LAB[k],
    fmt_ci(a$est, a$ci.lower, a$ci.upper), fmt_ci(b$est, b$ci.lower, b$ci.upper),
    fmt_ci(c2$est, c2$ci.lower, c2$ci.upper), c2$std.all, fmt_p(c2$pvalue))
}
mc <- rowsA2b[rowsA2b$rhs == "mcom_z", ]
w("\n\u5e73\u5747\u30b3\u30df\u30c3\u30c8\u30e1\u30f3\u30c8\uff081 SD\uff09: b = %s, beta = %.3f, p = %s\n",
  fmt_ci(mc$est, mc$ci.lower, mc$ci.upper), mc$std.all, fmt_p(mc$pvalue))
w("\u30af\u30e9\u30b9\u52b9\u679c\u306e\u540c\u6642\u691c\u5b9a: A2a Wald %.2f (df %d, p %s)\u3001A2b Wald %.2f (df %d, p %s)\n\n",
  waldA2a$stat, waldA2a$df, fmt_p(waldA2a$p.value),
  waldA2b$stat, waldA2b$df, fmt_p(waldA2b$p.value))

w("## 7. Model B\uff08\u591a\u7fa4\u3001\u5bb6\u65cf\u306e\u304c\u3093\u7d4c\u9a13\u5225\uff09\n\n")
w("| \u30af\u30e9\u30b9 | \u304c\u3093\u7d4c\u9a13\u306a\u3057 (n=1324) b [95%% CI] | beta | p | \u304c\u3093\u7d4c\u9a13\u3042\u308a (n=855) b [95%% CI] | beta | p |\n")
w("|---|---|---|---|---|---|---|\n")
for (k in 2:K) {
  b0 <- cls_row(rowsB0, k); b1 <- cls_row(rowsB1, k)
  w("| %d: %s | %s | %.3f | %s | %s | %.3f | %s |\n", k, CLASS_LAB[k],
    fmt_ci(b0$est, b0$ci.lower, b0$ci.upper), b0$std.all, fmt_p(b0$pvalue),
    fmt_ci(b1$est, b1$ci.lower, b1$ci.upper), b1$std.all, fmt_p(b1$pvalue))
}
w("\n2\u7fa4\u9593\u3067\u30af\u30e9\u30b9\u4fc2\u6570\u304c\u7b49\u3057\u3044: Wald chi-square %.2f, df %d, p %s\n\n",
  waldB_eq$stat, waldB_eq$df, fmt_p(waldB_eq$p.value))

w("### 7b. \u4e8c\u6b21\u56e0\u5b50 WB \u306e\u5206\u6563\uff08\u7fa4\u5225\uff09\n\n")
w("\u975e\u6a19\u6e96\u5316\u4fc2\u6570\u304c\u304c\u3093\u7d4c\u9a13\u7fa4\u3067\u5927\u304d\u3044\u306e\u304c\u3001\u5358\u306b\u305d\u306e\u7fa4\u3067 WB \u306e\u5206\u6563\u304c\u5927\u304d\u3044\u305b\u3044\u306a\u306e\u304b\u3092\u5224\u5225\u3059\u308b\u3002\n\n")
w("| \u7fa4 | WB \u5168\u5206\u6563\uff08\u30e2\u30c7\u30eb\u542b\u610f\uff09| WB \u6b8b\u5dee\u5206\u6563 | SD(WB) |\n|---|---|---|---|\n")
gl <- c("\u304c\u3093\u7d4c\u9a13\u306a\u3057 (n=1324)", "\u304c\u3093\u7d4c\u9a13\u3042\u308a (n=855)")
for (g in 1:2) w("| %s | %.5f | %.5f | %.4f |\n", gl[g], wb_var_total[g],
                 wb_var_res[g], wb_sd_total[g])
w("\nSD(WB) \u306e\u6bd4\uff08\u3042\u308a / \u306a\u3057\uff09= %.3f\u3002", wb_sd_total[2] / wb_sd_total[1])
w("\u30af\u30e9\u30b9\u5225\u306b\u898b\u308b\u3068:\n\n")
w("| \u30af\u30e9\u30b9 | b \u6bd4\uff08\u3042\u308a/\u306a\u3057\uff09| beta \u6bd4\uff08\u3042\u308a/\u306a\u3057\uff09|\n|---|---|---|\n")
for (k in 2:K) {
  b0 <- cls_row(rowsB0, k); b1 <- cls_row(rowsB1, k)
  w("| %d: %s | %.2f | %.2f |\n", k, CLASS_LAB[k],
    b1$est / b0$est, b1$std.all / b0$std.all)
}
w("\n\u5206\u6563\u306e\u9055\u3044\u3060\u3051\u304c\u539f\u56e0\u306a\u3089 b \u6bd4\u306f SD \u6bd4\u3067\u8aac\u660e\u3055\u308c\u3001beta \u6bd4\u306f 1 \u306b\u8fd1\u3065\u304f\u306f\u305a\u3067\u3042\u308b\u3002\n\n")

w("### 7c. \u7fa4\u9593\u6bd4\u306f\u4fa1\u5024\u30af\u30e9\u30b9\u306b\u7279\u7570\u7684\u304b\n\n")
w("\u304c\u3093\u7d4c\u9a13\u7fa4\u3067\u4fc2\u6570\u304c\u5927\u304d\u3044\u3053\u3068\u304c\u4fa1\u5024\u30af\u30e9\u30b9\u56fa\u6709\u306e\u73fe\u8c61\u306a\u306e\u304b\u3001\u305d\u306e\u7fa4\u3067\u306f")
w("\u3069\u306e\u4e88\u6e2c\u5909\u6570\u3082\u5927\u304d\u304f\u51fa\u308b\u306e\u304b\u3092\u898b\u308b\u3002\u5f8c\u8005\u306a\u3089\u3001\u3053\u306e\u7fa4\u5dee\u306f\u4fa1\u5024\u30af\u30e9\u30b9\u306b\u3064\u3044\u3066\u306e")
w("\u77e5\u898b\u3067\u306f\u306a\u304f\u7fa4\u306b\u3064\u3044\u3066\u306e\u77e5\u898b\u3067\u3042\u308b\u3002\n\n")
w("| \u9805 | \u306a\u3057 b | \u3042\u308a b | **b \u6bd4** | \u306a\u3057 beta | \u3042\u308a beta | beta \u6bd4 |\n")
w("|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(grp_cmp))) {
  r <- grp_cmp[i, ]
  w("| %s | %.4f | %.4f | **%.2f** | %.3f | %.3f | %.2f |\n", r$term,
    r$b_nocancer, r$b_cancer, r$b_ratio, r$beta_nocancer, r$beta_cancer, r$beta_ratio)
}
cls_ratios <- grp_cmp$b_ratio[grp_cmp$term %in% relab(CLSD)]
oth_ratios <- grp_cmp$b_ratio[!grp_cmp$term %in% relab(CLSD)]
w("\n\u30af\u30e9\u30b9\u30c0\u30df\u30fc\u306e b \u6bd4\u306f %s\u3001\u305d\u308c\u4ee5\u5916\u306e\u4e88\u6e2c\u5909\u6570\u306f %s\u3002\n\n",
  paste(sprintf("%.2f", cls_ratios), collapse = " / "),
  paste(sprintf("%.2f", oth_ratios), collapse = " / "))

w("## 8. K6 \u306e\u6709\u7121\uff08\u672c\u6587\u306b\u4e26\u8a18\uff09\n\n")
w("\u4e3b\u30e2\u30c7\u30eb\u306f K6 \u3092\u5165\u308c\u306a\u3044\u3002K6 \u5165\u308a\u3092\u4e26\u3079\u308b\u3002")
w("\u4ed5\u69d8\u306e\u9055\u3044\u306f K6 \u306e\u6709\u7121\u3068\u3001\u4e0a\u306e\u00a71 \u3067\u8ff0\u3079\u305f ResAll \u652a\u4e71\u9805\u306e\u6271\u3044\u306e 2 \u70b9\u3060\u3051\u3002\n\n")
w("K6 \u30e2\u30c7\u30eb\u9069\u5408: chi-square %.1f / %d df\u3001robust CFI %.3f\u3001RMSEA %.3f\u3001SRMR %.3f\n\n",
  fitK6_meas[["chisq"]], fitK6_meas[["df"]], fitK6_meas[["cfi.robust"]],
  fitK6_meas[["rmsea.robust"]], fitK6_meas[["srmr"]])
w("| \u9805 | K6 \u306a\u3057 b [95%% CI] | beta | p | K6 \u3042\u308a b [95%% CI] | beta | p | beta \u306e\u5909\u5316 |\n")
w("|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(k6_tab))) {
  r <- k6_tab[i, ]
  if (is.na(r$est_noK6)) {
    w("| %s | \u2014 | \u2014 | \u2014 | %s | %.3f | %s | \u2014 |\n", r$term,
      fmt_ci(r$est_K6, r$ci.lower_K6, r$ci.upper_K6), r$std.all_K6,
      fmt_p(r$pvalue_K6))
  } else {
    w("| %s | %s | %.3f | %s | %s | %.3f | %s | %+.3f |\n", r$term,
      fmt_ci(r$est_noK6, r$ci.lower_noK6, r$ci.upper_noK6), r$std.all_noK6,
      fmt_p(r$pvalue_noK6),
      fmt_ci(r$est_K6, r$ci.lower_K6, r$ci.upper_K6), r$std.all_K6,
      fmt_p(r$pvalue_K6), r$d_beta)
  }
}
w("\n\u30af\u30e9\u30b9\u52b9\u679c\u306e\u540c\u6642\u691c\u5b9a: K6 \u306a\u3057 Wald %.2f (df %d, p %s)\u3001K6 \u3042\u308a Wald %.2f (df %d, p %s)\n\n",
  waldA$stat, waldA$df, fmt_p(waldA$p.value),
  waldK6$stat, waldK6$df, fmt_p(waldK6$p.value))
w("\u30af\u30e9\u30b9 2 \u5bfe 3 \u306e\u5bfe\u6bd4\u306f K6 \u30e2\u30c7\u30eb\u3067 %s\uff08p %s\uff09\u3002\n\n",
  fmt_ci(con23$b[con23$model == "K6 model (n=2179)"],
         con23$ci.lower[con23$model == "K6 model (n=2179)"],
         con23$ci.upper[con23$model == "K6 model (n=2179)"], 4),
  fmt_p(con23$pvalue[con23$model == "K6 model (n=2179)"]))

w("## 9. \u4e8b\u524d\u6307\u5b9a\u306e\u4ea4\u4e92\u4f5c\u7528\uff08\u4fa1\u5024\u30af\u30e9\u30b9 \u00d7 \u5fa9\u8208\u72b6\u6cc1\uff09\n\n")
w("\u5168\u4f53\u306e Wald \u691c\u5b9a: chi-square %.2f, df %d, p %s\n\n",
  waldI$stat, waldI$df, fmt_p(waldI$p.value))
w("| \u30af\u30e9\u30b9 | \u88ab\u5bb3\u306a\u3057 | \u672a\u56de\u5fa9 | \u56de\u5fa9 |\n|---|---|---|---|\n")
for (k in seq_len(K)) w("| %d: %s | %d | %d | %d |\n", k, CLASS_LAB[k],
                        cells[k, "unexposed"], cells[k, "unrecovered"], cells[k, "recovered"])
w("\n\u6700\u5c0f\u30bb\u30eb n = %d\u3002\u63a2\u7d22\u7684\u3068\u660e\u793a\u3059\u308b\u3053\u3068\u3002\n\n", min(cells))

w("## 10. \u591a\u91cd\u5171\u7dda\u6027\uff08VIF\uff09\n\n")
w("Model A: %s\n\n", paste(sprintf("%s %.2f", names(vif_main), vif_main), collapse = ", "))
w("\u4ea4\u4e92\u4f5c\u7528\u30e2\u30c7\u30eb: \u7bc4\u56f2 %.2f \u2013 %.2f\n\n", min(vif_ix), max(vif_ix))

w("## 11. supplementary: \u4e00\u6b21\u56e0\u5b50\u3078\u306e\u30af\u30e9\u30b9\u52b9\u679c\n\n")
w("\u4e8c\u6b21\u56e0\u5b50\u3092\u7f6e\u304b\u305a\u3001NFE\uff08negative-free emotion\uff09\u3068 ResAll\uff08resource-oriented well-being\uff09\u3092")
w("\u305d\u308c\u305e\u308c\u540c\u3058\u4e88\u6e2c\u5909\u6570\u306b\u56de\u5e30\u3057\u3001\u652a\u4e71\u9805\u3092\u76f8\u95a2\u3055\u305b\u305f\u30e2\u30c7\u30eb\u3002\n\n")
w("\u9069\u5408: chi-square %.1f / %d df\u3001robust CFI %.3f\u3001RMSEA %.3f\u3001SRMR %.3f\n\n",
  fitFO_meas[["chisq"]], fitFO_meas[["df"]], fitFO_meas[["cfi.robust"]],
  fitFO_meas[["rmsea.robust"]], fitFO_meas[["srmr"]])
w("| \u30af\u30e9\u30b9 | NFE b [95%% CI] | beta | p | ResAll b [95%% CI] | beta | p |\n")
w("|---|---|---|---|---|---|---|\n")
for (k in 2:K) {
  n1 <- cls_row(rowsFO_n, k); r1 <- cls_row(rowsFO_r, k)
  w("| %d: %s | %s | %.3f | %s | %s | %.3f | %s |\n", k, CLASS_LAB[k],
    fmt_ci(n1$est, n1$ci.lower, n1$ci.upper), n1$std.all, fmt_p(n1$pvalue),
    fmt_ci(r1$est, r1$ci.lower, r1$ci.upper), r1$std.all, fmt_p(r1$pvalue))
}
w("\n\u30af\u30e9\u30b9\u52b9\u679c\u306e\u540c\u6642\u691c\u5b9a: NFE \u3078 Wald %.2f (df %d, p %s)\u3001ResAll \u3078 Wald %.2f (df %d, p %s)\n\n",
  waldFO_n$stat, waldFO_n$df, fmt_p(waldFO_n$p.value),
  waldFO_r$stat, waldFO_r$df, fmt_p(waldFO_r$p.value))
w("**\u3053\u306e\u30e2\u30c7\u30eb\u306e\u4ed5\u69d8\uff08\u91cd\u8981\uff09**\u3002\u4e8c\u968e\u56e0\u5b50\u306f\u7f6e\u3044\u3066\u3044\u306a\u3044\u3002\u6f5c\u5728\u5909\u6570\u306f NFE \u3068 ResAll \u306e")
w("2 \u3064\u3060\u3051\u3067\u3001`WB =~ NFE + ResAll` \u306f\u5b58\u5728\u305b\u305a\u3001ResAll \u306e\u652a\u4e71\u9805\u3082 0 \u306b\u56fa\u5b9a\u3057\u3066\u3044\u306a\u3044")
w("\uff08\u81ea\u7531\u63a8\u5b9a %.3f\u3001SE %.3f\u3001p %s\uff09\u3002\u4e21\u652a\u4e71\u9805\u306f\u76f8\u95a2\u3092\u8a31\u3057\u3066\u3044\u308b\uff08%.3f\u3001SE %.3f\uff09\u3002\n\n",
  fo_dist$est[1], fo_dist$se[1], fmt_p(fo_dist$pvalue[1]), fo_cov$est[1], fo_cov$se[1])
w("\u3053\u308c\u304c\u91cd\u8981\u306a\u306e\u306f\u3001\u3082\u3057\u4e8c\u968e\u56e0\u5b50\u30e2\u30c7\u30eb\u4e0a\u3067\u3001\u3057\u304b\u3082 ResAll \u306e\u652a\u4e71\u9805\u3092 0 \u306b\u56fa\u5b9a\u3057\u305f")
w("\u72b6\u614b\u3067\u63a8\u5b9a\u3057\u3066\u3044\u305f\u306a\u3089\u3001NFE \u3078\u306e\u7d4c\u8def\u306f\u8a2d\u5b9a\u4e0a\u4e8c\u968e\u56e0\u5b50\u7d4c\u7531\u306b\u5f37\u5236\u3055\u308c\u3001NFE \u304c")
w("null \u3068\u3044\u3046\u7d50\u679c\u306f\u4eba\u70ba\u7684\u306b\u306a\u308b\u304b\u3089\u3067\u3042\u308b\u3002\u305d\u3046\u3067\u306f\u306a\u3044\u3053\u3068\u3092\u30b9\u30af\u30ea\u30d7\u30c8\u5185\u3067")
w("`stopifnot` \u306b\u3088\u308a\u8868\u660e\u3057\u3066\u3044\u308b\uff08\u6f5c\u5728\u5909\u6570\u304c NFE \u3068 ResAll \u306e\u307f\u3001lhs == \"WB\" \u306e")
w("\u30d1\u30e9\u30e1\u30fc\u30bf\u884c\u304c 0 \u672c\u3001ResAll \u652a\u4e71\u9805\u306e SE > 0 \u304b\u3064\u63a8\u5b9a\u5024 > 0\uff09\u3002\n\n")
w("\u5168\u9805\u306f outputs/table6_firstorder.csv\u3002\n\n")

w("### 11b. 広さを曝露とした一次因子分解（完全結果）

")
w("適合: chi-square %.1f / %d df、robust CFI %.3f、RMSEA %.3f、SRMR %.3f

",
  fitFOns_meas[["chisq"]], fitFOns_meas[["df"]], fitFOns_meas[["cfi.robust"]],
  fitFOns_meas[["rmsea.robust"]], fitFOns_meas[["srmr"]])
w("**尺度の注意**。NFE のマーカーは二値項目、ResAll のマーカーは 0-10 項目なので、")
w("両因子の非標準化係数は相互に比較できず、Table 3 の二次因子係数とも比較できない。")
w("比較は標準化係数で行うこと。

")
w("| 支持数 | NFE b [95%% CI] | NFE beta | NFE p | ResAll b [95%% CI] | ResAll beta | ResAll p |
")
w("|---|---|---|---|---|---|---|
")
for (lv in 1:6) {
  a <- rowsFOns_n[rowsFOns_n$rhs == paste0("nsc", lv), ]
  b <- rowsFOns_r[rowsFOns_r$rhs == paste0("nsc", lv), ]
  w("| %s | %s | %.3f | %s | %s | %.3f | %s |
",
    c("1","2","3","4","5","6+")[lv],
    fmt_ci(a$est, a$ci.lower, a$ci.upper), a$std.all, fmt_p(a$pvalue),
    fmt_ci(b$est, b$ci.lower, b$ci.upper), b$std.all, fmt_p(b$pvalue))
}
w("
同時検定: NFE へ Wald %.2f (df %d, p %s)、ResAll へ Wald %.2f (df %d, p %s)

",
  waldFOns_n$stat, waldFOns_n$df, fmt_p(waldFOns_n$p.value),
  waldFOns_r$stat, waldFOns_r$df, fmt_p(waldFOns_r$p.value))
w("標準化解による NFE の係数と 95%% CI:

")
w("| 支持数 | NFE beta [95%% CI] | ResAll beta [95%% CI] |
|---|---|---|
")
for (lv in 1:6) {
  a <- stdFOns_n[stdFOns_n$rhs == paste0("nsc", lv), ]
  b <- stdFOns_r[stdFOns_r$rhs == paste0("nsc", lv), ]
  w("| %s | %s | %s |
", c("1","2","3","4","5","6+")[lv],
    fmt_ci(a$est.std, a$ci.lower, a$ci.upper),
    fmt_ci(b$est.std, b$ci.lower, b$ci.upper))
}
w("
**判定**: NFE の標準化区間の上限は最大で %.3f、", nfe_upper)
w("ResAll の標準化係数は %.3f 〜 %.3f。",
  min(res_std), max(res_std))
w("上限が ResAll の最小値を %s。

",
  if (nfe_upper < min(res_std)) "下回るので、NFE 側は ResAll と同等の効果量を排除できる"
  else "下回らないので、NFE 側は ResAll と同等の効果量を排除できない")
w("全予測変数の全項は outputs/table14_firstorder_full.csv。

")

w("### 11c. 信頼性指標

")
w("| 因子 | 指標 | 本データ | 第1論文公表値 | 一致 |
|---|---|---|---|---|
")
for (i in seq_len(nrow(rel_tab))) {
  r <- rel_tab[i, ]
  w("| %s | %s | %.3f | %.3f | %s |
", r$factor, r$statistic, r$value,
    r$paper1_published, if (r$matches_paper1) "OK" else "**不一致**")
}
w("
NFE は二値 2 項目の因子で、Spearman-Brown 係数は %.3f。

", rel_nfe_sb)

w("## 12. \u611f\u5ea6\u5206\u6790\uff08supplementary\uff09\n\n")
w("K6 \u5165\u308a\u306f\u672c\u6587 \u00a78 \u306b\u79fb\u3057\u305f\u306e\u3067\u3053\u3053\u306b\u306f\u306a\u3044\u3002\n\n")
w("| \u30af\u30e9\u30b9 | \u4e3b\uff08Model A\uff09| (ii) 215\u540d\u9664\u5916 | (iii) 3\u5024\u9806\u5e8fLCA | (iv) pp>=0.80 |\n|---|---|---|---|---|\n")
rS2 <- struct_rows(fitS2); rS3 <- struct_rows(fitS3); rS4a <- struct_rows(fitS4a)
for (k in 2:K) {
  a <- cls_row(rowsA, k); s2 <- cls_row(rS2, k)
  s3 <- rS3[rS3$rhs == paste0("vclassO", k), ]; s4 <- cls_row(rS4a, k)
  w("| %d | %s | %s | %s | %s |\n", k,
    fmt_ci(a$est, a$ci.lower, a$ci.upper), fmt_ci(s2$est, s2$ci.lower, s2$ci.upper),
    fmt_ci(s3$est, s3$ci.lower, s3$ci.upper), fmt_ci(s4$est, s4$ci.lower, s4$ci.upper))
}
w("\n**(ii) \u3068 A2a \u306f\u540c\u4e00\u306e\u63a8\u5b9a\u3067\u3042\u308b\u3002**215 \u540d\u306e\u96f6\u652f\u6301\u8005\u3092\u9664\u304f\u3068\u3044\u3046\u64cd\u4f5c\u306f\u3001")
w("\u652f\u6301\u9818\u57df\u304c 1 \u3064\u4ee5\u4e0a\u3042\u308b 1964 \u540d\u306b\u9650\u5b9a\u3059\u308b\u3053\u3068\u3068\u540c\u7fa9\u306a\u306e\u3067\u3001(ii) \u306f A2a \u305d\u306e\u3082\u306e\u3067\u3001")
w("\u72ec\u7acb\u3057\u305f\u691c\u8a3c\u3067\u306f\u306a\u3044\u3002\u4fc2\u6570\u306f %s\u3002",
  paste(sprintf("%.4f", sapply(2:K, function(k) cls_row(rowsA2a, k)$est)), collapse = " / "))
w("Model A\uff08\u5168 2179 \u540d\uff09\u306e %s \u3068\u307b\u3068\u3093\u3069\u5909\u308f\u3089\u306a\u3044\u306e\u3067\u3001",
  paste(sprintf("%.4f", sapply(2:K, function(k) cls_row(rowsA, k)$est)), collapse = " / "))
w("**\u96f6\u652f\u6301\u306e 215 \u540d\u304c\u30af\u30e9\u30b9 1 \u306e\u52b9\u679c\u3092\u99c6\u52d5\u3057\u3066\u3044\u308b\u308f\u3051\u3067\u306f\u306a\u3044**\u3068\u8a00\u3048\u308b\u3002")
w("\u6700\u5927\u5dee\u306f %.4f\u3002\n\n",
  max(abs(sapply(2:K, function(k) cls_row(rowsA2a, k)$est - cls_row(rowsA, k)$est))))
w("(iii) \u9806\u5e8f\u6307\u6a19\u3067\u652f\u6301\u304b\u3064\u91cd\u8981\u5ea6 1-2 \u306e\u30bb\u30eb\u306f %d \u500b\u3001\u30ec\u30d9\u30eb 2 \u306b\u542b\u3081\u305f\u3002", n_low_imp)
w("(iv) \u4fdd\u6301 %d / %d \u540d\u3002\u30af\u30e9\u30b9\u5225: %s\n\n", n_keep80, nrow(data),
  paste(sprintf("\u30af\u30e9\u30b9%d %d/%d", seq_len(K), keep80_by_class[, "TRUE"],
                rowSums(keep80_by_class)), collapse = "\u3001"))

w("## 13. supplementary: \u30af\u30e9\u30b9\u6570\u306e\u5217\u6319\u3068 K = 5 \u89e3\n\n")
w("| K | logLik | npar | AIC | BIC | aBIC | entropy | \u6700\u5c0f\u30af\u30e9\u30b9 n | \u6700\u5c0f\u30af\u30e9\u30b9 %% |\n")
w("|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(LC$fit))) {
  r <- LC$fit[i, ]
  w("| %d | %.1f | %d | %.1f | %.1f | %.1f | %s | %d | %.1f |\n",
    r$K, r$logLik, r$npar, r$AIC, r$BIC, r$aBIC,
    if (is.na(r$entropy)) "\u2014" else sprintf("%.3f", r$entropy),
    r$min_class_n, r$min_class_pct)
}
w("\nBIC \u6700\u5c0f\u306f K = %d\u3001aBIC \u6700\u5c0f\u306f K = %d\u3002\u4e71\u6570\u7a2e %d\u3001\u30b9\u30bf\u30fc\u30c8 %d \u56de\u3002\n\n",
  LC$fit$K[which.min(LC$fit$BIC)], LC$fit$K[which.min(LC$fit$aBIC)], LC$seed, LC$nrep)
w("### K = 5 \u306e\u30af\u30e9\u30b9\u5256\u9762\n\n")
w("| \u30af\u30e9\u30b9 | n | %% | \u5e73\u5747\u652f\u6301\u9818\u57df\u6570 | %s |\n", paste(DOM_LAB, collapse = " | "))
w("|---|---|---|---|%s\n", paste(rep("---|", 9), collapse = ""))
for (k in 1:K5) {
  w("| %d | %d | %.1f | %.2f | %s |\n", k, sum(m5$predclass == k),
    100 * sum(m5$predclass == k) / nrow(data), mean(d$nsup[m5$predclass == k]),
    paste(sprintf("%.2f", P5[k, ]), collapse = " | "))
}
w("\n**(a) BIC \u6700\u5c0f\u306f K = 5\uff08%.1f\uff09\u3001K = 4 \u306f %.1f\u3002**\n\n",
  min(LC$fit$BIC), LC$fit$BIC[LC$fit$K == K])
w("**(b) K = 5 \u306e\u7b2c5\u30af\u30e9\u30b9\u306f K = 4 \u89e3\u306e\u5206\u5272\u3067\u3042\u308b\u3002**\n\n")
w("| K=4 \\ K=5 | %s |\n|---|%s\n", paste(1:K5, collapse = " | "),
  paste(rep("---|", K5), collapse = ""))
for (k in seq_len(K)) w("| %d: %s | %s |\n", k, CLASS_LAB[k],
                        paste(xtab_45[k, ], collapse = " | "))
w("\n**(c) n = %d \u3067\u306f\u591a\u7fa4\u3068\u4ea4\u4e92\u4f5c\u7528\u3067\u30bb\u30eb\u304c\u500b\u4f4d\u6570\u306b\u306a\u308b\u3002**\n\n",
  min(as.integer(table(m5$predclass))))
w("- \u5bb6\u65cf\u306e\u304c\u3093\u7d4c\u9a13\u306a\u3057 %d \u540d\uff0f\u3042\u308a %d \u540d\n", cells5_fc[K5, 1], cells5_fc[K5, 2])
w("- \u88ab\u5bb3\u306a\u3057 %d \u540d\uff0f\u672a\u56de\u5fa9 %d \u540d\uff0f\u56de\u5fa9 %d \u540d\n\n",
  cells5_dis[K5, "unexposed"], cells5_dis[K5, "unrecovered"], cells5_dis[K5, "recovered"])

w("## 14. \u3053\u308c\u307e\u3067\u306b\u691c\u5b9a\u3057\u305f\u5bfe\u6bd4\u306e\u7dcf\u6570\n\n")
w("\u591a\u91cd\u6027\u306e\u8abf\u6574\u306f\u3057\u3066\u3044\u306a\u3044\u3002\u3053\u306e\u4e00\u89a7\u306f\u3001\u672c\u6587\u3067 confirmatory \u3068\u547c\u3079\u308b\u7bc4\u56f2\u3092")
w("\u610f\u56f3\u7684\u306b\u7dda\u5f15\u304d\u3059\u308b\u305f\u3081\u306e\u3082\u306e\u3067\u3042\u308b\u3002\n\n")
w("**\u691c\u5b9a\u7dcf\u6570 %d\u3001\u3046\u3061 p < 0.05 \u306f %d\u3002**\n\n", nrow(tests), n_sig)
w("| \u65cf | \u4e8b\u524d\u6307\u5b9a | \u5f8c\u304b\u3089\u8ffd\u52a0 | \u8a08 | p<0.05 |\n|---|---|---|---|---|\n")
for (fam in unique(tests$family)) {
  z <- tests[tests$family == fam, ]
  w("| %s | %d | %d | %d | %d |\n", fam,
    sum(z$status == "pre-specified"), sum(z$status != "pre-specified"),
    nrow(z), sum(z$p < 0.05, na.rm = TRUE))
}
w("| **\u5408\u8a08** | **%d** | **%d** | **%d** | **%d** |\n\n",
  sum(tests$status == "pre-specified"), sum(tests$status != "pre-specified"),
  nrow(tests), n_sig)
w("\u4e8b\u524d\u6307\u5b9a\u306f\u3001\u30af\u30e9\u30b9\u52b9\u679c\u306e\u540c\u6642\u691c\u5b9a\u3001\u5404\u30af\u30e9\u30b9\u5bfe\u53c2\u7167\u3001\u30b3\u30df\u30c3\u30c8\u30e1\u30f3\u30c8\u8abf\u6574\u3001")
w("\u7fa4\u9593\u6bd4\u8f03\u3001\u552f\u4e00\u306e\u4e8b\u524d\u6307\u5b9a\u4ea4\u4e92\u4f5c\u7528\u3002")
w("**\u5f8c\u304b\u3089\u8ffd\u52a0\u3057\u305f\u3082\u306e\uff08= \u672c\u6587\u3067\u63a2\u7d22\u7684\u3068\u660e\u793a\u3059\u3079\u304d\u7bc4\u56f2\uff09**\u306f\u3001")
w("\u5bfe\u30af\u30e9\u30b9\u5bfe\u6bd4\u3001b2 = b3 = b4 \u306e\u7b49\u5024\u691c\u5b9a\u3001\u4e00\u6b21\u56e0\u5b50\u5206\u89e3\u3001")
w("\u652f\u6301\u6570\u3092\u63c3\u3048\u305f\u5bfe\u6bd4\uff084c\uff09\u3001nsup \u8abf\u6574\u5f8c\u306e\u30af\u30e9\u30b9\u52b9\u679c\uff084d\uff09\u3002\n\n")
w("**1 \u9818\u57df\u3042\u305f\u308a\u306e\u50be\u304d\uff08\u00a75\uff09\u306f\u65b0\u305f\u306a\u691c\u5b9a\u3092\u884c\u3063\u3066\u3044\u306a\u3044**\u306e\u3067\u3001\u4e0a\u306e\u7dcf\u6570\u306b\u306f\u5165\u308c\u3066\u3044\u306a\u3044\u3002")
w("3 \u884c\u306f\u3044\u305a\u308c\u3082\u65e2\u51fa\u306e\u5bfe\u6bd4\uff08%s\uff09\u3092\u5b9a\u6570\u3067\u5272\u3063\u305f\u3082\u306e\u3067\u3042\u308a\u3001",
  paste(slope_dup, collapse = "\u3001"))
w("\u5225\u52d8\u5b9a\u306b\u3059\u308b\u3068\u540c\u3058\u691c\u5b9a\u3092\u4e8c\u91cd\u306b\u6570\u3048\u308b\u3053\u3068\u306b\u306a\u308b\u3002\n\n")
w("\u5168\u691c\u5b9a\u306e\u4e00\u89a7\u306f outputs/table9_test_inventory.csv\u3002\n\n")

w("## \u51fa\u5178\n\n")
w("\u5168\u6570\u5024\u306f `02_structural.R` \u304c\u751f\u6210\u3002\u8868\u306e\u5143\u30c7\u30fc\u30bf\u306f outputs/ \u306e ")
w("table1_background.csv\u3001table2_class_profile.csv\u3001table3_modelA.csv\u3001table4_modelB.csv\u3001")
w("table5_k6_contrast.csv\u3001table6_firstorder.csv\u3001table7_slope_per_domain.csv\u3001")
w("contrasts_class2_vs_3.csv\u3001sensitivity.csv\u3002\u30e2\u30c7\u30eb\u306f structural_fits.rds\u3002\n")
close(con)

saveRDS(list(K = K, labels = CLASS_LAB, fitA = fitA, fitA_free = fitA_free,
             fitA2a = fitA2a, fitA2b = fitA2b, fitB = fitB, fitI = fitI,
             fitK6 = fitK6, fitFO = fitFO, fitS3 = fitS3, fitS4a = fitS4a,
             fitS4b = fitS4b, lca_ord = lca_ord, vclass = data$vclass,
             pp_max = data$pp_max, entropy = ent_K, avepp = avepp,
             xtab_45 = xtab_45, wb_var_total = wb_var_total,
             wb_var_res = wb_var_res, con23 = con23, slope = slope,
             k6_tab = k6_tab, fo_tab = fo_tab,
             table1 = table1, table2 = t2, table3 = table3, table4 = table4,
             sens = sens, cells = cells, vif_main = vif_main, vif_ix = vif_ix),
        "structural_fits.rds")

cat(sprintf("done in %.1f s -- outputs/results_phase2.md\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

# =============================================================================
# 06_breadth_group_ratios.R -- covariate coefficients of the BREADTH multiple-
# group model (fit_ncat_g in 02_structural.R), by family cancer history.
#
# Adds no method. It refits the model that already produces Table 4 and reads
# off the covariate coefficients that 02_structural.R did not write out, so the
# Discussion can compare the amplification of the breadth terms against the
# amplification of the covariates within the SAME model.
#
#   output  outputs/table18_breadth_groupratio.csv   (new; touches nothing else)
# =============================================================================
suppressPackageStartupMessages({ library(lavaan); library(dplyr) })
options(stringsAsFactors = FALSE)
BASE <- "."  # the repository root; run every script from there
setwd(BASE)

d  <- readRDS("values_analysis.rds")
LC <- readRDS("lca_models.rds")
K  <- 4L
m  <- LC$models[[K]]

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
MEAS  <- paste("NFE =~ sadR + worryR", "ResAll =~ ls + gh + es + vital + comp",
               "WB =~ NFE + ResAll", "ls ~~ gh", "es ~~ vital",
               "ResAll ~~ 0*ResAll", sep = "\n")

data$nsupcat <- cut(pmin(data$nsup, 6), breaks = -1:6,
                    labels = c("0","1","2","3","4","5","6+"))
for (lv in 1:6) data[[paste0("nsc", lv)]] <- as.numeric(pmin(data$nsup, 6) == lv)
NSC_ALL <- paste0("nsc", 1:6)

M_ncat_g <- paste0(MEAS, "\nWB ~ ",
                   paste(sprintf("c(h0_%s, h1_%s)*%s", sub("nsc", "", NSC_ALL),
                                 sub("nsc", "", NSC_ALL), NSC_ALL), collapse = " + "),
                   " + ", COV_B, "\n")
fit <- sem(M_ncat_g, data, group = "fcancer2", estimator = "MLR", missing = "ML")

TMAP <- c(nsc1 = "breadth_1", nsc2 = "breadth_2", nsc3 = "breadth_3",
          nsc4 = "breadth_4", nsc5 = "breadth_5", nsc6 = "breadth_6plus",
          disaster1 = "disaster_unrecovered", disaster2 = "disaster_recovered",
          female1 = "female", age_z = "age_z", income_z = "income_z",
          edu22 = "edu_beyond_highschool", edu23 = "edu_declined",
          occupation22 = "occupation_employee", occupation23 = "occupation_other",
          marriage2 = "married", marriage3 = "divorced_widowed",
          comobidity1 = "comorbidity", fukushima1 = "fukushima_resident",
          health22 = "subjective_health_high")

pe <- parameterEstimates(fit, standardized = TRUE)
pe <- pe[pe$op == "~" & pe$lhs == "WB", ]
g0 <- pe[pe$group == 1, ]; g1 <- pe[pe$group == 2, ]
stopifnot(identical(g0$rhs, g1$rhs))
out <- data.frame(
  term        = unname(TMAP[g0$rhs]),
  b_nocancer  = round(g0$est, 4),   b_cancer  = round(g1$est, 4),
  b_ratio     = round(g1$est / g0$est, 2),
  beta_nocancer = round(g0$std.all, 3), beta_cancer = round(g1$std.all, 3),
  beta_ratio  = round(g1$std.all / g0$std.all, 2),
  p_nocancer  = signif(g0$pvalue, 3), p_cancer = signif(g1$pvalue, 3))
write.csv(out, "outputs/table18_breadth_groupratio.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
print(out, row.names = FALSE)

# within-group cell counts at each level of breadth, which table12 reports as
# overall counts for the two Model B blocks
print(table(data$nsupcat, data$fcancer2))

# ---------------------------------------------------------------------------
# The first-order model with breadth as the exposure (fitFO_ns in
# 02_structural.R). Two things that script computed but did not write out:
#   the standardised solution with its own confidence intervals, which Fig. 2a
#     needs, because rescaling an unstandardised interval by a point ratio is
#     not the same interval; and
#   the two joint Wald tests the Results reports, which the test inventory
#     lists in its class form only.
#
#   output  outputs/table20_firstorder_breadth_std.csv
#           outputs/table21_firstorder_breadth_tests.csv
# ---------------------------------------------------------------------------
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

ss <- standardizedSolution(fitFO_ns)
ss <- ss[ss$op == "~" & ss$lhs %in% c("NFE", "ResAll"), ]
LVLNAME <- c("1", "2", "3", "4", "5", "6+")
is_nsc <- grepl("^nsc[1-6]$", ss$rhs)
TMAP2 <- c(TMAP, fcancer21 = "family_cancer",
           mcom_z = "mean_commitment_z", k6_z = "k6_z")
std_out <- data.frame(
  factor = ifelse(ss$lhs == "NFE", "NFE (negative-free emotion)",
                  "ResAll (resource-oriented well-being)"),
  term = ifelse(ss$rhs %in% names(TMAP2), unname(TMAP2[ss$rhs]), ss$rhs),
  level = ifelse(is_nsc, LVLNAME[suppressWarnings(
    as.integer(sub("nsc", "", ss$rhs)))], ""),
  beta = round(ss$est.std, 4), se = round(ss$se, 4),
  ci_low = round(ss$ci.lower, 4), ci_high = round(ss$ci.upper, 4),
  p = signif(ss$pvalue, 4))
write.csv(std_out, "outputs/table20_firstorder_breadth_std.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

wn <- lavTestWald(fitFO_ns, constraints = paste(sprintf("fn%d == 0", 1:6),
                                                collapse = "\n"))
wr <- lavTestWald(fitFO_ns, constraints = paste(sprintf("fr%d == 0", 1:6),
                                                collapse = "\n"))
tests <- data.frame(
  family = "first-order, breadth as exposure",
  test = c("breadth on NFE, all levels = 0",
           "breadth on ResAll, all levels = 0"),
  statistic = round(c(wn$stat, wr$stat), 3),
  df = c(wn$df, wr$df),
  p = signif(c(wn$p.value, wr$p.value), 4),
  status = "added 2026-08-30")
write.csv(tests, "outputs/table21_firstorder_breadth_tests.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
print(std_out, row.names = FALSE)
print(tests, row.names = FALSE)

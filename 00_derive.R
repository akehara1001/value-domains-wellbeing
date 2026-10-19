# =============================================================================
# 00_derive.R -- analysis dataset for the VALUE DOMAINS x well-being paper
# =============================================================================
# Second paper from a 2019 web survey of adults in Japan (N=2179).
# The exposure side is the Q11 value-domain block; the well-being measurement
# model is carried over unchanged from paper 1, the companion analysis.
#
#   input   data/data.xlsx  sheet "data"    2179 x 177, raw item codes
#           data/data.xlsx  sheet "Format"  question and value labels
#   output  values_analysis.rds / .csv   analysis dataset, KEY retained
#           outputs/gate_00_derive.txt   verification gate log
#
# sem.xlsx is NEVER read: no script generating it exists, so its provenance
# cannot be audited. Everything here is rebuilt from the raw export.
#
# The covariate recodes reproduce the paper-1 dataset script, which was
# itself validated row-by-row against sem.xlsx, so that the paper-1 measurement
# model refits to the published fit indices. The value-domain variables are new.
#
# VERIFICATION GATE at the bottom: sample counts and moments from paper 1
# Table 1, the value-block distributions, and the measurement-model fit.
# Any failure aborts before anything is written.
#
# Run:  Rscript 00_derive.R
# =============================================================================

suppressPackageStartupMessages({library(readxl); library(dplyr); library(lavaan)})
options(stringsAsFactors = FALSE)

BASE   <- "."     # the repository root; run every script from there
RAWDIR <- "data"  # data/data.xlsx, the raw survey export
setwd(BASE)
dir.create("outputs", showWarnings = FALSE)

t_start <- Sys.time()
GATE <- file("outputs/gate_00_derive.txt", open = "wt", encoding = "UTF-8")
say <- function(...) { txt <- sprintf(...); cat(txt); cat(txt, file = GATE) }

raw <- as.data.frame(read_excel(file.path(RAWDIR, "data.xlsx"), sheet = "data"))
stopifnot(nrow(raw) == 2179, ncol(raw) == 177)
q <- function(nm) raw[[nm]]

# -----------------------------------------------------------------------------
# 1. Value domains (Q11 block) -- the new exposure
# -----------------------------------------------------------------------------
# Q11_1[1..9]  multi-select, "which value domains matter in your life"; 0/1 with
#              no missingness. Q11_1[10] "none of these" is the exact complement
#              of endorsing nothing (checked in the gate), so it carries no
#              information beyond the nine indicators and is not an indicator.
# Q11_3..Q11_6 importance / commitment / intention to advance / consistency over
#              the past 2-3 months. ROUTED: asked only about endorsed domains,
#              hence NA for every non-endorsed domain and for all nine among the
#              215 who endorsed none. Never used as LCA indicators.
# Q11_7..Q11_9 retrospective change around the family member's cancer diagnosis;
#              asked only of respondents with a family cancer history, and only
#              about endorsed domains.
DOM <- c("family", "friends", "romance", "work", "education",
         "leisure", "spiritual", "community", "health")
DOM_JA <- c("\u5bb6\u65cf\u95a2\u4fc2",
            "\u53cb\u4eba\uff0f\u793e\u4f1a\u95a2\u4fc2",
            "\u604b\u4eba\uff0f\u604b\u611b\u95a2\u4fc2",
            "\u4ed5\u4e8b\uff0f\u30ad\u30e3\u30ea\u30a2",
            "\u6559\u80b2\uff0f\u500b\u4eba\u7684\u6210\u9577",
            "\u30ec\u30af\u30ea\u30a8\u30fc\u30b7\u30e7\u30f3\uff0f\u30ec\u30b8\u30e3\u30fc\uff0f\u30b9\u30dd\u30fc\u30c4",
            "\u30b9\u30d4\u30ea\u30c1\u30e5\u30a2\u30ea\u30c6\u30a3\uff0f\u5b97\u6559",
            "\u5730\u57df\u6027(\u30b3\u30df\u30e5\u30cb\u30c6\u30a3)\uff0f\u56fd\u6c11\u6027",
            "\u5065\u5eb7\uff0f\u8eab\u4f53\u7684\u30a6\u30a7\u30eb\u30d3\u30fc\u30a4\u30f3\u30b0")

sup <- sapply(1:9, function(j) as.integer(q(sprintf("Q11_1[%d]", j))))
colnames(sup) <- paste0("v_", DOM)
none <- as.integer(q("Q11_1[10]"))

# the routed 1-5 blocks, NA preserved
grab <- function(prefix) {
  m <- sapply(1:9, function(j) as.numeric(q(sprintf("%s_%d", prefix, j))))
  colnames(m) <- paste0(prefix, "_", DOM)
  m
}
imp     <- grab("Q11_3")   # importance,  1 not at all .. 5 extremely
com     <- grab("Q11_4")   # commitment,  1 not at all .. 5 extremely
int     <- grab("Q11_5")   # intention,   1 not at all .. 5 extremely
con     <- grab("Q11_6")   # consistency, 1 = 0-20% .. 5 = 81-100%   BANDED
chg_imp <- grab("Q11_7")   # change in importance,  3 = unchanged
chg_com <- grab("Q11_8")   # change in commitment,  3 = unchanged
chg_int <- grab("Q11_9")   # change in intention,   3 = unchanged

# Importance (Q11_3, a 5-point Likert) and consistency (Q11_6, five 20-point
# percentage bands) are not on a common scale. Their raw difference is not a
# value-behaviour discrepancy and is deliberately never computed anywhere.

# Person-level means across ENDORSED domains only, NA for the 215 who endorsed
# nothing. rowMeans(na.rm = TRUE) on an all-NA row returns NaN, so those rows
# are set to NA explicitly.
rowmean_endorsed <- function(m) {
  ifelse(rowSums(!is.na(m)) > 0, rowMeans(m, na.rm = TRUE), NA_real_)
}

# -----------------------------------------------------------------------------
# 2. Well-being items and covariates -- coding identical to paper 1
# -----------------------------------------------------------------------------
d <- tibble(
  KEY = raw$KEY,
  row = seq_len(nrow(raw)),

  ## demographics
  age    = q("Q0_1"),
  female = as.integer(q("Q0_2") == 2),               # Q0_2: 1 male, 2 female
  income = q("QW15"),                                # 14 ordered bands
  # QW3 is the panel-supplied prefecture code; 7 = Fukushima. This is residence
  # in 2019, not in 2011. Equals QBD in {1,2} on all 2179 rows (checked).
  fukushima = as.integer(q("QW3") == 7),

  ## family cancer history. Q0_3 1-6 = relationship of the diagnosed member,
  ## 7 = none. Code 8, "prefer not to say", does not occur in these data.
  fcancer2 = as.integer(q("Q0_3") != 7),

  ## FAMCARE-10 satisfaction sum; asked only of the 855 with a family history
  sf = rowSums(sapply(paste0("Q2_", 1:10), function(v) q(v) - 1)),

  ## K6: items shifted to 0-4 and summed. Q9_7 (impairment) is not part of K6.
  k6 = rowSums(sapply(paste0("Q9_", 1:6), function(v) q(v) - 1)),

  ## well-being items -- the seven in the paper-1 measurement model.
  ## Q4/Q5 are stored REVERSED (code 1 = 10 points ... code 11 = 0 points), so
  ## 11 - code restores the intended 0-10 metric.
  ls = 11 - q("Q4"),                                 # life satisfaction, 0-10
  gh = 11 - q("Q5"),                                 # general happiness, 0-10
  # Q6_5 sadness and Q6_6 worry are 1 = yes, 2 = no. Reversed to ABSENCE
  # (1 = free of the state) so all seven items point the same way.
  sadR   = 1 - ifelse(q("Q6_5") == 2, 0, q("Q6_5")),
  worryR = 1 - ifelse(q("Q6_6") == 2, 0, q("Q6_6")),
  vital = q("Q7_1_1"),                               # 1-4, full of vitality
  es    = q("Q7_1_2"),                               # 1-4, calm and peaceful
  comp  = q("Q7_2_7"),                               # 1-5, sense of accomplishment

  ## covariates, coding conventions carried over from paper 1 verbatim
  health2     = ifelse(q("Q8") %in% c(1, 2), 2, 1),  # 2 = high, 1 = low
  comobidity  = ifelse(q("Q10[11]") == 1, 0, 1),     # Q10[11] = none of these
  occupation2 = ifelse(q("Q12_1") %in% 1:5, 2,
                ifelse(q("Q12_1") %in% c(6, 7), 1, 3)),  # 1 = self-employed, ref
  marriage    = ifelse(q("Q12_2") == 1, 1,
                ifelse(q("Q12_2") %in% c(2, 3), 2, 3)),  # 1 = never married, ref
  edu2        = ifelse(q("Q12_4") %in% c(1, 2), 1,
                ifelse(q("Q12_4") %in% 3:8, 2, 3))       # 1 = junior/high, ref
)

## Disaster recovery status. Q3_1 in {1,2} = no damage, so unexposed. Otherwise
## Q3_2 in {1,2} -> recovered, Q3_2 in {3,4} -> unrecovered. Paper 1 codes
## 0 = unexposed, 1 = UNRECOVERED (n = 95), 2 = RECOVERED (n = 372). That
## ordering is counter-intuitive but is what every published table uses.
d_ex2 <- ifelse(q("Q3_1") %in% c(1, 2), 0, 1)
d_re2 <- ifelse(q("Q3_2") %in% c(1, 2), 1, ifelse(q("Q3_2") %in% c(3, 4), 0, NA))
d$disaster <- ifelse(d_ex2 == 0, 0, ifelse(d_re2 == 1, 2, 1))

## Income collapsed for the descriptive table only; models use the 1-14 band.
d$income3 <- ifelse(d$income %in% 1:3, 1, ifelse(d$income %in% 4:6, 2, 3))

# -----------------------------------------------------------------------------
# 3. Bind the value block on
# -----------------------------------------------------------------------------
d <- bind_cols(d,
  as.data.frame(sup),
  tibble(v_none = none, nsup = as.integer(rowSums(sup))),
  as.data.frame(imp), as.data.frame(com), as.data.frame(int), as.data.frame(con),
  as.data.frame(chg_imp), as.data.frame(chg_com), as.data.frame(chg_int),
  tibble(m_importance  = rowmean_endorsed(imp),
         m_commitment  = rowmean_endorsed(com),
         m_intention   = rowmean_endorsed(int),
         m_consistency = rowmean_endorsed(con)))
d <- as.data.frame(d)

# =============================================================================
# VERIFICATION GATE
# =============================================================================
# Every expected value below comes from paper 1 Table 1 / Table 2 or from the
# instruction sheet. One failure aborts: no downstream script may run on a
# dataset that does not reproduce the published sample description.
say("=============================================================\n")
say(" 00_derive.R verification gate   %s\n", format(t_start, "%Y-%m-%d %H:%M:%S"))
say("=============================================================\n\n")

fails <- 0L
chk <- function(label, obs, exp, tol = 0) {
  ok <- !is.na(obs) && abs(obs - exp) <= tol
  if (!ok) fails <<- fails + 1L
  say("%-48s obs %9.4f  exp %9.4f  %s\n", label, obs, exp,
      if (ok) "OK" else "<<< FAIL")
  invisible(ok)
}
chkv <- function(label, obs, exp) {
  ok <- length(obs) == length(exp) && all(obs == exp)
  if (!ok) fails <<- fails + 1L
  say("%-48s %s\n", label, if (ok) "OK" else "<<< FAIL")
  if (!ok) {
    say("    obs: %s\n", paste(obs, collapse = " "))
    say("    exp: %s\n", paste(exp, collapse = " "))
  }
  invisible(ok)
}

say("-- sample and covariates ------------------------------------\n")
chk("N",                               nrow(d), 2179)
chk("family cancer: yes",              sum(d$fcancer2 == 1), 855)
chk("family cancer: no",               sum(d$fcancer2 == 0), 1324)
chk("Fukushima resident",              sum(d$fukushima == 1), 1077)
chk("outside Fukushima",               sum(d$fukushima == 0), 1102)
chk("disaster: unexposed",             sum(d$disaster == 0), 1712)
chk("disaster: recovered",             sum(d$disaster == 2), 372)
chk("disaster: unrecovered",           sum(d$disaster == 1), 95)
chk("  Q3_2 = mostly recovered",       sum(q("Q3_2") == 1, na.rm = TRUE), 306)
chk("  Q3_2 = recovery foreseeable",   sum(q("Q3_2") == 2, na.rm = TRUE), 66)
chk("  Q3_2 = not foreseeable",        sum(q("Q3_2") == 3, na.rm = TRUE), 35)
chk("  Q3_2 = no choice but to accept", sum(q("Q3_2") == 4, na.rm = TRUE), 60)
chkv("income, 3 bands",                as.integer(table(d$income3)), c(608, 841, 730))
chkv("education",                      as.integer(table(d$edu2)), c(808, 1347, 24))
chkv("marital status",                 as.integer(table(d$marriage)), c(686, 1309, 184))
chkv("occupation self/employee/other", as.integer(table(d$occupation2)), c(159, 1003, 1017))
chk("subjective health: high",         sum(d$health2 == 2), 803)
chk("subjective health: low",          sum(d$health2 == 1), 1376)
chk("comorbidity: yes",                sum(d$comobidity == 1), 1006)
chk("comorbidity: no",                 sum(d$comobidity == 0), 1173)
chk("sex: male",                       sum(d$female == 0), 1203)
chk("sex: female",                     sum(d$female == 1), 976)
chk("age mean",                        mean(d$age), 46.7, 0.05)
chk("age SD",                          sd(d$age), 13.4, 0.05)
chk("age min",                         min(d$age), 20)
chk("age max",                         max(d$age), 69)
chk("K6 mean",                         mean(d$k6), 6.5, 0.05)
chk("K6 SD",                           sd(d$k6), 5.7, 0.05)
chk("K6 min",                          min(d$k6), 0)
chk("K6 max",                          max(d$k6), 24)

say("\n-- value domains --------------------------------------------\n")
chkv("domains endorsed, distribution 0..9",
     as.integer(table(factor(d$nsup, levels = 0:9))),
     c(215, 367, 391, 398, 302, 208, 127, 113, 48, 10))
exp_rate <- c(77.6, 53.4, 21.6, 35.3, 21.3, 24.4, 4.1, 12.6, 53.3)
for (j in 1:9) {
  chk(sprintf("endorsement pct -- %s", DOM[j]),
      100 * mean(sup[, j]), exp_rate[j], 0.05)
}
chk("Q11_1[10] 'none of these'",       sum(none == 1), 215)
chk("  disagreements with nsup == 0",  sum((none == 1) != (d$nsup == 0)), 0)
routed_ok <- all(sapply(1:9, function(j)
  all(is.na(imp[, j]) == (sup[, j] == 0)) &&
  all(is.na(com[, j]) == (sup[, j] == 0)) &&
  all(is.na(int[, j]) == (sup[, j] == 0)) &&
  all(is.na(con[, j]) == (sup[, j] == 0))))
if (!routed_ok) fails <- fails + 1L
say("%-48s %s\n", "Q11_3..6 observed iff domain endorsed",
    if (routed_ok) "OK" else "<<< FAIL")
chk("mean commitment defined for n",   sum(!is.na(d$m_commitment)), 1964)
chkv("Q11_7 observed n by domain",
     as.integer(colSums(!is.na(chg_imp))),
     c(677, 484, 185, 322, 185, 223, 37, 112, 456))

say("\n-- well-being items -----------------------------------------\n")
chk("ls min", min(d$ls), 0); chk("ls max", max(d$ls), 10)
chk("gh min", min(d$gh), 0); chk("gh max", max(d$gh), 10)
# Item means by family cancer status, paper 1 Table 2. The comparison is made at
# the precision the published table actually carries: the two binary items are
# printed to 3 decimals, the five continuous ones to 2. Comparing a 2-decimal
# published figure against an unrounded mean at a 4-decimal tolerance would fail
# on display precision alone, so the observed mean is rounded to the published
# number of decimals and must then match exactly.
T2 <- list(sadR   = c(0.714, 0.673, 3), worryR = c(0.477, 0.426, 3),
           ls     = c(5.88,  5.91,  2), gh     = c(6.17,  6.21,  2),
           es     = c(2.41,  2.39,  2), vital  = c(2.12,  2.15,  2),
           comp   = c(2.84,  2.89,  2))
for (v in names(T2)) {
  dp <- T2[[v]][3]
  chk(sprintf("%-6s mean, no family cancer", v),
      round(mean(d[[v]][d$fcancer2 == 0]), dp), T2[[v]][1])
  chk(sprintf("%-6s mean, family cancer", v),
      round(mean(d[[v]][d$fcancer2 == 1]), dp), T2[[v]][2])
}

# -----------------------------------------------------------------------------
# Measurement model refit -- must reproduce paper 1's published fit
# -----------------------------------------------------------------------------
# This is fit_main of the paper-1 measurement script: the 7-item, two-first-
# order-factor, second-order model, MULTIGROUP by family cancer status, MLR,
# missing = "ML", WITH paper 1's structural regression (K6 included). The three
# published fit indices are produced by that whole specification and not by a
# bare single-group CFA, so the whole specification is what is refitted here.
# K6 appears in this gate for that reason only; the models of paper 2 exclude it.
say("\n-- paper 1 measurement model refit ---------------------------\n")
md <- d
md$age    <- as.numeric(scale(md$age))
md$income <- as.numeric(scale(md$income))
md$k6     <- as.numeric(scale(md$k6))
md <- md |> mutate(across(c(disaster, female, edu2, occupation2, marriage,
                            comobidity, fukushima, health2), factor))
dum <- model.matrix(~ disaster + female + edu2 + occupation2 + marriage +
                      comobidity + fukushima + health2, md)[, -1]
md <- cbind(md, dum)
STRUCT1 <- paste("disaster1+disaster2+female1 + k6 + age + income + edu22+edu23",
                 "+ occupation22 +occupation23+ marriage2 +marriage3+ comobidity1",
                 "+ fukushima1+health22")
MODEL7 <- paste0("
  NFE    =~ sadR + worryR
  ResAll =~ ls + gh + es + vital + comp
  WB     =~ NFE + ResAll
  ls ~~ gh
  es ~~ vital
  WB ~ ", STRUCT1, "\n")
fit1 <- sem(MODEL7, md, group = "fcancer2", estimator = "MLR", missing = "ML")
fm <- fitMeasures(fit1, c("chisq", "df", "cfi.robust", "rmsea.robust", "srmr"))
say("chi-square %.1f on %d df\n", fm[["chisq"]], fm[["df"]])
chk("robust CFI",   fm[["cfi.robust"]],   0.952, 0.0005)
chk("robust RMSEA", fm[["rmsea.robust"]], 0.044, 0.0005)
chk("SRMR",         fm[["srmr"]],         0.026, 0.0005)

say("\n=============================================================\n")
if (fails > 0L) {
  say(" GATE FAILED: %d check(s) did not match. Nothing was written.\n", fails)
  say("=============================================================\n")
  close(GATE)
  stop(sprintf("verification gate failed (%d checks)", fails))
}
say(" GATE PASSED: every check matches. Dataset written.\n")

saveRDS(d, "values_analysis.rds")
write.csv(d, "values_analysis.csv", row.names = FALSE, fileEncoding = "UTF-8")

cb <- data.frame(domain_index = 1:9, name = DOM, label_ja = DOM_JA,
                 item = sprintf("Q11_1[%d]", 1:9),
                 n_endorsed = as.integer(colSums(sup)),
                 pct_endorsed = round(100 * colMeans(sup), 1))
write.csv(cb, "outputs/domain_codebook.csv", row.names = FALSE, fileEncoding = "UTF-8")

say(" values_analysis.rds  %d rows x %d cols\n", nrow(d), ncol(d))
say(" elapsed %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs")))
say("=============================================================\n")
close(GATE)

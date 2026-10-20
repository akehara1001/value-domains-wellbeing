# =============================================================================
# 00_load_release.R -- analysis dataset from the released data file
# =============================================================================
# Builds values_analysis.rds, the dataset every later script reads, from the
# released file data/data_release_v53.csv (2179 rows, 92 columns, described in
# codebook_v53.csv). The released file already carries the coding the models
# use; this script only adds the variables derived from it and checks the
# sample against the published description.
#
#   input   data/data_release_v53.csv   not in this repository; see README
#   output  values_analysis.rds         analysis dataset
#           outputs/gate_00_load_release.txt   verification gate log
#           outputs/domain_codebook.csv  endorsement by domain, as 00_derive.R
#
# 00_derive.R builds the same dataset from the raw survey export, which is not
# public. Run one of the two, not both.
#
# Derived here, exactly as in 00_derive.R:
#   row      = id (the serial number, in the row order of the survey export)
#   nsup     = number of endorsed domains, the sum of the nine v_ indicators
#   v_none   = "none of these"; equal to nsup == 0 on every row of the raw
#              export (checked in 00_derive.R), so it is rebuilt from nsup
#   income3  = income in three bands, for the descriptive table only
#   m_importance, m_commitment, m_intention, m_consistency
#            = person-level means across the endorsed domains only, NA for
#              respondents who endorsed none
# The FAMCARE-10 satisfaction score (sf) is not released. 03_value_change.R
# skips its FAMCARE tabulation when sf is absent; no article table uses it.
#
# Run:  Rscript 00_load_release.R
# =============================================================================

suppressPackageStartupMessages({library(dplyr); library(lavaan)})
options(stringsAsFactors = FALSE)

BASE    <- "."     # the repository root; run every script from there
RELEASE <- file.path("data", "data_release_v53.csv")
setwd(BASE)
source("round_half_up.R")  # round() and signif() round half up from here on (see that file)
dir.create("outputs", showWarnings = FALSE)

t_start <- Sys.time()
GATE <- file("outputs/gate_00_load_release.txt", open = "wt", encoding = "UTF-8")
say <- function(...) { txt <- sprintf(...); cat(txt); cat(txt, file = GATE) }

if (!file.exists(RELEASE))
  stop("data/data_release_v53.csv not found; see the Data section of README.md")
rel <- read.csv(RELEASE, na.strings = "", fileEncoding = "UTF-8")

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
COLS <- c("id", "age", "female", "income", "edu2", "occupation2", "marriage",
          "fukushima", "health2", "comobidity", "fcancer2", "disaster", "k6",
          "ls", "gh", "sadR", "worryR", "vital", "es", "comp",
          paste0("v_", DOM),
          unlist(lapply(c("Q11_3", "Q11_4", "Q11_5", "Q11_6",
                          "Q11_7", "Q11_8", "Q11_9"),
                        function(p) paste0(p, "_", DOM))))
stopifnot(identical(names(rel), COLS), nrow(rel) == 2179,
          identical(rel$id, seq_len(2179)))

blk <- function(prefix) {
  m <- as.matrix(rel[, paste0(prefix, "_", DOM)])
  storage.mode(m) <- "double"
  m
}
rowmean_endorsed <- function(m) {
  ifelse(rowSums(!is.na(m)) > 0, rowMeans(m, na.rm = TRUE), NA_real_)
}
sup <- as.matrix(rel[, paste0("v_", DOM)])
storage.mode(sup) <- "integer"
imp <- blk("Q11_3"); com <- blk("Q11_4"); int <- blk("Q11_5"); con <- blk("Q11_6")
nsup <- as.integer(rowSums(sup))

# same columns, order and types as the dataset written by 00_derive.R, less the
# panel key and sf
num <- function(v) as.numeric(rel[[v]])
d <- data.frame(
  row = rel$id,
  age = num("age"), female = as.integer(rel$female), income = num("income"),
  fukushima = as.integer(rel$fukushima), fcancer2 = as.integer(rel$fcancer2),
  k6 = num("k6"), ls = num("ls"), gh = num("gh"),
  sadR = num("sadR"), worryR = num("worryR"),
  vital = num("vital"), es = num("es"), comp = num("comp"),
  health2 = num("health2"), comobidity = num("comobidity"),
  occupation2 = num("occupation2"), marriage = num("marriage"),
  edu2 = num("edu2"), disaster = num("disaster"))
d$income3 <- ifelse(d$income %in% 1:3, 1, ifelse(d$income %in% 4:6, 2, 3))
d <- cbind(d, as.data.frame(sup),
           data.frame(v_none = as.integer(nsup == 0), nsup = nsup),
           as.data.frame(imp), as.data.frame(com), as.data.frame(int),
           as.data.frame(con),
           as.data.frame(blk("Q11_7")), as.data.frame(blk("Q11_8")),
           as.data.frame(blk("Q11_9")),
           data.frame(m_importance  = rowmean_endorsed(imp),
                      m_commitment  = rowmean_endorsed(com),
                      m_intention   = rowmean_endorsed(int),
                      m_consistency = rowmean_endorsed(con)))

# =============================================================================
# VERIFICATION GATE -- the checks of 00_derive.R that the released file allows
# =============================================================================
say("=============================================================\n")
say(" 00_load_release.R verification gate   %s\n", format(t_start, "%Y-%m-%d %H:%M:%S"))
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
  invisible(ok)
}

say("-- sample and covariates ------------------------------------\n")
chk("N",                               nrow(d), 2179)
chk("family cancer: yes",              sum(d$fcancer2 == 1), 855)
chk("Fukushima resident",              sum(d$fukushima == 1), 1077)
chk("disaster: unexposed",             sum(d$disaster == 0), 1712)
chk("disaster: recovered",             sum(d$disaster == 2), 372)
chk("disaster: unrecovered",           sum(d$disaster == 1), 95)
chkv("income, 3 bands",                as.integer(table(d$income3)), c(608, 841, 730))
chkv("education",                      as.integer(table(d$edu2)), c(808, 1347, 24))
chkv("marital status",                 as.integer(table(d$marriage)), c(686, 1309, 184))
chkv("occupation self/employee/other", as.integer(table(d$occupation2)), c(159, 1003, 1017))
chk("subjective health: high",         sum(d$health2 == 2), 803)
chk("comorbidity: yes",                sum(d$comobidity == 1), 1006)
chk("sex: female",                     sum(d$female == 1), 976)
chk("age mean",                        mean(d$age), 46.7, 0.05)
chk("age SD",                          sd(d$age), 13.4, 0.05)
chk("K6 mean",                         mean(d$k6), 6.5, 0.05)
chk("K6 SD",                           sd(d$k6), 5.7, 0.05)

say("\n-- value domains --------------------------------------------\n")
chkv("domains endorsed, distribution 0..9",
     as.integer(table(factor(d$nsup, levels = 0:9))),
     c(215, 367, 391, 398, 302, 208, 127, 113, 48, 10))
routed_ok <- all(sapply(1:9, function(j)
  all(is.na(imp[, j]) == (sup[, j] == 0)) && all(is.na(com[, j]) == (sup[, j] == 0)) &&
  all(is.na(int[, j]) == (sup[, j] == 0)) && all(is.na(con[, j]) == (sup[, j] == 0))))
if (!routed_ok) fails <- fails + 1L
say("%-48s %s\n", "Q11_3..6 observed iff domain endorsed",
    if (routed_ok) "OK" else "<<< FAIL")
chk("mean commitment defined for n",   sum(!is.na(d$m_commitment)), 1964)
chkv("Q11_7 observed n by domain",
     as.integer(colSums(!is.na(blk("Q11_7")))),
     c(677, 484, 185, 322, 185, 223, 37, 112, 456))

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

cb <- data.frame(domain_index = 1:9, name = DOM, label_ja = DOM_JA,
                 item = sprintf("Q11_1[%d]", 1:9),
                 n_endorsed = as.integer(colSums(sup)),
                 pct_endorsed = round(100 * colMeans(sup), 1))
write.csv(cb, "outputs/domain_codebook.csv", row.names = FALSE, fileEncoding = "UTF-8")

say(" values_analysis.rds  %d rows x %d cols\n", nrow(d), ncol(d))
say(" elapsed %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs")))
say("=============================================================\n")
close(GATE)

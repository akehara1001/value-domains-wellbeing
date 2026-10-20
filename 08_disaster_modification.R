# =============================================================================
# 08_disaster_modification.R -- does the domain-count gradient differ between
# respondents who reported disaster damage and those who did not?
#
# The symmetric partner of the family-cancer multiple-group model in
# 07_reference_one_domain.R. Same measurement model, same adjustment set, same
# reference category (one endorsed domain), same df 6 equality test between
# groups.
#
# TWO GROUPS, NOT THREE. Damage judged not repaired is reported by 95 people.
# Split three ways, cells at a level of the domain count fall into single
# digits, so exposure is dichotomised and recovery status is carried as a
# covariate inside the exposed group, where it varies.
#
# WHY THE GROUPS ARE FITTED SEPARATELY. Recovery status has no variance among
# the unexposed, and lavaan refuses any observed variable that is constant in
# a group, whether or not it is treated as fixed. Nothing is constrained across
# groups in the step 6 model either, so the two groups are independent and the
# multiple-group fit factorises into two single-group fits. The check below
# refits step 6 this way and asserts that the coefficients come out identical
# to the published ones and the equality test to three decimals.
#
#   input   values_analysis.rds   (00_derive.R)
#           lca_models.rds        (01_lca.R)
#           outputs/table22_breadth_ref1.csv  (07, for the reproduction check)
#           outputs/table23_ref1_wald.csv     (07, for the reproduction check)
#   output  outputs/table31_disaster_group.csv    coefficients, both groups
#           outputs/table31b_disaster_precise.csv the same at full precision
#           outputs/table32_disaster_wald.csv     equality and joint tests
#           outputs/table33_fourgroup.csv         family cancer x disaster
#           outputs/table34_fourgroup_wald.csv    its test
#           outputs/table35_cells.csv             the four cell counts
# =============================================================================
suppressPackageStartupMessages({ library(lavaan) })
options(stringsAsFactors = FALSE)
BASE <- "."  # the repository root; run every script from there
setwd(BASE)
source("round_half_up.R")  # round() and signif() round half up from here on (see that file)

d  <- readRDS("values_analysis.rds")
LC <- readRDS("lca_models.rds")
K  <- 4L
m  <- LC$models[[K]]

# ------------------------------------------------------------------- data
# Built exactly as in 07_reference_one_domain.R so that every quantity this
# script shares with that one comes out identical.
data <- as.data.frame(d)
data$vclass   <- factor(m$predclass, levels = seq_len(K))
data$age_z    <- as.numeric(scale(data$age))
data$income_z <- as.numeric(scale(data$income))
data$k6_z     <- as.numeric(scale(data$k6))
data$mcom_z   <- as.numeric(scale(data$m_commitment))
FCTS <- c("vclass", "disaster", "fcancer2", "female", "edu2", "occupation2",
          "marriage", "comobidity", "fukushima", "health2")
for (v in FCTS) data[[v]] <- factor(data[[v]])
data <- cbind(data, model.matrix(~ vclass + disaster + fcancer2 + female + edu2 +
                                   occupation2 + marriage + comobidity +
                                   fukushima + health2, data)[, -1])
stopifnot(sum(data$disaster1) == 95, sum(data$disaster2) == 372)

LVL <- c("0", "1", "2", "3", "4", "5", "6+")
data$nsupcat <- cut(pmin(data$nsup, 6), breaks = -1:6, labels = LVL)
for (lv in 0:6) data[[paste0("nsc", lv)]] <- as.numeric(pmin(data$nsup, 6) == lv)

# disaster exposure: any reported damage from a major disaster in the previous
# ten years. Recovery status inside that group: 1 = damage judged not repaired.
data$dexp2 <- factor(as.integer(data$disaster != "0"), levels = c(0, 1))
data$unrec <- as.numeric(data$disaster == "1")
stopifnot(sum(data$dexp2 == 1) == 467, sum(data$dexp2 == 0) == 1712,
          sum(data$unrec[data$dexp2 == 1]) == 95,
          sum(data$unrec[data$dexp2 == 0]) == 0)

MEAS <- paste("NFE =~ sadR + worryR",
              "ResAll =~ ls + gh + es + vital + comp",
              "WB =~ NFE + ResAll", "ls ~~ gh", "es ~~ vital",
              "ResAll ~~ 0*ResAll", sep = "\n")
REF1 <- paste0("nsc", c(0, 2:6))
sfx  <- sub("nsc", "", REF1)

COV_ALL <- c("female1", "age_z", "income_z", "edu22", "edu23", "occupation22",
             "occupation23", "marriage2", "marriage3", "comobidity1",
             "fukushima1", "health22")

model_for <- function(extra) {
  paste0(MEAS, "\nWB ~ ",
         paste(sprintf("a%s*%s", sfx, REF1), collapse = " + "), " + ",
         paste(c(COV_ALL, extra), collapse = " + "), "\n")
}

fitg <- function(dat, extra)
  sem(model_for(extra), dat, estimator = "MLR", missing = "ML")

grab <- function(fit) {
  pe <- parameterEstimates(fit, standardized = TRUE)
  pe <- pe[pe$op == "~" & pe$lhs == "WB" & pe$rhs %in% REF1, ]
  pe <- pe[order(as.integer(sub("nsc", "", pe$rhs))), ]
  idx <- match(pe$label, colnames(vcov(fit)))
  list(pe = pe, b = pe$est, V = vcov(fit)[idx, idx, drop = FALSE])
}

# difference of two independent coefficient vectors, tested jointly
wald2 <- function(g1, g2) {
  dif <- g2$b - g1$b
  V   <- g1$V + g2$V
  s   <- as.numeric(t(dif) %*% solve(V) %*% dif)
  c(s, length(dif), pchisq(s, length(dif), lower.tail = FALSE))
}
waldK <- function(gs) {                    # all groups equal to the first
  dif <- unlist(lapply(gs[-1], function(g) g$b - gs[[1]]$b))
  k   <- length(gs) - 1L
  p   <- length(gs[[1]]$b)
  V   <- matrix(0, k * p, k * p)
  for (i in seq_len(k)) for (j in seq_len(k)) {
    blk <- gs[[1]]$V + if (i == j) gs[[i + 1]]$V else 0
    V[((i - 1) * p + 1):(i * p), ((j - 1) * p + 1):(j * p)] <- blk
  }
  s <- as.numeric(t(dif) %*% solve(V) %*% dif)
  c(s, k * p, pchisq(s, k * p, lower.tail = FALSE))
}

# ------------------------------------------------- reproduction check, step 6
# The two disaster dummies are covariates here, as they are in step 6; family
# cancer history defines the groups and so is dropped from the adjustment set.
fc0 <- grab(fitg(data[data$fcancer2 == "0", ], c("disaster1", "disaster2")))
fc1 <- grab(fitg(data[data$fcancer2 == "1", ], c("disaster1", "disaster2")))
pub  <- read.csv("outputs/table22_breadth_ref1.csv")
pubw <- read.csv("outputs/table23_ref1_wald.csv")
p0 <- pub$b[pub$model == "No family cancer history (ref = 1 domain)"]
p1 <- pub$b[pub$model == "Family cancer history (ref = 1 domain)"]
wfc <- wald2(fc0, fc1)
pw  <- pubw[pubw$test ==
              "Multiple-group model, breadth coefficients equal across groups", ]
stopifnot(identical(round(fc0$b, 4), p0), identical(round(fc1$b, 4), p1),
          identical(round(wfc[1], 3), pw$statistic), wfc[2] == pw$df)
cat("check: step 6 reproduced from two separate fits --",
    "coefficients identical to 4 dp, Wald", round(wfc[1], 3), "df", wfc[2],
    "(published", pw$statistic, "df", pw$df, ")\n")

# ---------------------------------------- step 7: modification by disaster
# family cancer history stays a covariate; recovery status enters only where it
# varies, that is among the exposed
de0 <- grab(fitg(data[data$dexp2 == "0", ], "fcancer21"))
de1 <- grab(fitg(data[data$dexp2 == "1", ], c("fcancer21", "unrec")))
nscat_n_dexp <- table(data$nsupcat, data$dexp2)

pack <- function(g, label, n_src) {
  r <- g$pe
  lv <- LVL[as.integer(sub("nsc", "", r$rhs)) + 1]
  data.frame(model = label, level = lv, n = as.integer(n_src[lv]),
             b = round(r$est, 4), se = round(r$se, 4),
             ci_low = round(r$ci.lower, 4), ci_high = round(r$ci.upper, 4),
             beta = round(r$std.all, 3), p = signif(r$pvalue, 4))
}
pack_raw <- function(g, label) {
  r <- g$pe
  data.frame(model = label,
             level = LVL[as.integer(sub("nsc", "", r$rhs)) + 1],
             b = r$est, ci_low = r$ci.lower, ci_high = r$ci.upper, p = r$pvalue)
}

dex <- rbind(pack(de0, "No disaster damage (ref = 1 domain)",
                  nscat_n_dexp[, 1]),
             pack(de1, "Disaster damage (ref = 1 domain)", nscat_n_dexp[, 2]))
write.csv(dex, "outputs/table31_disaster_group.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(rbind(pack_raw(de0, "No disaster damage (ref = 1 domain)"),
                pack_raw(de1, "Disaster damage (ref = 1 domain)")),
          "outputs/table31b_disaster_precise.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

jointw <- function(g) {
  s <- as.numeric(t(g$b) %*% solve(g$V) %*% g$b)
  c(s, length(g$b), pchisq(s, length(g$b), lower.tail = FALSE))
}
grp_d  <- wald2(de0, de1)
join_0 <- jointw(de0)
join_1 <- jointw(de1)
tests_d <- data.frame(
  test = c("Breadth coefficients equal across disaster exposure groups",
           "No disaster damage, all levels = 0",
           "Disaster damage, all levels = 0"),
  statistic = round(c(grp_d[1], join_0[1], join_1[1]), 3),
  df = c(grp_d[2], join_0[2], join_1[2]),
  p = signif(c(grp_d[3], join_0[3], join_1[3]), 4))
write.csv(tests_d, "outputs/table32_disaster_wald.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
# the manuscript prints Wald statistics to two decimals, so the rounding is
# done once from the full-precision value rather than twice from the table
write.csv(data.frame(test = tests_d$test,
                     statistic = c(grp_d[1], join_0[1], join_1[1]),
                     df = tests_d$df,
                     p = c(grp_d[3], join_0[3], join_1[3])),
          "outputs/table32b_disaster_wald_precise.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
print(tests_d)

# --------------------------------------- the two modifiers taken together
cells <- as.data.frame(table(family_cancer = data$fcancer2,
                             disaster_damage = data$dexp2))
names(cells)[3] <- "n"
write.csv(cells, "outputs/table35_cells.csv", row.names = FALSE,
          fileEncoding = "UTF-8")
print(cells)
MIN_CELL  <- min(cells$n)
THRESHOLD <- 100L        # fixed in the instructions before the cells were seen
cat("smallest of the four cells:", MIN_CELL, "| threshold:", THRESHOLD,
    "| four-group model run:", MIN_CELL >= THRESHOLD, "\n")

if (MIN_CELL >= THRESHOLD) {
  key <- paste0(as.character(data$fcancer2), as.character(data$dexp2))
  LAB <- c("00" = "No family cancer history, no disaster damage",
           "01" = "No family cancer history, disaster damage",
           "10" = "Family cancer history, no disaster damage",
           "11" = "Family cancer history, disaster damage")
  data$cell4 <- factor(LAB[key], levels = LAB)
  stopifnot(!any(is.na(data$cell4)))
  nscat_n_cell <- table(data$nsupcat, data$cell4)
  gs <- lapply(names(LAB), function(k) {
    sub <- data[key == k, ]
    grab(fitg(sub, if (substr(k, 2, 2) == "1") "unrec" else NULL))
  })
  four <- do.call(rbind, lapply(seq_along(gs), function(i)
    pack(gs[[i]], unname(LAB[i]), nscat_n_cell[, i])))
  write.csv(four, "outputs/table33_fourgroup.csv", row.names = FALSE,
            fileEncoding = "UTF-8")
  eq4 <- waldK(gs)
  tests_4 <- data.frame(
    test = "Breadth coefficients equal across the four cells",
    statistic = round(eq4[1], 3), df = eq4[2], p = signif(eq4[3], 4))
  write.csv(tests_4, "outputs/table34_fourgroup_wald.csv", row.names = FALSE,
            fileEncoding = "UTF-8")
  write.csv(data.frame(test = tests_4$test, statistic = eq4[1], df = eq4[2],
                       p = eq4[3]),
            "outputs/table34b_fourgroup_wald_precise.csv", row.names = FALSE,
            fileEncoding = "UTF-8")
  print(tests_4)
} else {
  cat("four-group model not run: smallest cell below the threshold\n")
}

cat("done -- outputs/table31 to table35\n")

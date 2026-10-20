# =============================================================================
# 03_value_change.R -- retrospective value change around a family cancer
#                      diagnosis (Phase 3, secondary and DESCRIPTIVE)
# =============================================================================
#   Rscript 03_value_change.R <K>
#
#   input   values_analysis.rds, lca_models.rds
#   output  outputs/results_phase3.md
#           outputs/phase3_by_domain.csv
#           outputs/phase3_by_class.csv
#           outputs/phase3_by_famcare.csv   (only when sf is present)
#
# The FAMCARE-10 score sf exists only in the dataset built by 00_derive.R from
# the raw export. The released data do not carry it, so with 00_load_release.R
# the FAMCARE tabulation is skipped; no article table uses it.
#
# Sample: the 855 respondents with a family cancer history, and within them only
# the domains each respondent endorsed (Q11_7..Q11_9 are routed on both).
#
# Q11_7 importance / Q11_8 commitment / Q11_9 intention are SINGLE-ITEM
# RETROSPECTIVE change ratings centred on 3 = unchanged. They ask the respondent
# in 2019 to recall a state before the diagnosis, so they are contaminated by
# the current state and cannot support a causal reading. Nothing here is
# modelled: proportions and cross-tabulations only, per the instruction.
# =============================================================================

suppressPackageStartupMessages({library(dplyr)})
options(stringsAsFactors = FALSE)

BASE <- "."  # the repository root; run every script from there
setwd(BASE)
source("round_half_up.R")  # round() and signif() round half up from here on (see that file)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || is.na(suppressWarnings(as.integer(args[1]))))
  stop("usage: Rscript 03_value_change.R <K>")
K <- as.integer(args[1])

t_start <- Sys.time()
d  <- readRDS("values_analysis.rds")
LC <- readRDS("lca_models.rds")
m  <- LC$models[[K]]
DOM <- LC$domains; DOM_LAB <- LC$dom_lab

d$vclass <- factor(m$predclass, levels = seq_len(K))

# FAMCARE tertiles, cut exactly as in paper 1 (quantile breaks on the observed
# sf, include.lowest). Paper 1 reports 297 / 413 / 145; the bands are uneven
# because sf is a discrete sum.
HAS_SF <- "sf" %in% names(d)
if (HAS_SF) {
  d$sf_cat <- cut(d$sf, breaks = quantile(d$sf, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                  include.lowest = TRUE, labels = c("low", "mid", "high"))
  tt <- table(d$sf_cat)
  stopifnot(as.integer(tt) == c(297, 413, 145))
} else {
  d$sf_cat <- factor(NA, levels = c("low", "mid", "high"))
}

# -----------------------------------------------------------------------------
# Long format: one row per (respondent, endorsed domain) among family-cancer
# respondents. Each row carries the three change ratings.
# -----------------------------------------------------------------------------
ITEMS <- c(Q11_7 = "Importance", Q11_8 = "Commitment", Q11_9 = "Intention")
long <- do.call(rbind, lapply(seq_along(DOM), function(j) {
  data.frame(row = d$row, vclass = d$vclass, sf_cat = d$sf_cat,
             fcancer2 = d$fcancer2, endorsed = d[[paste0("v_", DOM[j])]],
             domain = DOM[j], domain_lab = DOM_LAB[j],
             Q11_7 = d[[paste0("Q11_7_", DOM[j])]],
             Q11_8 = d[[paste0("Q11_8_", DOM[j])]],
             Q11_9 = d[[paste0("Q11_9_", DOM[j])]])
}))
long <- long[long$fcancer2 == 1 & long$endorsed == 1, ]
stopifnot(!anyNA(long$Q11_7), !anyNA(long$Q11_8), !anyNA(long$Q11_9))
n_by_dom <- as.integer(table(factor(long$domain, levels = DOM)))
stopifnot(n_by_dom == c(677, 484, 185, 322, 185, 223, 37, 112, 456))

# 1-2 = decreased, 3 = unchanged, 4-5 = increased
band <- function(x) factor(ifelse(x <= 2, "decreased",
                           ifelse(x == 3, "unchanged", "increased")),
                           levels = c("decreased", "unchanged", "increased"))
for (it in names(ITEMS)) long[[paste0("b_", it)]] <- band(long[[it]])

pct3 <- function(b) {
  n <- length(b); t <- table(b)
  c(n = n, dec = 100 * t[["decreased"]] / n, unc = 100 * t[["unchanged"]] / n,
    inc = 100 * t[["increased"]] / n)
}

# -----------------------------------------------------------------------------
# (1) by domain
# -----------------------------------------------------------------------------
by_dom <- do.call(rbind, lapply(names(ITEMS), function(it)
  do.call(rbind, lapply(DOM, function(dm) {
    z <- long[long$domain == dm, paste0("b_", it)]
    data.frame(item = ITEMS[[it]], domain = dm,
               domain_lab = DOM_LAB[match(dm, DOM)], t(pct3(z)))
  }))))
by_dom[, 4:7] <- lapply(by_dom[, 4:7], function(x) round(x, 1))
by_dom$n <- as.integer(by_dom$n)
write.csv(by_dom, "outputs/phase3_by_domain.csv", row.names = FALSE, fileEncoding = "UTF-8")

# -----------------------------------------------------------------------------
# (2) by value class, (3) by FAMCARE tertile
# -----------------------------------------------------------------------------
# Two summaries are given for each stratum:
#   response level  every (person, endorsed domain) rating counts once. People
#                   who endorsed more domains contribute more rows, and the rows
#                   are clustered within people. Descriptive only.
#   person level    each person's share of their OWN endorsed domains rated
#                   "increased", averaged over people in the stratum. One value
#                   per person, so no clustering.
strat_tab <- function(stratvar, levs) {
  do.call(rbind, lapply(names(ITEMS), function(it) {
    do.call(rbind, lapply(levs, function(lv) {
      sel <- long[!is.na(long[[stratvar]]) & long[[stratvar]] == lv, ]
      if (!nrow(sel)) return(NULL)
      b <- sel[[paste0("b_", it)]]
      per <- tapply(b == "increased", sel$row, mean)
      data.frame(item = ITEMS[[it]], stratum = as.character(lv),
                 n_persons = length(unique(sel$row)), n_responses = nrow(sel),
                 dec_pct = round(100 * mean(b == "decreased"), 1),
                 unc_pct = round(100 * mean(b == "unchanged"), 1),
                 inc_pct = round(100 * mean(b == "increased"), 1),
                 person_mean_inc_pct = round(100 * mean(per), 1))
    }))
  }))
}
by_class   <- strat_tab("vclass", levels(long$vclass))
write.csv(by_class,   "outputs/phase3_by_class.csv",   row.names = FALSE, fileEncoding = "UTF-8")
if (HAS_SF) {
  by_famcare <- strat_tab("sf_cat", c("low", "mid", "high"))
  write.csv(by_famcare, "outputs/phase3_by_famcare.csv", row.names = FALSE, fileEncoding = "UTF-8")
}

# domain x class and domain x tertile detail
detail <- do.call(rbind, lapply(names(ITEMS), function(it)
  do.call(rbind, lapply(DOM, function(dm) {
    do.call(rbind, lapply(levels(long$vclass), function(k) {
      sel <- long[long$domain == dm & long$vclass == k, ]
      if (!nrow(sel)) return(NULL)
      b <- sel[[paste0("b_", it)]]
      data.frame(item = ITEMS[[it]], domain = dm, class = k, n = nrow(sel),
                 dec_pct = round(100 * mean(b == "decreased"), 1),
                 unc_pct = round(100 * mean(b == "unchanged"), 1),
                 inc_pct = round(100 * mean(b == "increased"), 1))
    }))
  }))))
write.csv(detail, "outputs/phase3_domain_by_class.csv", row.names = FALSE, fileEncoding = "UTF-8")

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
con <- file("outputs/results_phase3.md", open = "wt", encoding = "UTF-8")
w <- function(...) cat(sprintf(...), file = con)

w("# Phase 3 \u7d50\u679c \u2014 \u304c\u3093\u7f79\u60a3\u5bb6\u65cf\u306b\u304a\u3051\u308b\u4fa1\u5024\u306e\u5909\u5316\uff08\u526f\u6b21\u30fb\u8a18\u8ff0\uff09\n\n")
w("\u751f\u6210: `03_value_change.R %d`\uff08%s\uff09\n\n", K, format(t_start, "%Y-%m-%d %H:%M"))
w("**\u56e0\u679c\u7684\u306b\u8aad\u307e\u306a\u3044\u3053\u3068\u3002** Q11_7\u301c9 \u306f 3 = \u5909\u308f\u3089\u306a\u3044 \u3092\u4e2d\u5fc3\u3068\u3059\u308b\u5358\u4e00\u9805\u76ee\u306e\u56de\u9867\u7684\u5909\u5316\u8a55\u4fa1\u3067\u3042\u308a\u3001")
w("2019\u5e74\u6642\u70b9\u306e\u56de\u7b54\u8005\u304c\u8a3a\u65ad\u524d\u306e\u72b6\u614b\u3092\u601d\u3044\u51fa\u3057\u3066\u7b54\u3048\u3066\u3044\u308b\u3002\u73fe\u5728\u306e\u72b6\u614b\u306b\u5f15\u304d\u305a\u3089\u308c\u308b\u3002")
w("\u56de\u5e30\u306f\u7d44\u3093\u3067\u3044\u306a\u3044\u3002\n\n")
w("\u5bfe\u8c61: \u5bb6\u65cf\u306e\u304c\u3093\u7d4c\u9a13\u3042\u308a %d \u540d\u306e\u3046\u3061\u3001\u5f53\u8a72\u9818\u57df\u3092\u652f\u6301\u3057\u305f\u8005\u3002",
  sum(d$fcancer2 == 1))
w("\u5ef6\u3079\u56de\u7b54\u6570 %d\uff08\u4eba\u00d7\u652f\u6301\u9818\u57df\uff09\u3002\n\n", nrow(long))
w("1-2 = \u4f4e\u4e0b\u30013 = \u4e0d\u5909\u30014-5 = \u4e0a\u6607\u3002\n\n")

w("## 1. \u9818\u57df\u5225\n\n")
for (it in names(ITEMS)) {
  w("### %s (%s)\n\n", ITEMS[[it]], it)
  w("| \u9818\u57df | n | \u4f4e\u4e0b %% | \u4e0d\u5909 %% | \u4e0a\u6607 %% |\n|---|---|---|---|---|\n")
  z <- by_dom[by_dom$item == ITEMS[[it]], ]
  z <- z[match(DOM, z$domain), ]
  for (i in seq_len(nrow(z)))
    w("| %s | %d | %.1f | %.1f | %.1f |\n", z$domain_lab[i], z$n[i],
      z$dec[i], z$unc[i], z$inc[i])
  w("\n")
}

w("## 2. \u4fa1\u5024\u30af\u30e9\u30b9\u5225\n\n")
w("\u5ef6\u3079\u56de\u7b54\u30ec\u30d9\u30eb\u306e\u5272\u5408\u306f\u3001\u652f\u6301\u9818\u57df\u306e\u591a\u3044\u4eba\u307b\u3069\u91cd\u304f\u306a\u308a\u3001\u540c\u4e00\u4eba\u5185\u3067\u76f8\u95a2\u3059\u308b\u3002")
w("\u305d\u306e\u305f\u3081\u3001\u4eba\u5358\u4f4d\uff08\u81ea\u5206\u304c\u652f\u6301\u3057\u305f\u9818\u57df\u306e\u3046\u3061\u4e0a\u6607\u3068\u7b54\u3048\u305f\u5272\u5408\u306e\u500b\u4eba\u5e73\u5747\uff09\u3082\u4e26\u8a18\u3059\u308b\u3002\n\n")
for (it in names(ITEMS)) {
  w("### %s\n\n", ITEMS[[it]])
  w("| \u30af\u30e9\u30b9 | \u4eba\u6570 | \u5ef6\u3079\u56de\u7b54 | \u4f4e\u4e0b %% | \u4e0d\u5909 %% | \u4e0a\u6607 %% | \u4eba\u5358\u4f4d\u4e0a\u6607 %% |\n")
  w("|---|---|---|---|---|---|---|\n")
  z <- by_class[by_class$item == ITEMS[[it]], ]
  for (i in seq_len(nrow(z)))
    w("| %s | %d | %d | %.1f | %.1f | %.1f | %.1f |\n", z$stratum[i],
      z$n_persons[i], z$n_responses[i], z$dec_pct[i], z$unc_pct[i],
      z$inc_pct[i], z$person_mean_inc_pct[i])
  w("\n")
}

if (HAS_SF) {
w("## 3. FAMCARE \u4e09\u5206\u4f4d\u5225\n\n")
w("\u4e09\u5206\u4f4d\u306e\u5207\u65ad\u306f\u7b2c1\u8ad6\u6587\u3068\u540c\u4e00\uff08low %d / mid %d / high %d\uff09\u3002sf \u304c\u96e2\u6563\u548c\u306e\u305f\u3081\u5e45\u306f\u4e0d\u5747\u7b49\u3002\n\n",
  tt[["low"]], tt[["mid"]], tt[["high"]])
for (it in names(ITEMS)) {
  w("### %s\n\n", ITEMS[[it]])
  w("| \u4e09\u5206\u4f4d | \u4eba\u6570 | \u5ef6\u3079\u56de\u7b54 | \u4f4e\u4e0b %% | \u4e0d\u5909 %% | \u4e0a\u6607 %% | \u4eba\u5358\u4f4d\u4e0a\u6607 %% |\n")
  w("|---|---|---|---|---|---|---|\n")
  z <- by_famcare[by_famcare$item == ITEMS[[it]], ]
  for (i in seq_len(nrow(z)))
    w("| %s | %d | %d | %.1f | %.1f | %.1f | %.1f |\n", z$stratum[i],
      z$n_persons[i], z$n_responses[i], z$dec_pct[i], z$unc_pct[i],
      z$inc_pct[i], z$person_mean_inc_pct[i])
  w("\n")
}
}

w("## \u51fa\u5178\n\n")
w("\u5168\u6570\u5024\u306f `03_value_change.R` \u304c\u751f\u6210\u3002\u9818\u57df\u00d7\u30af\u30e9\u30b9\u306e\u7d30\u90e8\u306f ")
w("outputs/phase3_domain_by_class.csv\u3002\u9818\u57df\u5225\u306e n \u304c\u5c0f\u3055\u3044\u30bb\u30eb")
w("\uff08spiritual \u306f\u5168\u4f53\u3067 %d\uff09\u3092\u5272\u5408\u3060\u3051\u3067\u8aad\u307e\u306a\u3044\u3053\u3068\u3002\n",
  n_by_dom[7])
close(con)

cat(sprintf("done in %.1f s -- outputs/results_phase3.md\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

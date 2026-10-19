# =============================================================================
# 01_lca.R -- latent class enumeration on the nine value-domain indicators
# =============================================================================
# PHASE 1. This script ENUMERATES 1..6 class solutions and STOPS. It does not
# choose a number of classes and nothing downstream depends on it choosing one:
# the class count is the authors' decision, taken from outputs/lca_enumeration.md.
#
#   input   values_analysis.rds  (built and gated by 00_derive.R)
#   output  lca_models.rds             all six fitted poLCA objects, reordered
#           outputs/lca_enumeration.md the enumeration report
#           outputs/lca_fit_table.csv  the fit table as data
#           outputs/lca_profiles.csv   conditional probabilities, all solutions
#           outputs/fig_profiles_K*.png  one profile plot per solution
#
# Indicators are the nine Q11_1 endorsements, all observed for all 2179
# respondents. The 215 who endorsed nothing are retained as an all-zero row.
# Q11_3..Q11_6 are routed items and are NOT indicators.
#
# Run:  Rscript 01_lca.R
# =============================================================================

suppressPackageStartupMessages({library(poLCA)})
options(stringsAsFactors = FALSE)

BASE <- "."  # the repository root; run every script from there
setwd(BASE)
dir.create("outputs", showWarnings = FALSE)

SEED  <- 20260829L   # fixed and reported; changing it changes nothing below
NREP  <- 100L        # random starts per solution
MAXIT <- 50000L
KMAX  <- 6L

t_start <- Sys.time()
d <- readRDS("values_analysis.rds")
stopifnot(nrow(d) == 2179)

DOM <- c("family", "friends", "romance", "work", "education",
         "leisure", "spiritual", "community", "health")
DOM_LAB <- c("Family", "Friends", "Romance", "Work", "Education",
             "Leisure", "Spiritual", "Community", "Health")
IND <- paste0("v_", DOM)

# poLCA requires manifest variables coded as consecutive positive integers, so
# the 0/1 endorsements become 1 = not endorsed, 2 = endorsed. Every reference to
# "level 2" below is therefore the endorsement probability.
lc <- d[, IND]
stopifnot(!anyNA(lc), all(unlist(lc) %in% c(0, 1)))
lc <- as.data.frame(lapply(lc, function(x) x + 1L))
f <- as.formula(paste0("cbind(", paste(IND, collapse = ", "), ") ~ 1"))

n <- nrow(lc)
marg <- colMeans(d[, IND])           # overall endorsement rates
mean_nsup_all <- mean(d$nsup)

# -----------------------------------------------------------------------------
# Fit 1..KMAX
# -----------------------------------------------------------------------------
entropy_R2 <- function(post) {
  K <- ncol(post)
  if (K < 2) return(NA_real_)
  p <- post
  p[p <= 0] <- 1e-300
  1 + sum(post * log(p)) / (n * log(K))
}

# Estimation is expensive (K=5 needs ~7000 EM iterations per start, so the full
# enumeration is over an hour). The fitted models are therefore CACHED and are
# written to disk the moment fitting finishes, before any reporting code runs:
# a formatting error further down must never cost the estimation again.
CACHE <- "lca_models_raw.rds"
cache_key <- list(seed = SEED, nrep = NREP, maxiter = MAXIT, kmax = KMAX,
                  ind = IND, n = nrow(lc))
cached <- if (file.exists(CACHE)) readRDS(CACHE) else NULL
if (!is.null(cached) && identical(cached$key, cache_key)) {
  cat("reusing cached fits from", CACHE, "\n")
  models <- cached$models
} else {
  set.seed(SEED)
  models <- vector("list", KMAX)
  cat("fitting:")
  for (K in 1:KMAX) {
    cat(" K=", K, sep = "")
    models[[K]] <- poLCA(f, lc, nclass = K, nrep = if (K == 1) 1L else NREP,
                         maxiter = MAXIT, tol = 1e-10, na.rm = FALSE,
                         graphs = FALSE, verbose = FALSE)
  }
  cat("\n")
  saveRDS(list(key = cache_key, models = models), CACHE)
}

# Reorder classes by descending size, so class 1 is always the largest. poLCA's
# own class order is an artefact of the winning random start and is not stable
# across K; every table and figure below uses the size-ordered labelling.
reorder_classes <- function(m) {
  if (length(m$P) == 1) return(m)
  o <- order(m$P, decreasing = TRUE)
  m$P <- m$P[o]
  m$probs <- lapply(m$probs, function(pr) {
    pr <- pr[o, , drop = FALSE]; rownames(pr) <- paste0("class ", seq_len(nrow(pr))); pr
  })
  m$posterior <- m$posterior[, o, drop = FALSE]
  m$predclass <- apply(m$posterior, 1, which.max)
  m
}
models <- lapply(models, reorder_classes)

# -----------------------------------------------------------------------------
# Fit table
# -----------------------------------------------------------------------------
fit_rows <- lapply(1:KMAX, function(K) {
  m <- models[[K]]
  post <- m$posterior
  if (is.null(dim(post))) post <- matrix(1, nrow = n, ncol = 1)
  cls <- m$predclass
  sizes <- as.integer(table(factor(cls, levels = seq_along(m$P))))
  att <- m$attempts
  best_reps <- if (is.null(att)) NA_integer_ else sum(abs(att - max(att)) < 0.01)
  # sample-size adjusted BIC uses the Sclove n* = (n + 2) / 24
  abic <- -2 * m$llik + m$npar * log((n + 2) / 24)
  data.frame(
    K = K, logLik = m$llik, npar = m$npar,
    AIC = m$aic, BIC = m$bic, aBIC = abic,
    entropy = entropy_R2(post),
    min_class_n = min(sizes), min_class_pct = 100 * min(sizes) / n,
    converged = isTRUE(m$numiter < MAXIT),
    best_of_nrep = best_reps,
    nrep = if (K == 1) 1L else NREP)
})
fit <- do.call(rbind, fit_rows)
write.csv(fit, "outputs/lca_fit_table.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Profiles: conditional endorsement probability by class, plus the count check
# -----------------------------------------------------------------------------
# The dominant axis in these data is simply how many domains a person endorses
# (0 for 215 people, 9 for 10). A class solution that only re-expressed that
# count would be of no interest, so for every solution two things are recorded:
#
#   share      the fraction of a class's expected endorsements falling on each
#              domain, P(j|k) / sum_j P(j|k), against the same share overall.
#              Purely descriptive, no model.
#   tilt resid the residual from a one-parameter-per-class benchmark in which
#              classes differ ONLY in overall propensity:
#                  logit P(j|k) = beta_j + alpha_k
#              beta_j is the overall domain logit, and alpha_k is solved so the
#              benchmark reproduces that class's expected number of endorsements
#              exactly. Whatever remains is composition that a count gradient
#              cannot produce.
#
# Neither quantity is an inferential statistic; both are descriptions used to
# answer "is this class structure just a count gradient?".
logit  <- function(p) log(p / (1 - p))
expit  <- function(x) 1 / (1 + exp(-x))
beta_j <- logit(pmin(pmax(marg, 1e-6), 1 - 1e-6))

tilt_alpha <- function(target_count) {
  # alpha such that sum_j expit(beta_j + alpha) == target_count
  fn <- function(a) sum(expit(beta_j + a)) - target_count
  if (target_count <= 1e-8) return(-Inf)
  if (target_count >= 9 - 1e-8) return(Inf)
  uniroot(fn, c(-30, 30), tol = 1e-10)$root
}

prof_list <- list()
for (K in 1:KMAX) {
  m <- models[[K]]
  P <- sapply(IND, function(v) m$probs[[v]][, 2])          # K x 9 (or vector)
  if (K == 1) P <- matrix(P, nrow = 1, dimnames = list(NULL, IND))
  cls <- m$predclass
  for (k in seq_len(K)) {
    pk <- P[k, ]
    exp_count <- sum(pk)
    a <- tilt_alpha(exp_count)
    bench <- if (is.finite(a)) expit(beta_j + a) else rep(exp_count / 9, 9)
    prof_list[[length(prof_list) + 1]] <- data.frame(
      K = K, class = k,
      class_share = m$P[k],
      n_modal = sum(cls == k),
      domain = DOM,
      p_endorse = as.numeric(pk),
      exp_n_domains = exp_count,
      mean_nsup_modal = mean(d$nsup[cls == k]),
      share = as.numeric(pk / exp_count),
      share_overall = as.numeric(marg / mean_nsup_all),
      tilt_benchmark = as.numeric(bench),
      tilt_resid = as.numeric(pk - bench))
  }
}
prof <- do.call(rbind, prof_list)
write.csv(prof, "outputs/lca_profiles.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Figures: one profile plot per solution
# -----------------------------------------------------------------------------
pal <- c("#1f4e79", "#c00000", "#2e7d32", "#e69f00", "#7b3294", "#0097a7")
for (K in 2:KMAX) {
  m <- models[[K]]
  P <- sapply(IND, function(v) m$probs[[v]][, 2])
  png(sprintf("outputs/fig_profiles_K%d.png", K), width = 1500, height = 950, res = 150)
  par(mar = c(7.5, 4.5, 3, 1))
  plot(NA, xlim = c(1, 9), ylim = c(0, 1), xaxt = "n", xlab = "",
       ylab = "Conditional probability of endorsement",
       main = sprintf("%d-class solution", K))
  axis(1, at = 1:9, labels = DOM_LAB, las = 2)
  abline(h = seq(0, 1, 0.2), col = "grey90")
  lines(1:9, marg, col = "grey45", lwd = 2, lty = 3, type = "b", pch = 4)
  for (k in 1:K) lines(1:9, P[k, ], col = pal[k], lwd = 2.2, type = "b", pch = 19)
  legend("topright", bty = "n", cex = 0.8,
         legend = c(sprintf("Class %d (%.1f%%)", 1:K, 100 * m$P), "Overall rate"),
         col = c(pal[1:K], "grey45"), lwd = 2, lty = c(rep(1, K), 3),
         pch = c(rep(19, K), 4))
  dev.off()
}

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
con <- file("outputs/lca_enumeration.md", open = "wt", encoding = "UTF-8")
w <- function(...) cat(sprintf(...), file = con)
nl <- function() cat("\n", file = con)

w("# \u4fa1\u5024\u9818\u57df\u306e\u6f5c\u5728\u30af\u30e9\u30b9\u5206\u6790 \u2014 \u30af\u30e9\u30b9\u6570\u306e\u5217\u6319\n\n")
w("\u751f\u6210: `01_lca.R`\uff08%s\uff09  \u4e71\u6570\u7a2e %d\uff0f\u30e9\u30f3\u30c0\u30e0\u30b9\u30bf\u30fc\u30c8 %d\u56de\uff0fmaxiter %d\n\n",
  format(t_start, "%Y-%m-%d %H:%M"), SEED, NREP, MAXIT)
w("\u6307\u6a19: Q11_1[1]\u301c[9] \u306e\u4e8c\u5024\u652f\u6301\u3002N = %d\u3001\u6b20\u6e2c\u306a\u3057\u3002",  n)
w("\u652f\u6301\u9818\u57df\u306a\u3057\u306e %d \u540d\u306f9\u6307\u6a19\u3059\u3079\u3066 0 \u306e\u884c\u3068\u3057\u3066\u542b\u3081\u3066\u3044\u308b\u3002\n\n",
  sum(d$nsup == 0))
w("**\u3053\u306e\u6bb5\u968e\u3067\u306f\u30af\u30e9\u30b9\u6570\u3092\u6c7a\u3081\u3066\u3044\u306a\u3044\u3002**\u4ee5\u4e0b\u306f\u5224\u65ad\u6750\u6599\u3067\u3042\u308a\u3001\u6c7a\u5b9a\u306f\u8457\u8005\u304c\u884c\u3046\u3002\n\n")

w("## 1. \u9069\u5408\u6307\u6a19\n\n")
w("| K | logLik | npar | AIC | BIC | aBIC | entropy | \u6700\u5c0f\u30af\u30e9\u30b9 n | \u6700\u5c0f\u30af\u30e9\u30b9 %% | \u53ce\u675f | \u6700\u5c24\u89e3\u306e\u518d\u73fe |\n")
w("|---|---|---|---|---|---|---|---|---|---|---|\n")
for (i in 1:KMAX) {
  r <- fit[i, ]
  w("| %d | %.1f | %d | %.1f | %.1f | %.1f | %s | %d | %.1f | %s | %s |\n",
    r$K, r$logLik, r$npar, r$AIC, r$BIC, r$aBIC,
    if (is.na(r$entropy)) "\u2014" else sprintf("%.3f", r$entropy),
    r$min_class_n, r$min_class_pct,
    if (r$converged) "\u3042\u308a" else "**\u306a\u3057**",
    if (r$K == 1) "\u2014" else sprintf("%d / %d", r$best_of_nrep, r$nrep))
}
nl()
bBIC <- which.min(fit$BIC); baBIC <- which.min(fit$aBIC); bAIC <- which.min(fit$AIC)
w("BIC \u6700\u5c0f = %d \u30af\u30e9\u30b9\u3001aBIC \u6700\u5c0f = %d \u30af\u30e9\u30b9\u3001AIC \u6700\u5c0f = %d \u30af\u30e9\u30b9\u3002\n\n",
  bBIC, baBIC, bAIC)

w("## 2. \u5206\u985e\u306e\u78ba\u304b\u3089\u3057\u3055\uff08\u5e73\u5747\u4e8b\u5f8c\u78ba\u7387\u306e\u5bfe\u89d2\u8981\u7d20\uff09\n\n")
w("\u5404\u30af\u30e9\u30b9\u306b\u6700\u5927\u4e8b\u5f8c\u78ba\u7387\u3067\u5272\u308a\u5f53\u3066\u3089\u308c\u305f\u8005\u306e\u3001\u305d\u306e\u30af\u30e9\u30b9\u306e\u4e8b\u5f8c\u78ba\u7387\u5e73\u5747\u3002\n\n")
w("| K | %s |\n", paste(sprintf("\u30af\u30e9\u30b9%d", 1:KMAX), collapse = " | "))
w("|---|%s\n", paste(rep("---|", KMAX), collapse = ""))
for (K in 2:KMAX) {
  m <- models[[K]]; post <- m$posterior; cls <- m$predclass
  app <- sapply(1:K, function(k) mean(post[cls == k, k]))
  cells <- c(sprintf("%.3f", app), rep("", KMAX - K))
  w("| %d | %s |\n", K, paste(cells, collapse = " | "))
}
nl()

w("## 3. \u30af\u30e9\u30b9\u5225\u30d7\u30ed\u30d5\u30a1\u30a4\u30eb\uff08\u6761\u4ef6\u4ed8\u304d\u652f\u6301\u78ba\u7387\uff09\n\n")
w("\u56f3\u306f `outputs/fig_profiles_K2.png` \u301c `fig_profiles_K6.png`\u3002\u5168\u4f53\u306e\u652f\u6301\u7387\u306f ")
w("%s\u3002\n\n", paste(sprintf("%s %.1f%%", DOM_LAB, 100 * marg), collapse = "\u3001"))
for (K in 2:KMAX) {
  w("### K = %d\n\n", K)
  m <- models[[K]]
  pk <- prof[prof$K == K, ]
  # n and "%" must be the SAME quantity. n is the modal-assignment count, so the
  # percentage beside it is the modal share; the model's own class share (the
  # mixing proportion) is a different number and gets its own column.
  w("| \u30af\u30e9\u30b9 | n\uff08\u6700\u5927\u4e8b\u5f8c\u5272\u5f53\uff09| \u5272\u5f53 %% | \u30e2\u30c7\u30eb\u30af\u30e9\u30b9\u6bd4 %% | \u5e73\u5747\u652f\u6301\u6570\uff08\u5272\u5f53\uff09| \u30e2\u30c7\u30eb\u671f\u5f85\u652f\u6301\u6570 | %s |\n",
    paste(DOM_LAB, collapse = " | "))
  w("|---|---|---|---|---|---|%s\n", paste(rep("---|", 9), collapse = ""))
  for (k in 1:K) {
    z <- pk[pk$class == k, ]
    w("| %d | %d | %.1f | %.1f | %.2f | %.2f | %s |\n", k, z$n_modal[1],
      100 * z$n_modal[1] / n, 100 * z$class_share[1],
      z$mean_nsup_modal[1], z$exp_n_domains[1],
      paste(sprintf("%.2f", z$p_endorse), collapse = " | "))
  }
  nl()
}

w("## 4. \u300c\u652f\u6301\u6570\u306e\u591a\u5be1\u3060\u3051\u300d\u3067\u8aac\u660e\u3067\u304d\u308b\u304b\u306e\u70b9\u691c\n\n")
w("\u3053\u306e\u30c7\u30fc\u30bf\u306e\u6700\u5927\u306e\u5206\u6563\u8ef8\u306f\u652f\u6301\u9818\u57df\u6570\u305d\u306e\u3082\u306e\u3067\u3042\u308b\uff080 \u304c %d \u540d\u30019 \u304c %d \u540d\uff09\u3002",
  sum(d$nsup == 0), sum(d$nsup == 9))
w("\u3057\u305f\u304c\u3063\u3066\u3001\u30af\u30e9\u30b9\u304c\u5358\u306b\u652f\u6301\u6570\u306e\u5c64\u3092\u8a00\u3044\u63db\u3048\u305f\u3060\u3051\u3067\u306a\u3044\u3053\u3068\u3092\u78ba\u304b\u3081\u308b\u5fc5\u8981\u304c\u3042\u308b\u3002\n\n")
w("\u57fa\u6e96\u30e2\u30c7\u30eb: `logit P(\u9818\u57dfj | \u30af\u30e9\u30b9k) = beta_j + alpha_k`\u3002")
w("beta_j \u306f\u5168\u4f53\u306e\u9818\u57df\u5225\u30ed\u30b8\u30c3\u30c8\u3001alpha_k \u306f\u305d\u306e\u30af\u30e9\u30b9\u306e\u671f\u5f85\u652f\u6301\u9818\u57df\u6570\u3092")
w("\u518d\u73fe\u3059\u308b\u3088\u3046\u306b\u4e00\u610f\u306b\u5b9a\u3081\u308b\u3002\u3059\u306a\u308f\u3061\u3053\u306e\u57fa\u6e96\u306f\u300c\u30af\u30e9\u30b9\u306f\u652f\u6301\u306e\u3057\u3084\u3059\u3055\u3060\u3051\u304c\u9055\u3046\u300d")
w("\u3068\u3044\u3046\u4eee\u5b9a\u305d\u306e\u3082\u306e\u3067\u3042\u308a\u3001\u6b8b\u5dee\u306f\u652f\u6301\u6570\u306e\u52fe\u914d\u3067\u306f\u4f5c\u308c\u306a\u3044\u9818\u57df\u69cb\u6210\u3092\u8868\u3059\u3002")
w("\u691c\u5b9a\u3067\u306f\u306a\u304f\u8a18\u8ff0\u3067\u3042\u308b\u3002\n\n")
for (K in 2:KMAX) {
  pk <- prof[prof$K == K, ]
  w("### K = %d\n\n", K)
  w("| \u30af\u30e9\u30b9 | \u671f\u5f85\u652f\u6301\u6570 | \u6b8b\u5dee\u6700\u5927\u7d76\u5bfe\u5024 | \u9818\u57df\u5225\u6b8b\u5dee\uff08\u89b3\u6e2c \u2212 \u57fa\u6e96\uff09|\n")
  w("|---|---|---|---|\n")
  for (k in 1:K) {
    z <- pk[pk$class == k, ]
    ord <- order(-abs(z$tilt_resid))
    top <- head(ord, 4)
    w("| %d | %.2f | %.2f | %s |\n", k, z$exp_n_domains[1], max(abs(z$tilt_resid)),
      paste(sprintf("%s %+.2f", DOM_LAB[top], z$tilt_resid[top]), collapse = ", "))
  }
  nl()
  # pairs of classes with similar expected counts but different composition
  ec <- tapply(pk$exp_n_domains, pk$class, function(x) x[1])
  sh <- do.call(rbind, lapply(1:K, function(k) pk$share[pk$class == k]))
  pairs_txt <- character(0)
  for (a in 1:(K - 1)) for (b in (a + 1):K) {
    if (abs(ec[a] - ec[b]) <= 1.0) {
      dmax <- max(abs(sh[a, ] - sh[b, ]))
      j <- which.max(abs(sh[a, ] - sh[b, ]))
      pairs_txt <- c(pairs_txt, sprintf(
        "- \u30af\u30e9\u30b9%d\uff08\u671f\u5f85%.2f\uff09\u3068\u30af\u30e9\u30b9%d\uff08\u671f\u5f85%.2f\uff09: \u652f\u6301\u6570\u306f\u8fd1\u3044\u304c\u3001\u69cb\u6210\u6bd4\u306e\u6700\u5927\u5dee\u306f %s \u3067 %.3f\n",
        a, ec[a], b, ec[b], DOM_LAB[j], dmax))
    }
  }
  if (length(pairs_txt)) {
    w("\u652f\u6301\u9818\u57df\u6570\u304c 1.0 \u4ee5\u5185\u3067\u63c3\u3046\u30af\u30e9\u30b9\u5bfe\u306e\u3046\u3061\u3001\u69cb\u6210\u6bd4\uff08\u305d\u306e\u30af\u30e9\u30b9\u306e\u652f\u6301\u304c\u5404\u9818\u57df\u306b\u5206\u914d\u3055\u308c\u308b\u5272\u5408\uff09\u304c\u9055\u3046\u3082\u306e:\n\n")
    for (s in pairs_txt) cat(s, file = con)
  } else {
    w("\u652f\u6301\u9818\u57df\u6570\u304c 1.0 \u4ee5\u5185\u3067\u63c3\u3046\u30af\u30e9\u30b9\u5bfe\u306f\u306a\u3044\uff08= \u3053\u306e\u89e3\u306e\u30af\u30e9\u30b9\u306f\u652f\u6301\u6570\u3067\u4e26\u3093\u3067\u3044\u308b\uff09\u3002\n")
  }
  nl()
  # class x nsup cross-tabulation
  m <- models[[K]]
  ct <- table(m$predclass, d$nsup)
  w("\u30af\u30e9\u30b9 \u00d7 \u652f\u6301\u9818\u57df\u6570\u306e\u30af\u30ed\u30b9\u96c6\u8a08:\n\n")
  w("| \u30af\u30e9\u30b9 | %s |\n", paste(colnames(ct), collapse = " | "))
  w("|---|%s\n", paste(rep("---|", ncol(ct)), collapse = ""))
  for (k in 1:K) w("| %d | %s |\n", k, paste(ct[k, ], collapse = " | "))
  nl()
}

w("## 5. \u6b21\u306e\u4e00\u624b\n\n")
w("\u30af\u30e9\u30b9\u6570\u304c\u6c7a\u307e\u3063\u305f\u3089 `02_structural.R` \u3092\u8d70\u3089\u305b\u308b\u3002")
w("\u3053\u306e\u30b9\u30af\u30ea\u30d7\u30c8\u306f\u30af\u30e9\u30b9\u6570\u3092\u9078\u3070\u305a\u3001\u5168\u89e3\u3092 `lca_models.rds` \u306b\u4fdd\u5b58\u3057\u3066\u7d42\u308f\u308b\u3002\n")
close(con)

saveRDS(list(models = models, fit = fit, prof = prof, seed = SEED,
             nrep = NREP, indicators = IND, domains = DOM, dom_lab = DOM_LAB),
        "lca_models.rds")

cat(sprintf("\ndone in %.1f s -- see outputs/lca_enumeration.md\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
print(fit[, c("K", "logLik", "npar", "BIC", "aBIC", "entropy",
              "min_class_n", "min_class_pct", "best_of_nrep")], row.names = FALSE)

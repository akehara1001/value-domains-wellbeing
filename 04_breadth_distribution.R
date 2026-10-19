# 04_breadth_distribution.R -- frequency and percentage of breadth of endorsement.
#
# The counts were already reported in outputs/results_phase2.md; the percentages
# were not, and the manuscript now reports both. This script derives the table
# from the analysis dataset so that every figure in the text traces to an
# output file rather than to arithmetic done by hand.
#
#   Rscript 04_breadth_distribution.R

d <- readRDS("values_analysis.rds")
stopifnot(nrow(d) == 2179, !any(is.na(d$nsup)))

lev <- pmin(d$nsup, 6)
tab <- table(factor(lev, levels = 0:6))
pct <- round(100 * as.numeric(tab) / nrow(d), 1)

out <- data.frame(
  level = c("0", "1", "2", "3", "4", "5", "6+"),
  n     = as.integer(tab),
  pct   = pct,
  stringsAsFactors = FALSE)

stopifnot(sum(out$n) == 2179)
write.csv(out, "outputs/table16_breadth_distribution.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("Breadth of endorsement, N =", nrow(d), "\n")
print(out, row.names = FALSE)
cat("sum n =", sum(out$n), "  sum pct =", sum(out$pct), "\n")

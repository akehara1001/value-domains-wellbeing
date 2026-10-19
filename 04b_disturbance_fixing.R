# 04b_disturbance_fixing.R -- effect of fixing the ResAll disturbance at zero.
#
# The manuscript's Limitations state what fixing the second-order disturbance
# cost. results_phase2.md reports the free estimate and the class coefficients
# rounded to three decimals; this script recovers the exact change so the text
# quotes a figure that exists in an output file.
#
#   Rscript 04b_disturbance_fixing.R

suppressMessages(library(lavaan))
f <- readRDS("structural_fits.rds")

key <- function(d) paste(d$lhs, d$op, d$rhs)
fx <- standardizedSolution(f$fitA)
fr <- standardizedSolution(f$fitA_free)
m <- merge(data.frame(k = key(fx), fixed = fx$est.std),
           data.frame(k = key(fr), free = fr$est.std), by = "k")
m$abs_diff <- abs(m$fixed - m$free)

cls <- m[grepl("^WB ~ vclass", m$k), ]
cls <- cls[order(-cls$abs_diff), ]

out <- data.frame(
  term      = sub("^WB ~ ", "", cls$k),
  std_fixed = round(cls$fixed, 6),
  std_free  = round(cls$free, 6),
  abs_diff  = round(cls$abs_diff, 6),
  stringsAsFactors = FALSE)

write.csv(out, "outputs/table17_disturbance_fixing.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("Standardised class coefficients, disturbance fixed at 0 vs freely estimated\n")
print(out, row.names = FALSE)
cat("\nlargest change, 4 dp:", sprintf("%.4f", max(cls$abs_diff)), "\n")
cat("identical to 3 dp:",
    all(round(cls$fixed, 3) == round(cls$free, 3)), "\n")
cat("free disturbance estimate: ",
    sprintf("%.3f (SE %.3f, p = %.3f)",
            fr$est.std[key(fr) == "ResAll ~~ ResAll"][1],
            fr$se[key(fr) == "ResAll ~~ ResAll"][1],
            fr$pvalue[key(fr) == "ResAll ~~ ResAll"][1]), "\n")

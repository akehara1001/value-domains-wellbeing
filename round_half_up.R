# =============================================================================
# round_half_up.R -- the rounding used for every table these scripts write
# =============================================================================
# R's round() and signif() send a value lying exactly half way to the even
# digit (63/112 = 56.25% becomes 56.2 at one decimal). The article rounds half
# up (56.3). Each script that writes tables sources this file near its top;
# from then on round() and signif() in that script round half up, and every
# sprintf("%.Nf") that writes a table cell is given a value already rounded by
# round(x, N), so the printed digits are those of half-up rounding.
#
# A value is treated as exactly half way when it is within 5e-8 of the
# rounding unit of a tie after scaling, which removes binary representation
# noise (2.675 is stored as 2.67499999999999982236431605997495353221893310546875
# and is rounded to 2.68). Package code is not affected: lavaan, poLCA and
# base functions keep base::round.
#
# 00_derive.R (run only on the raw export) and 01_lca.R (the latent class
# models, whose table outputs are written at full precision) do not source it.
# =============================================================================

round_half_up <- function(x, digits = 0) {
  if (is.data.frame(x)) {
    x[] <- lapply(x, round_half_up, digits = digits)
    return(x)
  }
  if (!is.numeric(x)) return(base::round(x, digits))
  p <- 10^digits
  z <- base::round(abs(x) * p, 7)   # binary noise off; 7 decimals of the scaled value
  out <- sign(x) * floor(z + 0.5) / p
  miss <- is.na(x)
  out[miss] <- x[miss]
  out
}

signif_half_up <- function(x, digits = 6) {
  if (!is.numeric(x)) return(base::signif(x, digits))
  digits <- pmax(1, digits)
  ok <- is.finite(x) & x != 0
  e <- ifelse(ok, ceiling(log10(abs(x))), 0)
  out <- round_half_up(x, digits - e)
  out[!ok] <- x[!ok]
  out
}

round  <- round_half_up
signif <- signif_half_up

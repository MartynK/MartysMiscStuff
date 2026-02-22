#' Generate SONOMIND allocations (4 sites) via stratified permuted blocks
#'
#' @param n_total    Total sample size (default 40)
#' @param ratio      Allocation ratios c(High, Low, Sham) (default c(2,1,1))
#' @param block_sz   Block size per stratum (must be a multiple of sum(ratio)); default 4
#' @param site_prob  Numeric vector of length 3 for Sites A–C; Site D = 1 - sum(site_prob)
#' @return           Data.frame with severity, site, and treatment columns
#' @examples
#' alloc <- strat_block_sonomind_4sites(site_prob = c(0.4, 0.3, 0.2))
#' table(alloc$Trt, alloc$site)
#' @export
strat_block_sonomind_4sites <- function(n_total   = 40,
                                        ratio     = c(2, 1, 1),
                                        block_sz  = 4,
                                        site_prob = c(0.4, 0.3, 0.2)) {

  ## 1. Compute site probabilities ----------------------------------
  if (length(site_prob) != 3) {
    stop("site_prob must be length 3")
  }
  p4 <- 1 - sum(site_prob)
  if (p4 < 0) {
    stop("sum(site_prob) must be <= 1")
  }
  all_sites <- c("SiteA", "SiteB", "SiteC", "SiteD")
  probs     <- c(site_prob, p4)

  ## 2. Build covariates data frame --------------------------------
  severity <- factor(sample(
    c("moderate", "severe"),
    size    = n_total,
    replace = TRUE,
    prob    = c(0.6, 0.4)
  ))
  site <- factor(sample(
    all_sites,
    size    = n_total,
    replace = TRUE,
    prob    = probs
  ), levels = all_sites)

  dat <- data.frame(
    severity = severity,
    site     = site,
    stringsAsFactors = TRUE
  )

  ## 3. Run stratified permuted blocks ------------------------------
  library(carat)
  out1 <- StrPBR(
    data  = dat,
    bsize = block_sz
  )$assignments

  out2 <- StrPBR(
    data  = dat,
    bsize = block_sz
  )$assignments

  # recombine; AA and BB are the same, AB and BA are concatenated as SHAM
  out  <- data.frame(
    Trt = factor(
      paste0(out1, out2)
    # ,
    #   levels = c("AA", "AB", "BA", "BB"),
    #   labels = c("High", "Sham", "Sham", "Low")
    ),
    severity = dat$severity,
    site     = dat$site
  )

  return(out)
}


dat <- strat_block_sonomind_4sites()

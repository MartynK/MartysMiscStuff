######################  server.R  ######################
# -----------------------------------------------------
#  Pool ISCVs with a fully‑Bayesian model via cmdstanr
#
#  * Input expected:   per‑study intrasubject CV (unit = proportion, 0–1)
#                      and sample size N (≥3, two‑way crossover → df = N‑2)
#  * Prior on CV    :  CV ~ Beta(2, 2)   (support 0–1 ⇒ 0–100 %)
#  * Likelihood     :  s²_i ~ Inv‑Gamma(df_i/2, df_i * sigma² / 2)
#                      with sigma² = log(1 + CV²)
# -----------------------------------------------------

library(shiny)
library(ggplot2)
library(DT)
library(cmdstanr)
library(here)
options(mc.cores = parallel::detectCores())

# -----------------------------------------------------------------
# 1. Helper: get (and cache) the CmdStan model --------------------
# -----------------------------------------------------------------
get_model <- function() {
  dir.create("inst/iscv_pooling/stan", recursive = TRUE, showWarnings = FALSE)
  stan_src  <- here::here("inst", "iscv_pooling", "stan", "iscv_pool.stan")
  if (!file.exists(stan_src))
    stop("Stan file not found at ", stan_src)

  # cmdstanr handles caching via the locally built executable.
  # compile = FALSE → only compiles if missing.
  cmdstan_model(stan_src)  # compiles only if the executable is missing
}

# -----------------------------------------------------------------
# 2. Shiny server --------------------------------------------------
# -----------------------------------------------------------------
shinyServer(function(input, output, session) {

  # ----- 2·1  dynamic grid of study inputs -----------------------
  output$dynInputs <- renderUI({
    k <- min(max(1, as.integer(input$k)), 20L)
    tagList(lapply(seq_len(k), function(i) {
      fluidRow(
        column(6, numericInput(paste0("N_", i), paste("N", i),  24, min = 3, step = 1)),
        column(6, numericInput(paste0("CV_", i), paste("ISCV", i), 20, min = 1, step = 1))
      )
    }))
  })

  # ----- 2·2  run analysis when 'Update' pressed -----------------
  observeEvent(input$update, {
    # ---- 2·2·a  Collect & validate data -------------------------
    k   <- min(max(1, as.integer(input$k)), 20L)
    Ns  <- sapply(seq_len(k), function(i) as.numeric(input[[paste0("N_", i)]]))
    CVs <- sapply(seq_len(k), function(i) as.numeric(input[[paste0("CV_", i)]]))

    # --- sanity check: did user enter % instead of proportion? ---
    if (any(CVs > 1 & CVs <= 100, na.rm = TRUE)) {
      showNotification("Some CVs look like percentages; converting by ÷100.", type = "warning")
      CVs[CVs > 1 & CVs <= 100] <- CVs[CVs > 1 & CVs <= 100] / 100
    }

    idx <- which(!is.na(Ns) & !is.na(CVs) & Ns > 2 & CVs > 0 & CVs < 1)
    validate(need(length(idx) > 0,
                  "Enter at least one row with N > 2 and 0 < CV < 1 (i.e. 0–100 %)."))
    Ns  <- Ns[idx]
    CVs <- CVs[idx]

    # ---- 2·2·b  Prepare data for Stan ---------------------------
    s2   <- log1p(CVs^2)         # within‑subject variance estimates
    dfs  <- Ns - 2

    stan_dat <- list(M = length(s2), s2 = s2, df = dfs)

    # ---- 2·2·c  Fit Stan model via cmdstanr ---------------------
    n_draw   <- max(1000L, as.integer(input$ndraw))
    model    <- get_model()

    fit <- model$sample(
      data            = stan_dat,
      iter_warmup     = floor(n_draw / 2),
      iter_sampling   = n_draw,
      chains          = 4,
      parallel_chains = 4,
      seed            = 20250428,
      refresh         = 0
    )

    # cmdstanr returns an array (iterations × chains)
    # Extract the first column (one chain concatenated) as a vector
    post <- fit$draws("cv", format = "matrix")[, 1]

    # ---- 2·2·d  Draw prior for overlay --------------------------
    prior_draw <- rbeta(length(post), 2, 2)

    # ---- 2·2·e  Plots -------------------------------------------
    output$priorPlot <- renderPlot({
      ggplot(data.frame(cv = prior_draw), aes(cv)) +
        geom_histogram(aes(y = ..density..), bins = 100, fill = "grey70") +
        labs(title = "Prior:  CV ~ Beta(2,2)", x = "CV", y = "Density") +
        xlim(0, 1) + theme_minimal(base_size = 14)
    })

    output$postPlot <- renderPlot({
      ggplot(data.frame(cv = post), aes(cv)) +
        geom_histogram(aes(y = ..density..), bins = 100, fill = "steelblue", alpha = .75) +
        geom_vline(xintercept = CVs, linetype = "dashed", colour = "red") +
        labs(title = "Posterior distribution", x = "CV", y = "Density") +
        xlim(0, 1) + theme_minimal(base_size = 14)
    })

    # ---- 2·2·f  Summary table ----------------------------------
    qs <- quantile(post, probs = c(.05, .25, .5, .75, .95))
    out <- data.frame(p05 = qs[1], Q1 = qs[2], median = qs[3],
                      mean = mean(post), Q3 = qs[4], p95 = qs[5])
    output$postSum <- renderDT({
      datatable(round(out, 4), options = list(dom = "t"))
    })
  })
})


library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

# ---------------------------------------------------------------------------
# Create a tidy little 2×2 crossover example -------------------------------
# ---------------------------------------------------------------------------
set.seed(42)
n  <- 24                        # 24 subjects
mu <- c(A = 10, B = 12)         # mean response by treatment
sd_between <- 4                 # true between-subject SD
sd_within  <- 1                 # true within-subject SD

subjects <- tibble::tibble(
  subject = factor(seq_len(n)),
  period1 = sample(c("A", "B"), n, replace = TRUE)
) |>
  dplyr::mutate(
    period2 = ifelse(period1 == "A", "B", "A")
  ) |>
  tidyr::pivot_longer(
    cols      = starts_with("period"),
    names_to  = "period",
    values_to = "treatment",
    names_pattern = "period(\\d)"
  ) |>
  dplyr::arrange(subject, period) |>
  dplyr::group_by(subject) |>
  dplyr::mutate(
    # Between-subject random effect
    b_i = rnorm(1, 0, sd_between)
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    # Response = μ_treatment + between + within
    response = mu[treatment] + b_i + rnorm(dplyr::n(), 0, sd_within),
    period   = factor(period, levels = c("1", "2"))
  )



g1 <- ggplot(subjects,
             aes(x = period, y = response, group = subject)) +
  geom_point(size = 2) +
  geom_line(alpha = .4) +
  theme_classic() +
  labs(title = "Option 1 – Paired dot-plot",
       x = "Period", y = "Observed response")
print(g1)


##### Option 2


changes <- subjects %>%
  # keep only the bits we need for the delta
  select(subject, period, response) %>%
  pivot_wider(
    id_cols     = subject,      # <-- only subject groups the two rows
    names_from  = period,       # will give you columns `1` and `2`
    values_from = response
  ) %>%
  mutate(
    delta          = `2` - `1`,
    delta_centered = delta - mean(delta, na.rm = TRUE)
  )

# Left: demeaned raw
g_left <- subjects |>
  mutate(response_c = response - mean(response)) |>
  ggplot(aes(x = period, y = response_c, group = subject)) +
  geom_point(size = 2) +
  geom_line(alpha = .4) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_classic() +
  labs(x = "Period", y = "Centered observation")

# Right: changes
g_right <- ggplot(changes, aes(x = 1, y = delta)) +
  geom_point( size = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "Change in observation")

(g_left | g_right) +
  patchwork::plot_annotation(
    title = "Option 2 – Paired plot with centred change column")


##### Option 3

g3 <- subjects |>
  mutate(response_c = response - mean(response)) |>
  ggplot(aes(x = period, y = response_c, group = subject)) +
  geom_point(size = 2) +
  gghalves::geom_half_violin(position = position_nudge(x = -.05),
                             side = "l", alpha = .5, mapping = aes(group=NULL)) +
  geom_line(aes(group = subject), color="salmon4", alpha = .5) +
  geom_line(alpha = .4) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_classic() +
  labs(x = "Treatment", y = "Centered observations (2/subject)")


g3_right <- ggplot(changes, aes(x = 1, y = delta)) +
  geom_point( size = 2, color = "salmon4") +
  geom_hline(yintercept = 0, linetype = 2) +
  gghalves::geom_half_violin(position = position_nudge(x = -.05),
                             side = "l", alpha = .5) +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "Change in observation per subject")

(g3 | g3_right) +
  patchwork::plot_annotation(
    title = "Variability in a 2x2 crossover trial - changes within subject may have lower variability")


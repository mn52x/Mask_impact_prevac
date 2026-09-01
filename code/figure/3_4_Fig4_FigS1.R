suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(cowplot)
  library(posterior); library(stringr); library(ggdist); library(showtext)
})
showtext_auto()

base_pt      <- 10
country_name <- "JPN"

col_averted   <- "#d55e00"    
col_avertible <- "#0072b2"  

# =================
# read fit 
# 1-1 get standata from 2_2_prepare_runstan
# 1-2 run "3_1_read_csv.R" 
#=====================  

#===================
# Fig 3
#===================
#baseline result
fit_base <- make_csv_fit(
    dir("../JPN", 
        pattern    = "2_2_prepare_runstan-XXX-.*\\.csv$",
        full.names = TRUE))

#lower reporting coverage
fit_report_low <- make_csv_fit(
    dir("../SA_JPN_reportlow",
        pattern    = "2_2_prepare_runstan-XXX-.*\\.csv$",
        full.names = TRUE))

#higher reporting coverage 
fit_report_high <- make_csv_fit(
    dir("../SA_JPN_reporthigh",
        pattern    = "2_2_prepare_runstan-XXX-.*\\.csv$",
        full.names = TRUE))

#mask-effect lower
fit_exp_low <- make_csv_fit(
    dir("../SA_JPN_exp_lower",
        pattern    = "2_2_prepare_runstan-XXX-.*\\.csv$",
        full.names = TRUE))

#mask-effect upper 
fit_exp_high <- make_csv_fit(
    dir("../SA_JPN_exp_upper",
        pattern    = "2_2_prepare_runstan-XXX-.*\\.csv$",
        full.names = TRUE))
  

# -----------------------------------------------------------------------------
#　averted + avertible
# -----------------------------------------------------------------------------
# scenario_label: "base" / "CF1" / "CF2"
make_M_long <- function(fit, varname, scenario_label){
  dm <- fit$draws(varname, format = "draws_matrix")
  df <- as.data.frame(dm)  # 列名は M_hat_out[1,1] など
  
  df |>
    mutate(.draw = dplyr::row_number()) |>

    pivot_longer(
      cols = - .draw,
      names_to = "var",
      values_to = "mu_cases"  
    ) |>
    mutate(
      scenario = scenario_label,
      # 例: "M_hat_out[12,3]" -> t = 12, age = 3
      t   = as.integer(stringr::str_match(var, "\\[(\\d+),(\\d+)\\]")[, 2]),
      age = as.integer(stringr::str_match(var, "\\[(\\d+),(\\d+)\\]")[, 3])
    ) |>
    select(.draw, scenario, t, age, mu_cases)
}

calc_diff_summary <- function(fit, rate_label,
                              cum_start, cum_end, pop,
                              prefix_base = "M_hat_out",
                              prefix_cf1  = "M_hat_CF1_out",
                              prefix_cf2  = "M_hat_CF2_out") {
  M_base <- make_M_long(fit, prefix_base, "base")
  M_cf1  <- make_M_long(fit, prefix_cf1,  "CF1")
  M_cf2  <- make_M_long(fit, prefix_cf2,  "CF2")
  M_long <- bind_rows(M_base, M_cf1, M_cf2)
  rm(M_base, M_cf1, M_cf2); gc()

  cum_draws <- M_long %>%
    filter(t >= cum_start, t <= cum_end) %>%
    group_by(.draw, scenario) %>%
    summarise(cases_cum = sum(mu_cases, na.rm = TRUE), .groups = "drop") %>%  # ★修正
    pivot_wider(names_from = scenario, values_from = cases_cum) %>%
    mutate(
      averted_per100k   = (CF2  - base) / pop * 1e5,
      avertible_per100k = (base - CF1)  / pop * 1e5
    )

  bind_rows(
    cum_draws %>%
      select(.draw, value = averted_per100k) %>%
      mutate(quantity = "Averted"),
    cum_draws %>%
      select(.draw, value = avertible_per100k) %>%
      mutate(quantity = "Avertible")
  ) %>%
    mutate(rate = rate_label)
}

diff_draws_mhat <- bind_rows(
  calc_diff_summary(fit_report_low,  "Low reporting coverage (12.5%)",
                    cum_start, cum_end, pop_total),
  calc_diff_summary(fit_exp_low,     "Mask-wearing effect\n95% CrI lower (6%)",
                    cum_start, cum_end, pop_total),
  calc_diff_summary(fit_base,        "Base",
                    cum_start, cum_end, pop_total),
  calc_diff_summary(fit_report_high, "High reporting coverage (50%)",
                    cum_start, cum_end, pop_total),
  calc_diff_summary(fit_exp_high,    "Mask-wearing effect\n95% CrI upper (43%)",
                    cum_start, cum_end, pop_total)
) %>%
  mutate(
    rate = factor(rate, levels = c(
      "High reporting coverage (50%)",
      "Mask-wearing effect\n95% CrI upper (43%)",
      "Base",
      "Low reporting coverage (12.5%)",
      "Mask-wearing effect\n95% CrI lower (6%)"
    )),
    quantity = factor(quantity, levels = c("Avertible", "Averted")),
    source   = "M_hat"
  )
saveRDS(diff_draws_mhat, "fig4_diff_draws_mhat_effect.rds")
diff_draws_mhat <- readRDS("fig4_diff_draws_mhat_effect.rds")

# -----------------------------------------------------------------------------
# 6. A：averted 
# -----------------------------------------------------------------------------
diff_averted <- diff_draws_mhat %>%
  filter(quantity == "Averted")

p_halfeye_averted <- ggplot(diff_averted,
                            aes(x = value, y = rate)) +
  stat_halfeye(
    fill           = col_averted,
    color          = col_averted,
    .width         = c(0.50, 0.95),
    point_interval = "median_qi",
    slab_alpha     = 0.4,
    interval_size_range = c(0.4, 0.8),
    point_size     = 0.2,
    normalize      = "groups",
    scale = 0.6
  ) +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_discrete(limits = rev(levels(diff_averted$rate))) +
  labs(
    x = "Cumulative cases averted per 100,000 people",
    y = "Parameter values"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text.x    = element_text(size = base_pt * 0.48, color = "black"),
    axis.text.y   = element_text(size = base_pt * 0.5, color = "black"),
    axis.title.x = element_text(size = base_pt * 0.6, color = "black",
                                margin = margin(t = 6)),
    axis.title.y = element_text(size = base_pt * 0.6, color = "black",
                                margin = margin(r = 6)),
    axis.line    = element_line(color = "black", linewidth = 0.2),
    axis.ticks   = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 15, r = 6, b = 5, l = 10)
  )

# -----------------------------------------------------------------------------
# B：avertible
# -----------------------------------------------------------------------------
diff_avertible <- diff_draws_mhat %>%
  filter(quantity == "Avertible")

p_halfeye_avertible <- ggplot(diff_avertible,
                            aes(x = value, y = rate)) +
  stat_halfeye(
    fill           = col_avertible,
    color          = col_avertible,
    .width         = c(0.50, 0.95),
    interval_size_range = c(0.4, 0.8),
    point_interval = "median_qi",
    slab_alpha     = 0.4,
    point_size     = 0.2,
    normalize      = "groups",
    scale = 0.6
  ) +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_discrete(limits = rev(levels(diff_averted$rate))) +
  labs(
    x = "Cumulative cases avertible per 100,000 people",
    y = NULL
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text.x    = element_text(size = base_pt * 0.48, color = "black"),
    axis.text.y  = element_blank(),
    axis.title.x = element_text(size = base_pt * 0.6, color = "black",
                                margin = margin(t = 6)),
    axis.title.y = element_blank(),
    axis.line    = element_line(color = "black", linewidth = 0.2),
    axis.ticks   = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 15, r = 8, b = 5, l = 3)
  )

# --- merge ---
rw <- c(1.7, 1)
boundary <- rw[1] / sum(rw)

p_fig_4 <- cowplot::plot_grid(p_halfeye_averted, p_halfeye_avertible,
                             ncol = 2, rel_widths = rw,
                             align = "h"
) + theme(plot.margin = margin(0, 1, 0, 0))

p_fig_4 <- cowplot::ggdraw(p_fig_4) +
  cowplot::draw_plot_label(
    label    = c("A", "B"),
    x        = c(0.25, 0.62),
    y        = c(0.97, 0.97),
    size     = base_pt * 0.75,
    fontface = "bold",
    family   = "Arial",
    hjust    = c(0, 0),
    vjust    = c(1, 1)
  )


#===================
# FigS4
#==================
# serial inteval using lognormal distribution
fit_SI_lognormal <- make_csv_fit(
    dir("../SA_JPN_lognormal",
        pattern    = "2_2_prepare_runstan-XXX--.*\\.csv$",
        full.names = TRUE))
 # mask-wearing coverage using CTIS data 
fit_mask_CTIS <- make_csv_fit(
    dir("../SA_JPN_maskFB",
        pattern    = "2_2_prepare_runstan-XXX--.*\\.csv$",
        full.names = TRUE))
#wider prior distribution for MCMC sampling   
fit_widerprior <- make_csv_fit(
    dir("../SA_JPN_widerprior",
        pattern    = "2_2_prepare_runstan-XXX--.*\\.csv$",
        full.names = TRUE))

#lower wearer protection
fit_lowlambda <- make_csv_fit(
    dir("../SA_JPN_lambdalower",
        pattern    = "2_2_prepare_runstan-XXX--.*\\.csv$",
        full.names = TRUE))
#swapping wearer protection and source control effect values  
fit_swaplamda<- make_csv_fit(
    dir("../SA_JPN_switchlambda",
        pattern    = "2_2_prepare_runstan-XXX--.*\\.csv$",
        full.names = TRUE))

# -----------------------------------------------------------------------------
#　averted + avertible 
# -----------------------------------------------------------------------------
# scenario_label: "base" / "CF1" / "CF2"
make_M_long <- function(fit, varname, scenario_label){
  dm <- fit$draws(varname, format = "draws_matrix")
  df <- as.data.frame(dm) 
  
  df |>
    mutate(.draw = dplyr::row_number()) |>
    pivot_longer(
      cols = - .draw,
      names_to = "var",
      values_to = "mu_cases"  
    ) |>
    mutate(
      scenario = scenario_label,
      # 例: "M_hat_out[12,3]" -> t = 12, age = 3
      t   = as.integer(stringr::str_match(var, "\\[(\\d+),(\\d+)\\]")[, 2]),
      age = as.integer(stringr::str_match(var, "\\[(\\d+),(\\d+)\\]")[, 3])
    ) |>
    select(.draw, scenario, t, age, mu_cases)
}


calc_diff_summary <- function(fit, rate_label,
                              cum_start, cum_end, pop,
                              prefix_base = "M_hat_out",
                              prefix_cf1  = "M_hat_CF1_out",
                              prefix_cf2  = "M_hat_CF2_out") {
  M_base <- make_M_long(fit, prefix_base, "base")
  M_cf1  <- make_M_long(fit, prefix_cf1,  "CF1")
  M_cf2  <- make_M_long(fit, prefix_cf2,  "CF2")
  M_long <- bind_rows(M_base, M_cf1, M_cf2)
  rm(M_base, M_cf1, M_cf2); gc()

  cum_draws <- M_long %>%
    filter(t >= cum_start, t <= cum_end) %>%
    group_by(.draw, scenario) %>%
    summarise(cases_cum = sum(mu_cases, na.rm = TRUE), .groups = "drop") %>%  
    pivot_wider(names_from = scenario, values_from = cases_cum) %>%
    mutate(
      averted_per100k   = (CF2  - base) / pop * 1e5,
      avertible_per100k = (base - CF1)  / pop * 1e5
    )

  bind_rows(
    cum_draws %>%
      select(.draw, value = averted_per100k) %>%
      mutate(quantity = "Averted"),
    cum_draws %>%
      select(.draw, value = avertible_per100k) %>%
      mutate(quantity = "Avertible")
  ) %>%
    mutate(rate = rate_label)
}

diff_draws_mhat <- bind_rows(
  calc_diff_summary(fit_base,        "Base",
                    cum_start, cum_end, pop_total),
  calc_diff_summary(fit_mask_CTIS,     "Mask-wearing CTIS survey",
                    cum_start, cum_end, pop_total),
  calc_diff_summary(fit_lowlambda,     "Lower wearer protection (0.86)",
                    cum_start, cum_end, pop_total),
  calc_diff_summary(fit_swaplamda,     "Swapping two pathway values\nof mask wearing",
                    cum_start, cum_end, pop_total),
  calc_diff_summary(fit_SI_lognormal,   "Serial intervel\nlognormal distribution",
                    cum_start, cum_end, pop_total),
  calc_diff_summary(fit_widerprior, "Wider prior distribution",
                    cum_start, cum_end, pop_total),
  
) %>%
  mutate(
    rate = factor(rate, levels = c(
      "Base",
      "Mask-wearing CTIS survey",
      "Lower wearer protection (0.86)",
      "Swapping two pathway values\nof mask wearing",
      "Serial intervel\nlognormal distribution",
      "Wider prior distribution"
    )),
    quantity = factor(quantity, levels = c("Avertible", "Averted")),
    source   = "M_hat"
  )



saveRDS(diff_draws_mhat, "figS2_diff_draws.rds")

 diff_draws_mhat <- readRDS("figS2_diff_draws.rds")
 

# -----------------------------------------------------------------------------
#  A：averted 
# -----------------------------------------------------------------------------
diff_averted <- diff_draws_mhat %>%
  filter(quantity == "Averted")

p_halfeye_averted <- ggplot(diff_averted,
                            aes(x = value, y = rate)) +
  stat_halfeye(
    fill           = col_averted,
    color          = col_averted,
    .width         = c(0.50, 0.95),
    point_interval = "median_qi",
    slab_alpha     = 0.3,
    interval_size  = 0.6,
    point_size     = 0.3,
    normalize      = "groups",
    scale = 0.5
  ) +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_discrete(limits = rev(levels(diff_averted$rate))) +
  labs(
    x = "Cumulative cases averted per 100,000 people",
    y = "Parameter values"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text.x    = element_text(size = base_pt * 0.55, color = "black"),
    axis.text.y   = element_text(size = base_pt * 0.55, color = "black"),
    axis.title.x = element_text(size = base_pt * 0.6, color = "black",
                                margin = margin(t = 6)),
    axis.title.y = element_text(size = base_pt * 0.7, color = "black",
                                margin = margin(r = 6)),
    axis.line    = element_line(color = "black", linewidth = 0.2),
    axis.ticks   = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 15, r = 8, b = 5, l = 3)
  )

# -----------------------------------------------------------------------------
#  B：avertible
# -----------------------------------------------------------------------------
diff_avertible <- diff_draws_mhat %>%
  filter(quantity == "Avertible")

p_halfeye_avertible <- ggplot(diff_avertible,
                            aes(x = value, y = rate)) +
  stat_halfeye(
    fill           = col_avertible,
    color          = col_avertible,
    .width         = c(0.50, 0.95),
    point_interval = "median_qi",
    slab_alpha     = 0.3,
    interval_size  = 0.6,
    point_size     = 0.3,
    normalize      = "groups",
    scale = 0.5
  ) +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_discrete(limits = rev(levels(diff_averted$rate))) +
  labs(
    x = "Cumulative cases avertible per 100,000 people",
    y = NULL
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text.x    = element_text(size = base_pt * 0.55, color = "black"),
    axis.text.y  = element_blank(),
    axis.title.x = element_text(size = base_pt * 0.6, color = "black",
                                margin = margin(t = 6)),
    axis.title.y = element_blank(),
    axis.line    = element_line(color = "black", linewidth = 0.2),
    axis.ticks   = element_line(color = "black", linewidth = 0.2),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 15, r = 8, b = 5, l = 3)
  )

# ---merge --
rw <- c(1.5, 1)
boundary <- rw[1] / sum(rw)

p_fig_4 <- cowplot::plot_grid(p_halfeye_averted, p_halfeye_avertible,
                             ncol = 2, rel_widths = rw)

p_fig_4 <- cowplot::ggdraw(p_fig_4) +
  cowplot::draw_plot_label(
    label    = c("A", "B"),
    x        = c(0.23, 0.62),
    y        = c(0.99, 0.99),
    size     = base_pt * 0.75,
    fontface = "bold",
    family   = "Arial",
    hjust    = c(0, 0),
    vjust    = c(1, 1)
  )
p_fig_4

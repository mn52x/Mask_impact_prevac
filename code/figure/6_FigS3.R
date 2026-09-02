suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(posterior)
})
# ============================================================
# setting 
# 1. Run 2_2_prepare_runstan.r to get observed data
# 2. Run 3_1_read_csv.r to get fit data
# ============================================================
# ============================================================
# 1)utility 
# ============================================================
parse_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  suppressWarnings(as.Date(x, tryFormats = c(
    "%Y-%m-%d", "%Y/%m/%d", "%Y%m%d",
    "%m/%d/%Y", "%d/%m/%Y", "%d.%m.%Y"
  )))
}

rollmaybe <- function(x, window = NULL, align = "center") {
  if (is.null(window) || window <= 1) return(x)
  zoo::rollapply(x, width = window, FUN = mean, partial = TRUE, align = align)
}

# ============================================================
# 2) Observed data function
# ============================================================
make_obs_total_df <- function(data, window = NULL, align = "center",
                              date_col = "date",
                              age_cols = c("X19","X2039","X4059","X60over")) {
  stopifnot(all(c(date_col, age_cols) %in% names(data)))
  tibble::tibble(
    t        = seq_len(nrow(data)),
    date     = parse_date_safe(data[[date_col]]),
    observed = rollmaybe(rowSums(as.data.frame(data[age_cols])), window, align)
  )
}

make_obs_age_long <- function(data, window = NULL, align = "center",
                              date_col = "date",
                              age_cols = c("X19","X2039","X4059","X60over")) {
  stopifnot(all(c(date_col, age_cols) %in% names(data)))
  smoothed <- lapply(age_cols, function(col) rollmaybe(data[[col]], window, align))
  names(smoothed) <- age_cols
  dplyr::bind_cols(
    tibble::tibble(
      t    = seq_len(nrow(data)),
      date = parse_date_safe(data[[date_col]])
    ),
    tibble::as_tibble(smoothed)
  ) |>
    tidyr::pivot_longer(-c(t, date), names_to = "age", values_to = "observed")
}

# ============================================================
# 3) Posterior summary CI, Rt
# ============================================================
make_ci_total <- function(fit, prefix = "M_hat_out", window = NULL, align = "center") {
  fit$draws(prefix) |>
    posterior::as_draws_df() |>
    pivot_longer(starts_with(prefix),
                 names_to    = c("t","age"),
                 names_pattern = paste0("^", prefix, "\\[(\\d+),(\\d+)\\]$"),
                 values_to   = "mu") |>
    mutate(t = as.integer(t), age = as.integer(age)) |>
    group_by(.draw, t) |>
    summarise(mu_tot = sum(mu), .groups = "drop") |>
    arrange(.draw, t) |>
    group_by(.draw) |>
    mutate(mu_smooth = rollmaybe(mu_tot, window, align)) |>
    ungroup() |>
    group_by(t) |>
    summarise(
      med  = median(mu_smooth, na.rm = TRUE),
      low  = quantile(mu_smooth, 0.025, na.rm = TRUE),
      high = quantile(mu_smooth, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

make_ci_age <- function(fit, prefix = "M_hat_out", window = NULL, align = "center") {
  fit$draws(prefix) |>
    posterior::as_draws_df() |>
    pivot_longer(starts_with(prefix),
                 names_to    = c("t","age"),
                 names_pattern = paste0("^", prefix, "\\[(\\d+),(\\d+)\\]$"),
                 values_to   = "mu") |>
    mutate(t = as.integer(t), age = as.integer(age)) |>
    group_by(.draw, age) |>
    arrange(t, .by_group = TRUE) |>
    mutate(mu_smooth = rollmaybe(mu, window, align)) |>
    ungroup() |>
    group_by(t, age) |>
    summarise(
      med  = median(mu_smooth, na.rm = TRUE),
      low  = quantile(mu_smooth, 0.025, na.rm = TRUE),
      high = quantile(mu_smooth, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

make_Rt <- function(fit, rt_var = "Rt", window = NULL, align = "center") {
  fit$draws(rt_var) |>
    posterior::as_draws_df() |>
    pivot_longer(starts_with(rt_var),
                 names_to    = "t",
                 names_pattern = paste0("^", rt_var, "\\[(\\d+)\\]$"),
                 values_to   = "Rt") |>
    mutate(t = as.integer(t)) |>
    group_by(.draw) |>
    arrange(t, .by_group = TRUE) |>
    mutate(Rt_smooth = rollmaybe(Rt, window, align)) |>
    ungroup() |>
    group_by(t) |>
    summarise(
      med  = median(Rt_smooth, na.rm = TRUE),
      low  = quantile(Rt_smooth, 0.025, na.rm = TRUE),
      high = quantile(Rt_smooth, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

# ============================================================
# 4) save RDS
# ============================================================
save_country_summaries <- function(fit,
                                   data,
                                   fix_end,
                                   country_code,
                                   processed_root = "processed_supplementary",
                                   window_Rt = 7,
                                   window_inc = NULL,
                                   overwrite = FALSE) {
  
  out_dir <- file.path(processed_root, country_code)

  done_marker <- file.path(out_dir, "_done.txt")
  if (file.exists(done_marker) && !overwrite) {
    message(country_code, ": done (skip)")
    return(invisible(NULL))
  }
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  message(country_code, ": start...")
  
  # CI total
  message("  - CI total (Fit)")
  saveRDS(make_ci_total(fit, "M_hat_out",     window_inc),
          file.path(out_dir, "ci_total_fit.rds"))
  
  message("  - CI total (CF1)")
  saveRDS(make_ci_total(fit, "M_hat_CF1_out", window_inc),
          file.path(out_dir, "ci_total_cf1.rds"))
  
  message("  - CI total (CF2)")
  saveRDS(make_ci_total(fit, "M_hat_CF2_out", window_inc),
          file.path(out_dir, "ci_total_cf2.rds"))
  
  # CI age
  message("  - CI age (Fit)")
  saveRDS(make_ci_age(fit, "M_hat_out",     window_inc),
          file.path(out_dir, "ci_age_fit.rds"))
  
  message("  - CI age (CF1)")
  saveRDS(make_ci_age(fit, "M_hat_CF1_out", window_inc),
          file.path(out_dir, "ci_age_cf1.rds"))
  
  message("  - CI age (CF2)")
  saveRDS(make_ci_age(fit, "M_hat_CF2_out", window_inc),
          file.path(out_dir, "ci_age_cf2.rds"))
  
  # Rt
  message("  - Rt (Fit)")
  saveRDS(make_Rt(fit, "Rt",     window_Rt),
          file.path(out_dir, "rt_fit.rds"))
  
  message("  - Rt (CF1)")
  saveRDS(make_Rt(fit, "Rt_CF1", window_Rt),
          file.path(out_dir, "rt_cf1.rds"))
  
  message("  - Rt (CF2)")
  saveRDS(make_Rt(fit, "Rt_CF2", window_Rt),
          file.path(out_dir, "rt_cf2.rds"))
  
  # Observed data
  message("  - Observed data")
  saveRDS(make_obs_total_df(data, window = NULL),
          file.path(out_dir, "obs_tot.rds"))
  saveRDS(make_obs_age_long(data, window = NULL),
          file.path(out_dir, "obs_age.rds"))
  
  # Meta info
  saveRDS(list(country_code = country_code,
               fix_end      = fix_end,
               window_Rt    = window_Rt,
               window_inc   = window_inc,
               saved_at     = Sys.time()),
          file.path(out_dir, "meta.rds"))
  
  writeLines(as.character(Sys.time()), done_marker)
  message(country_code, ": done → ", out_dir)
  invisible(out_dir)
}

# ============================================================
# 5) save fit/data/fix_end 
# ============================================================
processed_root <- "../processed"

# country name
country_code <- country_name
save_country_summaries(
  fit            = fit,
  data           = data,
  fix_end        = fix_end,
  country_code   = country_code,
  processed_root = processed_root,
  window_Rt      = 7,
  window_inc     = NULL,
  overwrite      = TRUE   
)

gc()

# ============================================================
# ★ 
# ============================================================
BASE_SIZE   <- 8
BASE_FAMILY <- "Arial"   

size_axis_text   <- BASE_SIZE * 0.8
size_axis_title  <- BASE_SIZE * 0.79
size_strip_text  <- BASE_SIZE * 0.68
size_plot_title  <- BASE_SIZE * 0.79
size_plot_tag    <- BASE_SIZE * 1
size_legend_text <- BASE_SIZE * 0.78
size_main_title  <- BASE_SIZE * 1.1

processed_root <- "../processed"
out_root <- "../supplementary_figures"

# ============================================================
# 1)
# ============================================================
pt_to_mm <- function(pt) pt * 0.3528

win_label <- function(w) {
  if (is.null(w) || w <= 1) "Daily COVID-19 cases"
  else sprintf("COVID-19 cases (%d-day MA)", w)
}

fix_date_safe <- function(df, fix_end) {
  if (!is.numeric(fix_end) || length(fix_end) != 1) return(NA)
  if (fix_end < 1 || fix_end > nrow(df)) return(NA)
  df$date[fix_end]
}

fix_date_from_t <- function(df_long_or_wide, fix_end) {
  if (!is.numeric(fix_end) || length(fix_end) != 1) return(as.Date(NA))
  dd <- df_long_or_wide |>
    dplyr::select(t, date) |>
    dplyr::distinct(t, .keep_all = TRUE) |>
    dplyr::arrange(t)
  if (fix_end < 1 || fix_end > nrow(dd)) return(as.Date(NA))
  dd$date[fix_end]
}

scen_colors_fill <- c(
  "Fit" = "#ffd700",  
  "CF1" = "#65bbe9",   
  "CF2" = "#EF9F27"    
)

scen_colors_rt <- c(
  "Fit" = "#7a6a1c",
  "CF1" = "#0E4F8F",
  "CF2" = "#b4531e"
)

scen_labels <- c(
  "Fit" = "Actual",
  "CF2" = "Lower Mask: Asia 1.8%, Europe 0.2%",
  "CF1" = "Mask: 97%"
)

country_name_map <- c(
  AUS = "Australia",  BEL = "Belgium",       DEU = "Germany",
  DNK = "Denmark",    ESP = "Spain",         FRA = "France",
  GBR = "United Kingdom", HKG = "Hong Kong", ITA = "Italy",
  JPN = "Japan",      KOR = "South Korea",   MYS = "Malaysia",
  NLD = "Netherlands", PRT = "Portugal",     THA = "Thailand",
  TWN = "Taiwan"
)
country_cf2_rate <- c(
  # Asia (1.8%)
  AUS = 1.8, HKG = 1.8, JPN = 1.8, KOR = 1.8, MYS = 1.8, THA = 1.8, TWN = 1.8,
  # Europe (0.2%)
  BEL = 0.2, DEU = 0.2, DNK = 0.2, ESP = 0.2, FRA = 0.2,
  GBR = 0.2, ITA = 0.2, NLD = 0.2, PRT = 0.2 
)

XLIM_COMMON   <- as.Date(c("2020-02-01", "2021-03-31"))
XBREAKS_COMMON <- seq(as.Date("2020-02-01"), as.Date("2021-02-01"), by = "4 months")
                     
# ============================================================
# 1)read RDS
# ============================================================
load_country_summaries <- function(country_code,
                                   processed_root = "processed") {
  
  d <- file.path(processed_root, country_code)
  if (!dir.exists(d)) {
    stop("Processed directory not found for ", country_code, ": ", d)
  }
  
  list(
    country_code = country_code,
    ci_total_fit = readRDS(file.path(d, "ci_total_fit.rds")),
    ci_total_cf1 = readRDS(file.path(d, "ci_total_cf1.rds")),
    ci_total_cf2 = readRDS(file.path(d, "ci_total_cf2.rds")),
    ci_age_fit   = readRDS(file.path(d, "ci_age_fit.rds")),
    ci_age_cf1   = readRDS(file.path(d, "ci_age_cf1.rds")),
    ci_age_cf2   = readRDS(file.path(d, "ci_age_cf2.rds")),
    rt_fit       = readRDS(file.path(d, "rt_fit.rds")),
    rt_cf1       = readRDS(file.path(d, "rt_cf1.rds")),
    rt_cf2       = readRDS(file.path(d, "rt_cf2.rds")),
    obs_tot      = readRDS(file.path(d, "obs_tot.rds")),
    obs_age      = readRDS(file.path(d, "obs_age.rds")),
    meta         = readRDS(file.path(d, "meta.rds"))
  )
}

theme_supp <- function() {
  ggplot2::theme(
    text        = ggplot2::element_text(family = BASE_FAMILY),
    axis.text   = ggplot2::element_text(size = size_axis_text),
    axis.title  = ggplot2::element_text(size = size_axis_title),
    strip.text  = ggplot2::element_text(size = size_strip_text),
    plot.title  = ggplot2::element_text(size = size_plot_title),
    plot.tag    = ggplot2::element_text(size = size_plot_tag, face = "bold"),
    #legend.text = ggplot2::element_text(size = size_legend_text),
    plot.margin = ggplot2::margin(8, 12, 8, 12)
  )
}

theme_submit_like <- function(axis_lw_pt = 0.6,
                              tick_pt = 0.6,
                              tick_len_pt = 3) {
  ggplot2::theme_classic(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line  = ggplot2::element_line(color = "black",
                                         linewidth = pt_to_mm(axis_lw_pt)),
      axis.ticks = ggplot2::element_line(color = "black",
                                         linewidth = pt_to_mm(tick_pt)),
      axis.ticks.length = grid::unit(tick_len_pt, "pt"),
      plot.title = ggplot2::element_text(face = "plain"),
      plot.margin = ggplot2::margin(6, 8, 6, 6, unit = "pt")
    )
}

x_tweak_total <- function() list(
  ggplot2::scale_x_date(
    limits = XLIM_COMMON,
    breaks = XBREAKS_COMMON,
    expand = ggplot2::expansion(mult = c(0.02, 0.02)),
    labels = scales::label_date("%b\n%Y", locale = "en")
  )
)

age_x_tweak <- function() list(
  ggplot2::scale_x_date(
    limits = XLIM_COMMON,
    breaks = XBREAKS_COMMON,
    expand = ggplot2::expansion(mult = c(0.02, 0.02)),
    labels = scales::label_date("%b\n%Y", locale = "en")
  ),
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = size_axis_text * 0.85,
                                        lineheight = 0.9)
  )
)

age_strip_thin <- function(strip_lw_pt = 0.5) {
  ggplot2::theme(
    strip.background = ggplot2::element_rect(
      fill      = "white",
      color     = "black",
      linewidth = pt_to_mm(strip_lw_pt)
    ),
    strip.text = ggplot2::element_text(
      size   = size_strip_text,
      margin = ggplot2::margin(t = 2, b = 2)
    )
  )
}

plot_total_CI_scen_from_rds <- function(ci_total_df, rt_df, obs_df, fix_end,
                                        scenario = "Fit",
                                        window = NULL,
                                        rt_max_display = 4,
                                        y_max_override = NULL,
                                        show_obs = TRUE) {
  
  df <- ci_total_df |>
    dplyr::left_join(obs_df, by = "t") |>
    dplyr::left_join(rt_df, by = "t", suffix = c("", "_Rt")) |>
    dplyr::filter(!is.na(date)) |>
    dplyr::mutate(date = as.Date(date))
  
  y_max <- if (is.null(y_max_override)) max(df$high, na.rm = TRUE) else y_max_override
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  scale_factor <- y_max / max(rt_max_display, 1)
  
  vline_x <- fix_date_safe(obs_df, fix_end)
  
  col_ci <- scen_colors_fill[scenario]
  col_rt <- scen_colors_rt[scenario]
  
  # ★凡例用のキー→色の対応
  key_inc <- "Estimated"
  key_rt  <- "Rt"
  key_obs <- "Reported"
  
  #lgd_cols <- setNames(c(col_ci, col_rt, "gray10"),
                      # c(key_inc, key_rt, key_obs))
  keys <- c(key_inc, key_rt, if (show_obs) key_obs)  
  
  lgd_cols <- setNames(c(col_ci, col_rt, "gray10")[seq_along(keys)], keys)
  
  ov <- list(
    linetype  = ifelse(keys == key_obs, "blank", "solid"),
    shape     = ifelse(keys == key_obs, 16, NA),
    size      = ifelse(keys == key_obs, 0.6, 0),
    linewidth = ifelse(keys == key_obs, 0, 0.4)   
    #fill      = scales::alpha(lgd_cols, c(0.4, 0.3, 0)[seq_along(keys)]),  
    #alpha     = 1                               
  )
  
  lgd <- guide_legend(title = NULL, override.aes = ov)
  
  g <- ggplot(df, aes(x = date)) +
    geom_ribbon(aes(ymin = low, ymax = high, fill = key_inc),   
                alpha = 0.4) +
    geom_line(aes(y = med, color = key_inc),                  
              linewidth = 0.25)
  
  if (show_obs) {
    g <- g + geom_point(aes(y = observed, color = key_obs),     
                        shape = 16,
                        size = 0.35, stroke = 0, alpha = 0.6)
  }
  if (!is.na(vline_x)) {
    g <- g + geom_vline(xintercept = vline_x,
                        linetype = "dashed", color = "#009e73",
                        linewidth = 0.2)
  }
  
  g +
    geom_ribbon(aes(ymin = low_Rt * scale_factor,
                    ymax = high_Rt * scale_factor,
                    fill = key_rt),                             
                alpha = 0.3) +
    
    geom_line(aes(y = med_Rt * scale_factor, color = key_rt),   
              linewidth = 0.3, linetype = "solid") +
    
    geom_hline(yintercept = 1 * scale_factor,
               linetype = "dashed", color = "grey20",
               linewidth = 0.2) +
    
    scale_color_manual(values = lgd_cols, breaks = keys,
                       limits = keys, name = "lgd") +
    scale_fill_manual(values  = lgd_cols, breaks = keys,
                      limits = keys, name = "lgd") +
    guides(color = lgd, fill = lgd) +
    scale_y_continuous(
      name     = win_label(window),
      limits   = c(0, y_max),
      labels   = scales::label_number(accuracy = 1, big.mark = ","),
      breaks = function(limits) {
        hi <- limits[2]
        if (!is.finite(hi) || hi <= 0) return(c(0))
        step <- scales::breaks_pretty(n = 3)(c(0, hi))
        step <- diff(step)[1]
        seq(0, hi, by = step)
      },
      sec.axis = sec_axis(~ . / scale_factor,
                          name = "Effective reproduction number",
                          breaks = 0:rt_max_display)
    )  +
    labs(title = paste0("Total — ", scen_labels[scenario]),
         x = "Date") +
    theme_minimal(base_size = BASE_SIZE,
                  base_family = BASE_FAMILY) +
    theme(axis.title.y.right = element_text(color = col_rt,
                                            size = size_axis_title * 0.8,
                                            margin = margin(l = 8)
                                            ),
          plot.margin = margin(8, 15, 8, 12),
          legend.position        = "inside",
          legend.position.inside = c(0.9, 0.9),
          legend.justification   = c(1, 1),
          legend.direction       = "vertical",
          legend.title           = element_blank(),
          legend.text            = element_text(size = 2, margin = margin(l = 0)),
          legend.key.width       = unit(6, "pt"),
          legend.key.height      = unit(3, "pt"),
          legend.key.spacing.x   = unit(0.5, "pt"),
          legend.key.spacing.y   = unit(0.3, "pt"),
          legend.box.margin          = margin(0, 0, 0, 0),
          legend.margin     = margin(1.5, 2, 1.5, 2),
          legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                           colour = "gray10",
                                           linewidth = 0.15),
          legend.key             = element_rect(fill = "transparent", colour = NA)
    )
}

plot_age_CI_scen_from_rds <- function(ci_age_df, obs_age_df, fix_end,
                                      scenario = "Fit",
                                      window = NULL,
                                      show_obs = TRUE) {
  
  age_levels <- c("X19", "X2039", "X4059", "X60over")
  age_lab <- c(
    X19     = "Age 0–19 years",
    X2039   = "Age 20–39 years",
    X4059   = "Age 40–59 years",
    X60over = "Age ≥60 years"
  )
  
  ci <- ci_age_df |>
    dplyr::mutate(
      age = factor(age, levels = seq_along(age_levels),
                   labels = age_levels)
    )
  
  df <- dplyr::left_join(ci, obs_age_df, by = c("t","age")) |>
    dplyr::filter(!is.na(date)) |>
    dplyr::mutate(date = as.Date(date))
  
  vline_x <- fix_date_from_t(obs_age_df, fix_end)
  
  col_ci <- scen_colors_fill[scenario]
  
  g <- ggplot(df, aes(x = date)) +
    geom_ribbon(aes(ymin = low, ymax = high),
                fill = col_ci, alpha = 0.2) +
    geom_line(aes(y = med),
              color = col_ci, linewidth = 0.3)
  
  if (show_obs) {
    g <- g + geom_point(aes(y = observed),
                             shape = 16, color = "black",
                             size = 0.25, stroke = 0, alpha = 0.6)
  }
  if (!is.na(vline_x)) {
    g <- g + geom_vline(xintercept = vline_x,
                        linetype = "dashed", color = "#009e73",
                        linewidth = 0.1)
  }
  
  g +
    facet_wrap(~ age, scales = "fixed", ncol = 2,
               labeller = as_labeller(age_lab)) +
    scale_y_continuous(
      labels = scales::label_number(accuracy = 1, big.mark = ","),
      breaks = function(limits) {
        hi <- limits[2]
        if (!is.finite(hi) || hi <= 0) return(c(0))
        scales::breaks_pretty(n = 3)(c(0, hi))
      }
    ) +
    labs(title = paste0("By Age — ", scen_labels[scenario]),
         x = "Date", y = win_label(window)) +
    theme_minimal(base_size = BASE_SIZE,
                  base_family = BASE_FAMILY)
}

# ============================================================
# 6) 4 panels
# ============================================================
make_supplementary_for_country <- function(country_code,
                                           processed_root ="../processed",
                                           out_root ="../figures"
                                           ) 
  {country_full <- country_name_map[country_code]

  cf2_rate <- country_cf2_rate[country_code]
  scen_labels[["CF2"]] <<- if (!is.na(cf2_rate))
    sprintf("Mask: %.1f%%", cf2_rate) else "Mask: Asia 1.8%, Europe 0.2%"
  message("Generating supplementary: ", country_full)
  
  d <- load_country_summaries(country_code, processed_root)
  fix_end <- d$meta$fix_end
  
  # ---- Panel A: Fit ----
  p_fit_total <- plot_total_CI_scen_from_rds(
    d$ci_total_fit, d$rt_fit, d$obs_tot, fix_end,
    scenario = "Fit"
  ) + theme_lancet_like() + x_tweak_total() + theme_supp() +
    theme(legend.position        = "inside",
          legend.position.inside = c(0.99, 0.999),
          legend.justification   = c(1, 1),
          legend.text            = element_text(size = 5, margin = margin(l = 0)),
          legend.key.width       = unit(8, "pt"),
          legend.key.height      = unit(5, "pt"),
          legend.key.spacing.x   = unit(0.5, "pt"),
          legend.key.spacing.y   = unit(0.3, "pt"),
          legend.margin     = margin(1.5, 2, 1.5, 2),
          legend.box.margin          = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                           colour = "gray10",
                                           linewidth = 0.15),
          legend.key             = element_rect(fill = "transparent", colour = NA))
  
  p_fit_age <- plot_age_CI_scen_from_rds(
    d$ci_age_fit, d$obs_age, fix_end, scenario = "Fit"
  ) + theme_lancet_like() + age_x_tweak() +
    age_strip_thin(strip_lw_pt = 0.5) + theme_supp()
  
  panel_A <- (p_fit_total | p_fit_age) +
    patchwork::plot_layout(widths = c(1, 1.2))
  
  # ---- Panel B: CF2 (less mask) ----
  p_cf2_total <- plot_total_CI_scen_from_rds(
    d$ci_total_cf2, d$rt_cf2, d$obs_tot, fix_end,
    scenario = "CF2"
  ) + theme_lancet_like() + x_tweak_total() + theme_supp() +
    theme(legend.position        = "inside",
          legend.position.inside = c(0.9, 1),
          legend.justification   = c(1, 1),
          legend.text            = element_text(size = 5, margin = margin(l = 0)),
          legend.key.width       = unit(8, "pt"),
          legend.key.height      = unit(5, "pt"),
          legend.key.spacing.x   = unit(0.5, "pt"),
          legend.key.spacing.y   = unit(0.3, "pt"),
          legend.box.margin          = margin(0, 0, 0, 0),
          legend.margin     = margin(1.5, 2, 1.5, 1),
          legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                           colour = "gray10",
                                           linewidth = 0.1),
          legend.key             = element_rect(fill = "transparent", colour = NA))
  
  p_cf2_age <- plot_age_CI_scen_from_rds(
    d$ci_age_cf2, d$obs_age, fix_end, scenario = "CF2"
  ) + theme_submit_like() + age_x_tweak() +
    age_strip_thin(strip_lw_pt = 0.5) + theme_supp()
  
  panel_B <- (p_cf2_total | p_cf2_age) +
    patchwork::plot_layout(widths = c(1, 1.2))
  
  # ---- Panel C: CF1 (more mask) ----
  p_cf1_total <- plot_total_CI_scen_from_rds(
    d$ci_total_cf1, d$rt_cf1, d$obs_tot, fix_end,
    scenario = "CF1"
  ) + theme_lancet_like() + x_tweak_total() + theme_supp() +
    theme(legend.position        = "inside",
          legend.position.inside = c(0.3, 0.999),
          legend.justification   = c(1, 1),
          legend.text            = element_text(size = 5, margin = margin(l = 0)),
          legend.key.width       = unit(8, "pt"),
          legend.key.height      = unit(5, "pt"),
          legend.key.spacing.x   = unit(0.5, "pt"),
          legend.key.spacing.y   = unit(0.3, "pt"),
          legend.margin     = margin(1.5, 2, 1.5, 1),
          legend.box.margin          = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                           colour = "gray10",
                                           linewidth = 0.1),
          legend.key             = element_rect(fill = "transparent", colour = NA))
  
  p_cf1_age <- plot_age_CI_scen_from_rds(
    d$ci_age_cf1, d$obs_age, fix_end, scenario = "CF1"
  ) + theme_submit_like() + age_x_tweak() +
    age_strip_thin(strip_lw_pt = 0.5) + theme_supp()
  
  panel_C <- (p_cf1_total | p_cf1_age) +
    patchwork::plot_layout(widths = c(1, 1.2))
  
  # ---- Panel D: Rt  ----
  rt_all_df <- bind_rows(
    d$rt_fit |> mutate(scenario = "Fit"),
    d$rt_cf1 |> mutate(scenario = "CF1"),
    d$rt_cf2 |> mutate(scenario = "CF2")
  ) |>
    left_join(d$obs_tot |> select(t, date), by = "t") |>
    filter(!is.na(date)) |>
    mutate(scenario = factor(scenario, levels = c( "Fit","CF2", "CF1")))
  
  vline_x <- fix_date_safe(d$obs_tot, fix_end)
  
  panel_D <- ggplot(rt_all_df, aes(x = date, color = scenario, fill = scenario)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.15, color = NA) +
    geom_line(aes(y = med), linewidth = 0.25) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "grey10", linewidth = 0.2) +
    {if (!is.na(vline_x))
      geom_vline(xintercept = vline_x, linetype = "dashed",
                 color = "#36c78c", linewidth = 0.2)
    } +
    scale_color_manual(values = scen_colors_rt, labels = scen_labels,
                       breaks = c("Fit","CF2", "CF1")) +
    scale_fill_manual(values = scen_colors_rt, labels = scen_labels,
                      breaks = c("Fit","CF2", "CF1")) +
    scale_x_date(limits = XLIM_COMMON,
                 breaks = XBREAKS_COMMON,
                 expand = expansion(mult = c(0.02, 0.02)),
                 labels = scales::label_date("%b\n%Y", locale = "en")) +
    coord_cartesian(ylim = c(0, 4)) +
    labs(title = "Effective reproduction number",
         x = "Date", y = "Effective reproduction number") +
    theme_submit_like() + theme_supp() +
    theme(
      legend.position       = "inside",
      legend.position.inside  = c(0.9, 0.9),  
      legend.direction       = "vertical", 
      legend.text            = element_text(size = 5,  
                                            margin =margin(l=2)),      
      legend.key.size        = unit(6, "pt"),                  
      legend.key.width       = unit(8, "pt"),                  
      legend.key.height      = unit(5, "pt"),                  

      legend.title           = element_blank(),
      legend.box.margin      = margin(t = 0.5, r = 0, b = 0.5, l = 0),
      legend.margin     = margin(1.5, 2, 1.5, 1),
      legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                       colour = "gray10",
                                       linewidth = 0.1),
      legend.key             = element_rect(fill = "transparent", colour = NA),
      legend.key.spacing.y   = unit(1.5, "pt"),
      
      axis.title.y   = ggplot2::element_text(size = size_axis_title * 0.85),
     
    )
  
  # ----merge ----
  supplementary_grid <- panel_A / panel_B / panel_C / panel_D +
    patchwork::plot_annotation(
      title      = paste0("Country: ", country_full),
      tag_levels = "A",
      theme = theme(plot.title = element_text(size = size_main_title,
                                              face = "bold",
                                              family = BASE_FAMILY))
    ) +
    patchwork::plot_layout(heights = c(1, 1, 1, 0.8))
  

  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  
  ragg::agg_tiff(
    filename    = file.path(out_root, paste0(country_full, "_supplementary.tiff")),
    width       = 150, height = 190,
    units       = "mm", res = 300, compression = "lzw", scaling = 1
  )
  on.exit(if (dev.cur() != 1) dev.off(), add = TRUE)  
  print(supplementary_grid)
  dev.off()
  
  message("Saved: ", country_full)
  invisible(supplementary_grid)
}

make_supplementary_for_country("JPN") #set contry code

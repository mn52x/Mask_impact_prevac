
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble);library(readr)
  library(ggplot2); library(patchwork); library(cowplot);library(purrr)
  library(posterior); library(stringr); library(ggdist); library(showtext)
})
showtext_auto()


# ── 1. load CFR and population ─────────────────────────────────────────
CFR_path <- "../cfr.csv"
POP_path <- "../population.csv"

cfr_all <- read_csv(CFR_path,  show_col_types = FALSE)
pop_all <- read_csv(POP_path,  show_col_types = FALSE)

get_pop <- function(cc) {
  row <- pop_all |> filter(country == cc)
  if (nrow(row) != 1) stop("multiple rows: ", cc)
  row$total   
}


get_cfr_df <- function(cc) {
  row <- cfr_all |> filter(country == cc)
  if (nrow(row) != 1) stop("mulitiple cfrs: ", cc)
  tibble(
    age = 1:4,
    cfr = as.numeric(c(row$`19`, row$`2039`, row$`4059`, row$`60over`))
  )
}
# ── 2. summarize ────────────────────────────────────────────
#----2_1 make_csv_fit------
.find_header_row <- function(file, max_lines = 2000L){
  con <- file(file, "r"); on.exit(close(con), add = TRUE)
  i <- 0L
  repeat {
    ln <- readLines(con, n = 1L)
    if (length(ln) == 0L) break
    i <- i + 1L
    if (!startsWith(ln, "#") && nzchar(ln)) return(i)
    if (i >= max_lines) break
  }
  stop("cannnot get header: ", file)
}
.count_data_rows <- function(file, header_row) {
  lines <- readLines(file)
  data_lines <- lines[seq(header_row + 1, length(lines))]
  sum(!startsWith(data_lines, "#") & nzchar(trimws(data_lines)))
}

.read_header_names <- function(file, header_row){
  ln <- readLines(file, n = header_row)[header_row]
  trimws(strsplit(ln, ",", fixed = TRUE)[[1]])
}

.normalize_names <- function(nms){
  n1 <- sub("^([A-Za-z]\\w*)\\.(\\d+)\\.(\\d+)$", "\\1[\\2,\\3]", nms, perl = TRUE)
  n2 <- sub("^([A-Za-z]\\w*)\\.(\\d+)$", "\\1[\\2]", n1, perl = TRUE)
  n2
}

make_csv_fit <- function(csvs){
  stopifnot(length(csvs) >= 1, all(file.exists(csvs)))
  
  header_row <- .find_header_row(csvs[1])
  orig_names <- .read_header_names(csvs[1], header_row)   
  norm_names <- .normalize_names(orig_names)              
  
  orig_by_norm <- setNames(orig_names, norm_names)
  
  variables <- function() norm_names  
  
  .read_subset <- function(select_norm, format = c("draws_matrix","draws_df")){
    format <- match.arg(format)
    
    sel_norm <- intersect(select_norm, norm_names)
    if (!length(sel_norm)) stop("not found: ", paste(select_norm, collapse=", "))
    
    sel_orig <- unname(orig_by_norm[sel_norm])
    name2pos <- setNames(seq_along(orig_names), orig_names)
    sel_idx  <- unname(name2pos[sel_orig])
    
    dts <- lapply(seq_along(csvs), function(i){
      h_i    <- .find_header_row(csvs[i])
      n_data <- .count_data_rows(csvs[i], h_i)  
      dt <- data.table::fread(
        csvs[i],
        skip         = h_i,
        header       = FALSE,
        select       = sel_idx,
        nrows        = n_data,     
        showProgress = FALSE
      )
      data.table::setnames(dt, sel_norm)
      dt[, .chain     := i]
      dt[, .iteration := .I]
      dt
    })
    
    big <- data.table::rbindlist(dts, use.names = TRUE, fill = TRUE)
    
    draw_cols <- setdiff(names(big), c(".chain",".iteration"))
    mat <- as.matrix(big[, ..draw_cols])
    colnames(mat) <- draw_cols
    rownames(mat) <- NULL
    
    if (format == "draws_matrix") {
      posterior::as_draws_matrix(mat)
    } else {
      dm <- posterior::as_draws_matrix(mat)
      posterior::as_draws_df(dm) |>
        tibble::add_column(.chain = big$.chain, .iteration = big$.iteration, .before = 1)
    }
  }
  
  draws <- function(variables = NULL, format = c("draws_matrix","draws_df")){
    format <- match.arg(format)
    if (is.null(variables)) stop("set variables")
    
    if (length(variables) == 1 && !grepl("\\[", variables)) {
      pat  <- paste0("^", variables, "(\\[|$)")
      cols_norm <- grep(pat, norm_names, value = TRUE)
    } else {
      cols_norm <- variables
    }
    .read_subset(cols_norm, format = format)
  }
  
  structure(list(csvs = csvs, variables = variables, draws = draws),
            class = "csv_fit")
}
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
      t   = as.integer(stringr::str_match(var, "\\[(\\d+),(\\d+)\\]")[, 2]),
      age = as.integer(stringr::str_match(var, "\\[(\\d+),(\\d+)\\]")[, 3])
    ) |>
    select(.draw, scenario, t, age, mu_cases)
}


calc_diff_summary <- function(fit, rate_label,
                              cum_start, cum_end, pop) {

  M_base <- make_M_long(fit, "M_hat_out",     "base")
  M_cf1  <- make_M_long(fit, "M_hat_CF1_out", "CF1")
  M_cf2  <- make_M_long(fit, "M_hat_CF2_out", "CF2")
  M_long <- bind_rows(M_base, M_cf1, M_cf2)
  rm(M_base, M_cf1, M_cf2); gc()
  
  cum_draws <- M_long %>%
    filter(t >= cum_start, t <= cum_end) %>%
    group_by(.draw, scenario) %>%
    summarise(cases_cum = sum(mu_cases), .groups = "drop") %>%
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

#── 2_2 summary one country────────────────────────────────────────────
calc_country_draws <- function(cc, cum_start, cum_end) {
  
  message("Processing: ", cc)
  pop    <- get_pop(cc)
  cfr_df <- get_cfr_df(cc)
  
  csvs <- dir(
    file.path("../countries", cc),
    pattern    = "2_2_prepare_runstan-.*\\.csv$",
    full.names = TRUE
  )
  fit   <- make_csv_fit(csvs)
  
  M_base <- make_M_long(fit, "M_hat_out",     "base"); gc()
  M_cf1  <- make_M_long(fit, "M_hat_CF1_out", "CF1");  gc()
  M_cf2  <- make_M_long(fit, "M_hat_CF2_out", "CF2");  gc()
  M_long <- bind_rows(M_base, M_cf1, M_cf2)
  rm(M_base, M_cf1, M_cf2); gc()
  
  # cases
  cum_cases <- M_long |>
    filter(t >= cum_start, t <= cum_end) |>
    group_by(.draw, scenario) |>
    summarise(cases_cum = sum(mu_cases), .groups = "drop") |>
    pivot_wider(names_from = scenario, values_from = cases_cum) |>
    transmute(
      .draw,
      averted   = (CF2  - base) / pop * 1e5,
      avertible = (base - CF1)  / pop * 1e5
    ) |>
    pivot_longer(c(averted, avertible),
                 names_to = "scenario", values_to = "value") |>
    mutate(quantity = "cases_per100k")
  
  # deaths
  cum_deaths <- M_long |>
    filter(t >= cum_start, t <= cum_end) |>
    left_join(cfr_df, by = "age") |>
    mutate(deaths = mu_cases * cfr) |>
    group_by(.draw, scenario) |>
    summarise(deaths_cum = sum(deaths), .groups = "drop") |>
    pivot_wider(names_from = scenario, values_from = deaths_cum) |>
    transmute(
      .draw,
      averted   = (CF2  - base) / pop * 1e5,
      avertible = (base - CF1)  / pop * 1e5
    ) |>
    pivot_longer(c(averted, avertible),
                 names_to = "scenario", values_to = "value") |>
    mutate(quantity = "deaths_per100k")
  
  rm(M_long); gc()
  
  bind_rows(cum_cases, cum_deaths) |>
    mutate(country = cc)
}


# ── 3. 16 countries ──────────────────────────────────────────────────────────
countries_16 <- c(
  # row1
  "AUS","HKG","JPN","MYS",
  # row2 (Asia 3 + Europe 1)
  "KOR","THA","TWN","BEL",
  # row3
  "DNK","FRA","DEU","ITA",
  # row4
  "NLD","PRT","ESP","GBR"
)
stopifnot(length(countries_16) == 16)

df_fig3_draws <- map_dfr(
  countries_16,
  \(cc) calc_country_draws(cc, cum_start = 1, cum_end = 999)
)
saveRDS(df_fig3_draws, "df_fig3_draws_SIfix")
gc()

df_fig3_draws　<- readRDS(df_fig3_draws)


# =============================================================================
# Figure : averted/avertible by country
#   upper : cases
#   bottom :deaths
# =============================================================================

asia_oceania <- c("AUS","HKG","JPN","MYS","KOR","THA","TWN")
europe       <- c("BEL","DNK","FRA","DEU","ITA","NLD","PRT","ESP","GBR")
order_fig1   <- c(asia_oceania, europe)

country_labels <- c(
  AUS = "Australia",    HKG = "Hong Kong",     JPN = "Japan",
  KOR = "South Korea",  MYS = "Malaysia",      THA = "Thailand",
  TWN = "Taiwan",       BEL = "Belgium",       DEU = "Germany",
  DNK = "Denmark",      ESP = "Spain",         FRA = "France",
  GBR = "United Kingdom", ITA = "Italy",       NLD = "Netherlands",
  PRT = "Portugal"
)

# --- data frame ---
df_fig3_draws_plot <- df_fig3_draws |>
  mutate(
    country_label = factor(
      country_labels[country],
      levels = rev(country_labels[order_fig1]) 
    ),
    scenario = factor(scenario, levels = c("avertible", "averted")),
    quantity = factor(quantity,
                      levels = c("cases_per100k", "deaths_per100k"),
                      labels = c("Cumulative cases per 100,000 people",
                                 "Cumulative deaths per 100,000 people"))
  )

sep_y <- length(europe) + 0.65   # = 9.5
base_pt <- 11

col_averted   <- "#d55e00"    
col_avertible <- "#0072b2" 


# =============================================================================
# halfeye function
# =============================================================================
make_halfeye_single <- function(df, x_label, fill_color,
                                show_y = TRUE, show_zero = FALSE) {
  
  p <- ggplot(df, aes(x = value, y = country_label)) +
    scale_y_discrete()
  
  p <- p + geom_hline(yintercept = sep_y, linetype = "dashed",
                      color = "#CC79A7", linewidth = 0.35)
  
  if (show_zero) {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed",
                        color = "gray50", linewidth = 0.2)
  }
  
  p <- p +
    stat_halfeye(
      fill           = fill_color,
      color          = fill_color,
      .width         = c(0.50, 0.95),
      point_interval = "median_qi",
      slab_alpha     = 0.35,
      interval_size  = 0.4,
      point_size     = 0.6,
      normalize      = "groups",
      scale          = 0.5
    ) +
    scale_x_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0.01, 0.05))
    ) +
    coord_cartesian(clip = "off") +
    labs(x = x_label, y = NULL) +
    theme_classic(base_family = "Arial") +
    theme(
      axis.text.x  = element_text(size = base_pt * 0.6, color = "black"),
      axis.text.y  = if (show_y) element_text(size = base_pt * 0.7, color = "black",
                                              margin = margin(r = 4))
      else element_blank(),
      axis.title.x = element_text(size = base_pt * 0.75, color = "black",
                                  margin = margin(t = 5)),
      axis.line.x  = element_line(color = "black", linewidth = 0.4),
      axis.line.y  = if (show_y) element_line(color = "black", linewidth = 0.4) else element_blank(),
      axis.ticks.x = element_line(color = "black", linewidth = 0.4),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.margin  = margin(t = 15, r = 8, b = 7, l = 2)
    )
  p
}

# =============================================================================
# panels
# =============================================================================
#  Averted
p_A <- make_halfeye_single(
  df = df_fig3_draws_plot |>
    filter(as.character(quantity) == "Cumulative cases per 100,000 people",
           as.character(scenario) == "averted"),
  x_label = "Cumulative cases averted per 100,000 people",
  fill_color = col_averted, show_y = TRUE
)

p_B <- make_halfeye_single(
  df = df_fig3_draws_plot |>
    filter(as.character(quantity) == "Cumulative deaths per 100,000 people",
           as.character(scenario) == "averted"),
  x_label = "Cumulative deaths averted per 100,000 people",
  fill_color = col_averted, show_y = FALSE
) +
  scale_x_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.05)))

# Avertible
p_C <- make_halfeye_single(
  df = df_fig3_draws_plot |>
    filter(as.character(quantity) == "Cumulative cases per 100,000 people",
           as.character(scenario) == "avertible"),
  x_label = "Cumulative cases avertible per 100,000 people",
  fill_color = col_avertible, show_y = TRUE, show_zero = TRUE
)

p_D <- make_halfeye_single(
  df = df_fig3_draws_plot |>
    filter(as.character(quantity) == "Cumulative deaths per 100,000 people",
           as.character(scenario) == "avertible"),
  x_label = "Cumulative deaths avertible per 100,000 people",
  fill_color = col_avertible, show_y = FALSE, show_zero = TRUE
)

# =============================================================================
# merge
# =============================================================================
rw <- c(1.3, 1)
boundary <- rw[1] / sum(rw)

p_top    <- cowplot::plot_grid(p_A, p_B, ncol = 2, rel_widths = rw, align = "h")
p_bottom <- cowplot::plot_grid(p_C, p_D, ncol = 2, rel_widths = rw, align = "h")

p_fig3 <- cowplot::plot_grid(p_top, p_bottom, nrow = 2, align = "v") +
  theme(plot.margin = margin(0, 0, 0, 0))

p_fig3 <- cowplot::ggdraw(p_fig3) +
  cowplot::draw_plot_label(
    label    = c("A", "B", "C", "D"),
    x        = c(0.11, boundary, 0.11, boundary),
    y        = c(0.99, 0.99, 0.50, 0.50),
    size     = base_pt * 0.9,
    fontface = "bold", family = "Arial",
    hjust = 0, vjust = 1
  )

       dpi = 300, compression = "lzw",bg = "white")

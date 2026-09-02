suppressPackageStartupMessages({
  library(data.table);library(posterior);library(stringr)
  library(dplyr);library(tibble);library(dplyr)
  library(tidyr);library(readr);library(stringr)
})

CFR_path <- "../CFR.csv"
POP_path <- "../population.csv"


#----- run "3_1_read_csv"  to get make_csv_fit function 

calc_cum_from_prefix2d <- function(draws_df, prefix, method = c("sum","last"), by_age = TRUE){
  method <- match.arg(method)
  
  cols <- grep(paste0("^", prefix, "\\["), names(draws_df), value = TRUE)
  if (length(cols) == 0) stop("No time-series cols for prefix: ", prefix)
  
  m <- stringr::str_match(cols, "\\[(\\d+),(\\d+)\\]$")
  tt  <- as.integer(m[,2])
  age <- as.character(as.integer(m[,3]))   
  
  X <- as.matrix(draws_df[, cols, drop = FALSE])
  
  if (by_age) {
    ages <- sort(unique(age))
    out_list <- lapply(ages, function(a){
      idx <- which(age == a)
      idx <- idx[order(tt[idx])]
      Xa  <- X[, idx, drop = FALSE]
      cum <- if (method == "sum") rowSums(Xa, na.rm = TRUE) else Xa[, ncol(Xa)]
      tibble::tibble(draw_id = seq_len(nrow(draws_df)), age = a, cum = cum)
    })
    dplyr::bind_rows(out_list)
  } else {
    t_vals <- sort(unique(tt))
    Xt <- sapply(t_vals, function(ti){
      idx <- which(tt == ti)
      rowSums(X[, idx, drop = FALSE], na.rm = TRUE)
    })
    if (is.vector(Xt)) Xt <- matrix(Xt, ncol = length(t_vals))
    
    cum <- if (method == "sum") rowSums(Xt, na.rm = TRUE) else Xt[, ncol(Xt)]
    tibble::tibble(draw_id = seq_len(nrow(draws_df)), cum = cum)
  }
}
age_map <- c("1"="0019","2"="2039","3"="4059","4"="60over")


cfr_df <- readr::read_csv(CFR_path) %>%
  mutate(country = as.character(country))

get_cfr_vec <- function(cfr_df, country){
  row <- cfr_df %>% dplyr::filter(country == !!country)
  if (nrow(row) != 1) stop("CFR row not unique for: ", country)
  
  as.numeric(c(row$`19`, row$`2039`, row$`4059`, row$`60over`))
}

compute_averted_one_country_2d <- function(fit, country,
                                           prefix_out  = "M_hat_out",
                                           prefix_cf1  = "M_hat_CF1_out",
                                           prefix_cf2  = "M_hat_CF2_out",
                                           cum_method  = c("sum","last"),
                                           keep_by_age = TRUE){
  cum_method <- match.arg(cum_method)
  
  d_out <- fit$draws(variables = prefix_out, format="draws_df") %>% tibble::as_tibble()
  d_c1  <- fit$draws(variables = prefix_cf1, format="draws_df") %>% tibble::as_tibble()
  d_c2  <- fit$draws(variables = prefix_cf2, format="draws_df") %>% tibble::as_tibble()
  
  # ----total（sum age) ----
  cum_out_all <- calc_cum_from_prefix2d(d_out, prefix_out, cum_method, by_age = FALSE) %>% dplyr::rename(cum_out = cum)
  cum_c1_all  <- calc_cum_from_prefix2d(d_c1,  prefix_cf1, cum_method, by_age = FALSE) %>% dplyr::rename(cum_cf1 = cum)
  cum_c2_all  <- calc_cum_from_prefix2d(d_c2,  prefix_cf2, cum_method, by_age = FALSE) %>% dplyr::rename(cum_cf2 = cum)
  
  res_all <- cum_out_all %>%
    dplyr::left_join(cum_c1_all, by="draw_id") %>%
    dplyr::left_join(cum_c2_all, by="draw_id") %>%
    dplyr::mutate(
      country = country,
      age = "all",
      averted   = cum_cf2 - cum_out,
      avertible = cum_out - cum_cf1
    ) %>%
    dplyr::select(country, age, draw_id, cum_out, cum_cf1, cum_cf2, averted, avertible)
  
  if (!keep_by_age) return(res_all)
  
  # ---- by age（age=1..4） ----
  cum_out_age <- calc_cum_from_prefix2d(d_out, prefix_out, cum_method, by_age = TRUE) %>% dplyr::rename(cum_out = cum)
  cum_c1_age  <- calc_cum_from_prefix2d(d_c1,  prefix_cf1, cum_method, by_age = TRUE) %>% dplyr::rename(cum_cf1 = cum)
  cum_c2_age  <- calc_cum_from_prefix2d(d_c2,  prefix_cf2, cum_method, by_age = TRUE) %>% dplyr::rename(cum_cf2 = cum)
  
  res_age <- cum_out_age %>%
    dplyr::left_join(cum_c1_age, by=c("draw_id","age")) %>%
    dplyr::left_join(cum_c2_age, by=c("draw_id","age")) %>%
    dplyr::mutate(
      country = country,
      averted   = cum_cf2 - cum_out,
      avertible = cum_out - cum_cf1
    ) %>%
    dplyr::select(country, age, draw_id, cum_out, cum_cf1, cum_cf2, averted, avertible)
  
  dplyr::bind_rows(res_all, res_age)
}

calc_cum_deaths_from_prefix2d <- function(draws_df, prefix, cfr_vec, method=c("sum","last")){
  method <- match.arg(method)
  
  cols <- grep(paste0("^", prefix, "\\["), names(draws_df), value = TRUE)
  if (!length(cols)) stop("No cols for prefix: ", prefix)
  
  m <- stringr::str_match(cols, "\\[(\\d+),(\\d+)\\]$")
  tt  <- as.integer(m[,2])
  age <- as.integer(m[,3])
  
  X <- as.matrix(draws_df[, cols, drop = FALSE])
  
  deaths_total <- rep(0, nrow(X))
  for (a in sort(unique(age))) {
    idx <- which(age == a)
    idx <- idx[order(tt[idx])]
    Xa  <- X[, idx, drop = FALSE]
    cum_cases_a <- if (method == "sum") rowSums(Xa, na.rm=TRUE) else Xa[, ncol(Xa)]
    deaths_total <- deaths_total + cum_cases_a * cfr_vec[a]
  }
  
  tibble::tibble(draw_id = seq_len(nrow(draws_df)), deaths_cum = deaths_total)
}

compute_averted_deaths_one_country_2d <- function(
    fit, country, cfr_df,
    prefix_out="M_hat_out", prefix_cf1="M_hat_CF1_out", prefix_cf2="M_hat_CF2_out",
    cum_method=c("sum","last")
){
  cum_method <- match.arg(cum_method)
  
  cfr_vec <- get_cfr_vec(cfr_df, country)   # ★コード1と同じ
  
  d_out <- fit$draws(prefix_out, format="draws_df") %>% tibble::as_tibble()
  d_c1  <- fit$draws(prefix_cf1, format="draws_df") %>% tibble::as_tibble()
  d_c2  <- fit$draws(prefix_cf2, format="draws_df") %>% tibble::as_tibble()
  
  out <- calc_cum_deaths_from_prefix2d(d_out, prefix_out, cfr_vec, cum_method) %>%
    dplyr::rename(base = deaths_cum)
  c1  <- calc_cum_deaths_from_prefix2d(d_c1,  prefix_cf1, cfr_vec, cum_method) %>%
    dplyr::rename(CF1  = deaths_cum)
  c2  <- calc_cum_deaths_from_prefix2d(d_c2,  prefix_cf2, cfr_vec, cum_method) %>%
    dplyr::rename(CF2  = deaths_cum)
  out %>%
    dplyr::left_join(c1, by="draw_id") %>%
    dplyr::left_join(c2, by="draw_id") %>%
    dplyr::mutate(
      country = country,
      averted   = CF2 - base,
      avertible = base - CF1
    ) %>%
    select(country, draw_id, base, CF1, CF2, averted, avertible)
}


#====== folder setting ======

root_dir <- "../mcmc_result" # including all country result by country folder with country code
asia_oceania <- c("JPN","KOR","HKG","TWN","THA","MYS","AUS")
region_map <- function(country) ifelse(country %in% asia_oceania, "Asia/Oceania", "Europe")

country_dirs  <- list.dirs(root_dir, full.names=TRUE, recursive=FALSE)
country_codes <- basename(country_dirs)
keep <- stringr::str_detect(country_codes, "^[A-Z]{3}$")
country_dirs  <- country_dirs[keep]
country_codes <- country_codes[keep]

csv_pattern <- "maskage_Poisson_log-.*\\.csv$"
cum_method <- "sum"   

#====== create data ========
all_country_draws_cases <- purrr::map2_dfr(country_codes, country_dirs, function(cc, dd){
  message("Start cases: ", cc)
  csvs <- dir(dd, pattern = csv_pattern, full.names = TRUE)
  csvs <- csvs[!grepl("profile|diagnostic|metric", basename(csvs), ignore.case = TRUE)]
  if (length(csvs) == 0) return(NULL)
  
  fit <- make_csv_fit(csvs)
  
  compute_averted_one_country_2d(
    fit, cc,
    prefix_out = "M_hat_out",
    prefix_cf1 = "M_hat_CF1_out",
    prefix_cf2 = "M_hat_CF2_out",
    cum_method = cum_method,
    keep_by_age = TRUE
  ) %>%
    mutate(region = region_map(cc), outcome = "case")%>%
    transmute(
      country, region, outcome, age, draw_id,
      base = cum_out,
      CF1  = cum_cf1,
      CF2  = cum_cf2,
      averted, avertible
  )
})

saveRDS(all_country_draws_cases, "all_country_draws_cases_SIfix.rds")

cfr_df_wide <- readr::read_csv(CFR_path) %>%
  mutate(country = as.character(country))

all_country_draws_deaths <- purrr::map2_dfr(country_codes, country_dirs, function(cc, dd){
  message("Start deaths: ", cc)
  csvs <- dir(dd, pattern = csv_pattern, full.names = TRUE)
  csvs <- csvs[!grepl("profile|diagnostic|metric", basename(csvs), ignore.case = TRUE)]
  if (length(csvs) == 0) return(NULL)
  
  fit <- make_csv_fit(csvs)
  
  compute_averted_deaths_one_country_2d(
    fit, cc, cfr_df_wide,   # ★ここが重要：cfr_long ではなく wide を渡す
    prefix_out = "M_hat_out",
    prefix_cf1 = "M_hat_CF1_out",
    prefix_cf2 = "M_hat_CF2_out",
    cum_method = cum_method
  ) %>%
    transmute(
      country = cc,
      region  = region_map(cc),
      outcome = "death",
      age = "all",
      draw_id,
      base, CF1, CF2, averted, avertible
    )
})
saveRDS(all_country_draws_deaths, "all_country_draws_deaths_SIfix.rds")

#--summarize
draws_all <- bind_rows(
  all_country_draws_cases %>% filter(age == "all"),
  all_country_draws_deaths %>% filter(age == "all")
) %>%
  mutate(
    outcome = as.character(outcome),
    country = as.character(country),
    region  = as.character(region),
    draw_id = as.integer(draw_id)
  )

# =========
# foramtting
# =========
fmt_ci_int <- function(med, lo, hi){
  med_i <- round(med); lo_i <- round(lo); hi_i <- round(hi)
  paste0(
    format(med_i, big.mark=","),
    " (",
    format(lo_i, big.mark=","),
    "-",  # en dash
    format(hi_i, big.mark=","),
    ")"
  )
}

summarise_ci <- function(df, group_cols, value_col = "value"){
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      median = median(.data[[value_col]], na.rm = TRUE),
      q2.5 = quantile(.data[[value_col]], 0.025, na.rm = TRUE),
      q97.5  = quantile(.data[[value_col]], 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(disp = fmt_ci_int(median, q2.5, q97.5))
}

# =========
# by country table：outcome × country × metric
# =========
country_tbl <- draws_all %>%
  select(outcome, country, region, draw_id, averted, avertible) %>%
  pivot_longer(cols = c(averted, avertible), names_to = "metric", values_to = "value") %>%
  group_by(outcome, country, region, metric, draw_id) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  summarise_ci(group_cols = c("outcome","country","region","metric"), value_col = "value") %>%
  mutate(level = "country")

# =========
# region table ：outcome × region × metric
# =========
region_tbl <- draws_all %>%
  select(outcome, region, draw_id, base, CF1, CF2, averted, avertible) %>%
  group_by(outcome, region, draw_id) %>%
  summarise(
    base = sum(base, na.rm=TRUE),
    CF1  = sum(CF1,  na.rm=TRUE),
    CF2  = sum(CF2,  na.rm=TRUE),
    averted   = sum(averted,   na.rm=TRUE),
    avertible = sum(avertible, na.rm=TRUE),
    .groups="drop"
  ) %>%
  pivot_longer(cols = c(base, CF1, CF2, averted, avertible),
               names_to = "metric", values_to = "value") %>%
  summarise_ci(group_cols = c("outcome","region","metric"), value_col = "value") %>%
  mutate(level="region", country = NA_character_) %>%
  select(level, outcome, region, country, metric, median, q2.5, q97.5, disp)

# =========
# overall table : outcome × region × metric
# =========
overall_tbl <- draws_all %>%
  select(outcome, draw_id, base, CF1, CF2, averted, avertible) %>%
  group_by(outcome, draw_id) %>%                
  summarise(
    base = sum(base, na.rm=TRUE),
    CF1  = sum(CF1,  na.rm=TRUE),
    CF2  = sum(CF2,  na.rm=TRUE),
    averted   = sum(averted,   na.rm=TRUE),
    avertible = sum(avertible, na.rm=TRUE),
    .groups="drop"
  ) %>%
  pivot_longer(cols = c(base, CF1, CF2, averted, avertible),
               names_to = "metric", values_to = "value") %>%
  summarise_ci(group_cols = c("outcome","metric"), value_col = "value") %>%
  mutate(level = "overall", region = "Overall", country = NA_character_) %>%
  select(level, outcome, region, country, metric, median, q2.5, q97.5, disp)

region_tbl_all <- bind_rows(region_tbl, overall_tbl)

# =========
#　per 100k
# =========
#----- country per 100k------
pop_df <- read_csv(POP_path) %>%
  mutate(country = as.character(country)) %>%
  transmute(country, pop_total = as.numeric(total))  

country_tbl_per100k <- draws_all %>%
  left_join(pop_df, by="country") %>%
  select(outcome, country, region, draw_id, pop_total, base, CF1, CF2, averted, avertible) %>%
  pivot_longer(cols = c(base, CF1, CF2, averted, avertible),
               names_to="metric", values_to="value") %>%
  mutate(value = value / pop_total * 1e5) %>%
  group_by(outcome, country, region, metric, draw_id) %>%
  summarise(value = sum(value, na.rm=TRUE), .groups="drop") %>%
  summarise_ci(group_cols = c("outcome","country","region","metric"), value_col="value") %>%
  mutate(level="country_per100k")
# ----- region per 100k ------ 
region_pop <- tibble(country = country_codes) %>%
  mutate(region = region_map(country)) %>%
  left_join(pop_df, by="country") %>%
  group_by(region) %>%
  summarise(pop_total = sum(pop_total, na.rm=TRUE), .groups="drop")

region_tbl_per100k <- draws_all %>%
  group_by(outcome, region, draw_id) %>%
  summarise(
    base = sum(base, na.rm=TRUE),
    CF1  = sum(CF1,  na.rm=TRUE),
    CF2  = sum(CF2,  na.rm=TRUE),
    averted   = sum(averted,   na.rm=TRUE),
    avertible = sum(avertible, na.rm=TRUE),
    .groups="drop"
  ) %>%
  left_join(region_pop, by="region") %>%
  pivot_longer(cols = c(base, CF1, CF2, averted, avertible),
               names_to="metric", values_to="value") %>%
  mutate(value = value / pop_total * 1e5) %>%
  summarise_ci(group_cols = c("outcome","region","metric"), value_col="value") %>%
  mutate(level="region_per100k", country=NA_character_) %>%
  select(level, outcome, region, country, metric, median, q2.5, q97.5, disp)

#------ overall per 100k------
overall_pop <- pop_df %>%
  filter(country %in% country_codes) %>%
  summarise(pop_total = sum(pop_total, na.rm=TRUE)) %>%
  pull(pop_total)

overall_tbl_per100k <- draws_all %>%
  group_by(outcome, draw_id) %>%
  summarise(
    base = sum(base, na.rm=TRUE),
    CF1  = sum(CF1,  na.rm=TRUE),
    CF2  = sum(CF2,  na.rm=TRUE),
    averted   = sum(averted,   na.rm=TRUE),
    avertible = sum(avertible, na.rm=TRUE),
    .groups="drop"
  ) %>%
  pivot_longer(cols = c(base, CF1, CF2, averted, avertible),
               names_to="metric", values_to="value") %>%
  mutate(value = value / overall_pop * 1e5) %>%
  summarise_ci(group_cols = c("outcome","metric"), value_col="value") %>%
  mutate(level="overall_per100k", region="Overall", country=NA_character_) %>%
  select(level, outcome, region, country, metric, median, q2.5, q97.5, disp)


# counts
country_tbl2 <- country_tbl %>%
  select(level, outcome, region, country, metric, median, q2.5, q97.5, disp)

final_tbl <- bind_rows(overall_tbl, region_tbl, country_tbl2) %>%
  arrange(outcome, level, region, country, metric)

# per100k
country_tbl_per100k2 <- country_tbl_per100k %>%
  mutate(country = as.character(country)) %>%
  select(level, outcome, region, country, metric, median, q2.5, q97.5, disp)

final_tbl_per100k <- bind_rows(
  overall_tbl_per100k,
  region_tbl_per100k,
  country_tbl_per100k2
) %>% arrange(outcome, level, region, country, metric)

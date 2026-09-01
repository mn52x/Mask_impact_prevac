
library(tidyverse)
library(lubridate)
library(readr)
library(stringr)
library(zoo)

# === folder name==========================
input_dir <- "../mask-wearing_data"  # including all mask-wearing data in one folder 

# ===load file  ==========================================
files <- list.files(
  input_dir,
  pattern = "^[A-Z]{3}_mask\\.csv$",
  full.names = TRUE
)
stopifnot(length(files) > 0)

# ===target =============================================
date_start <- as.Date("2020-02-01")
date_end   <- as.Date("2021-03-31")
all_days   <- tibble(date = seq(date_start, date_end, by = "day"))

age_cols_target <- c("0_19", "20_39", "40_59", "60_plus", "ALL")

# === function ============================================
process_one <- function(fpath, n_min = 0, clip_eps = 0.0, use_bayes = FALSE, alpha = 2, beta = 2) {
  code <- stringr::str_sub(basename(fpath), 1, 3)
  
  date_start <- as.Date("2020-02-01")
  date_end   <- as.Date("2021-03-31")
  date_seq   <- seq(date_start, date_end, by = "day")
  
  df <- readr::read_csv(fpath, show_col_types = FALSE) %>%
    dplyr::mutate(date = as.Date(date)) %>%
    dplyr::filter(date >= date_start, date <= date_end) %>%
    dplyr::arrange(date)
  
  age_cols_target <- c("0_19","20_39","40_59","60_plus","ALL")
  age_cols <- intersect(names(df), age_cols_target)
  if (length(age_cols) == 0) {
    warning(sprintf("[%s] cannnot find age column: %s", code, fpath)); return(NULL)
  }
  
  n_cols_all <- paste0("N_", age_cols)
  n_cols <- intersect(names(df), n_cols_all)
  
  mask_long <- df %>%
    dplyr::select(date, dplyr::all_of(age_cols)) %>%
    tidyr::pivot_longer(-date, names_to = "age", values_to = "mask")
  
  if (length(n_cols) > 0) {
    n_long <- df %>%
      dplyr::select(date, dplyr::all_of(n_cols)) %>%
      tidyr::pivot_longer(-date, names_to = "ageN", values_to = "N") %>%
      dplyr::mutate(age = stringr::str_remove(ageN, "^N_")) %>%
      dplyr::select(-ageN)
    
    mask_long <- dplyr::left_join(mask_long, n_long, by = c("date","age"))
  } else {
    mask_long <- dplyr::mutate(mask_long, N = NA_integer_)
  }

  mask_long <- mask_long %>%
    dplyr::mutate(
      mask      = as.numeric(mask),
      mask_orig = mask,
      mask      = dplyr::if_else(!is.na(N) & N < n_min, NA_real_, mask)
    )
  
  out_long <- mask_long %>%
    dplyr::group_by(age) %>%
    tidyr::complete(date = date_seq) %>%
    dplyr::arrange(age, date) %>%
    dplyr::mutate(
      mask = zoo::na.approx(mask, x = as.numeric(date), na.rm = FALSE, rule = 2),
      mask = pmin(pmax(mask, 0), 1)
    ) %>%
    dplyr::ungroup()
  
  out_wide_interp <- out_long %>%
    dplyr::select(date, age, mask) %>%
    tidyr::pivot_wider(names_from = age, values_from = mask) %>%
    dplyr::arrange(date)
  
  out_wide_orig <- mask_long %>%        
    dplyr::select(date, age, mask_orig) %>%
    tidyr::pivot_wider(names_from = age, values_from = mask_orig,
                       names_glue = "{age}_orig") %>%
    dplyr::arrange(date)
  
  out_wide <- dplyr::left_join(out_wide_interp, out_wide_orig, by = "date")
  
  out_file <- file.path(out_dir, sprintf("%s_mask_age_interporate.csv", code))
  readr::write_csv(out_wide, out_file)
  
  out_long %>%
    dplyr::mutate(country = code, region = country_to_region(code))
}

# ==plot =============================================
plot_data <- map(files, process_one) %>%
  list_rbind() 

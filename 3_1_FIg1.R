suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr);library(lubridate);library(slider);
  library(ggplot2); library(scales); library(patchwork); library(showtext); library(cowplot);
  library(grid)
})

# =========================
# Settings
# =========================
policy_csv <- "../mask_policydate.csv"
pop_csv    <- "../population_2020_4age.csv"

start_date <- as.Date("2020-02-01")
end_date   <- as.Date("2021-03-31")
who_date   <- as.Date("2020-03-21")

asia_oceania <- c("AUS","JPN","KOR","HKG","TWN","THA","MYS")
europe       <- c("BEL","DEU","DNK","ESP","FRA","GBR","ITA","NLD","PRT")

use_smooth <- FALSE  

# Continent-specific scaling settings
mask_q    <- 1    
headroom  <- 1.15    

# =========================
# Colors 
# =========================
age_cols_base <- c(
  age0019   = "#0072B2",
  age2039   = "#009E73",
  age4059   = "#E69F00",
  age60over = "#D55E00"
)

mix_with_white <- function(col, w = 0.35){
  rgb_mat <- col2rgb(col)
  rgb((rgb_mat[1,]*(1-w) + 255*w)/255,
      (rgb_mat[2,]*(1-w) + 255*w)/255,
      (rgb_mat[3,]*(1-w) + 255*w)/255)
}
age_cols_fill <- vapply(age_cols_base, mix_with_white, character(1), w = 0.05) # bars: lighter
age_cols_line <- vapply(age_cols_base, mix_with_white, character(1), w = 0.3) # lines: slightly lighter

age_levels <- c("age60over","age4059","age2039","age0019") 

age_labels <- c(
  age0019 = "0–19 years" ,
  age2039 = "20–39 years",
  age4059 = "40–59 years",
  age60over  = "60+ years "
)

# =========================
# Population (country, 0019,2039,4059,60,total) -> long
# =========================
pop_long <- read_csv(pop_csv, show_col_types = FALSE) %>%
  rename(country = 1) %>%
  pivot_longer(cols = -country, names_to = "age_raw", values_to = "pop") %>%
  mutate(
    age_group = case_when(
      age_raw == "19"    ~ "age0019",
      age_raw == "2039"  ~ "age2039",
      age_raw == "4059"  ~ "age4059",
      age_raw == "60"    ~ "age60over",
      age_raw == "total" ~ "all",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_group)) %>%
  dplyr::select(country, age_group, pop)

# =========================
# Load files
# =========================
raw_all_df <- read.csv("../DatasetsS1_observed_case_mask.csv")
raw_all_df <- raw_all_df%>%
  mutate(
  continent = case_when(
    country %in% asia_oceania ~ "Asia–Oceania",
    country %in% europe       ~ "Europe",
    TRUE ~ "Other"
  ))

colnames(raw_all_df)
# =========================
# Mask time series (long)
# =========================
mask_ts<- raw_all_df %>%
  {
    mask_cols <- c("mask_age_00_19",  "mask_age_20_39",  "mask_age_40_59",  "mask_age_60over")#"mask_all")
    mask_cols <- mask_cols[mask_cols %in% names(.)]
    select(., country, continent, date, all_of(mask_cols))
  }  %>%
  mutate(date = as.Date(date), country = country) %>%
  pivot_longer(
    cols = starts_with("mask_"),
    names_to = "age_group",
    values_to = "mask"
  ) %>%
  mutate(age_group = recode(age_group,
                            mask_age_00_19  = "age0019",
                            mask_age_20_39  = "age2039",
                            mask_age_40_59  = "age4059",
                            mask_age_60over = "age60over"
                            #mask_all    = "all"
  )) %>%
  mutate(age_group = factor(age_group, levels = age_levels))

# =========================
# Case time series (long) + per100k
# =========================
case_ts<- raw_all_df %>%
  {
    case_cols <- c("case_age_00_19",  "case_age_20_39",  "case_age_40_59",  "case_age_60over")
    case_cols <- case_cols[case_cols %in% names(.)]
    select(., country, continent, date, any_of(case_cols))
  } %>%
  mutate(date = as.Date(date), country = country) %>%
  pivot_longer(
    cols = -c(country, continent, date),
    names_to = "age_raw",
    values_to = "cases"
  ) %>%
  mutate(
    age_raw = sub("^X", "", age_raw),
    age_group = case_when(
      age_raw == "case_age_00_19"  ~ "age0019",
      age_raw == "case_age_20_39"  ~ "age2039",
      age_raw == "case_age_40_59"  ~ "age4059",
      age_raw == "case_age_60over"    ~ "age60over",
      #age_raw == "total" ~ "all",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_group)) %>%
  left_join(pop_long, by = c("country","age_group")) %>%
  mutate(
    inc_per100k = (cases / pop) * 1e5
  ) %>%
  group_by(country, continent, age_group) %>%
  arrange(date) %>%
  mutate(
    # computed but NOT used unless use_smooth == TRUE
    inc7_per100k = slide_dbl(inc_per100k, mean, .before = 6, .complete = FALSE)
  ) %>%
  ungroup()%>%
  mutate(age_group = factor(age_group, levels = age_levels))

# =========================
# Policy dates
# =========================
policy_df <- read_csv(policy_csv, show_col_types = FALSE) %>%
  transmute(country = country,
            policy_date = as.Date(mask_policydate))

# =========================
# Continent-specific dual-axis scales
# =========================
get_ylim_upper <- function(case_ts_df, continent_name, use_smooth, headroom = 1.15){
  total <- case_ts %>%
    filter(continent == continent_name, age_group != "all") %>%
    mutate(y_inc = if (use_smooth) inc7_per100k else inc_per100k) %>%
    group_by(country, date) %>%
    summarise(total_inc = sum(y_inc, na.rm = TRUE), .groups = "drop")
  max(total$total_inc, na.rm = TRUE) * headroom
}

trans_asia <- function(x) log1p(x)

get_ylim_continent <- function(continent_name, use_smooth, headroom = 1.05){
  total <- case_ts %>%
    dplyr::filter(continent == continent_name, age_group %in% age_levels) %>%
    dplyr::mutate(y_inc_raw = if (use_smooth) inc7_per100k else inc_per100k) %>%
    dplyr::group_by(country, date) %>%
    dplyr::summarise(
      total_raw  = sum(y_inc_raw, na.rm = TRUE),
      total_plot = if (continent_name == "Asia–Oceania") {
        log1p(total_raw)   # ★ 合計してからlog1p（旧: Σlog1p）
      } else {
        total_raw
      },
      .groups = "drop"
    )
  max(total$total_plot, na.rm = TRUE) * headroom
}

ylim_asia <- get_ylim_continent("Asia–Oceania", use_smooth, headroom=1.05)
ylim_eur  <- get_ylim_continent("Europe",       use_smooth, headroom=1.05)
ylim_eur  <- 750

# =========================
# Plot function (dual axis, continent-specific scaling)
# =========================
countries_16 <- c(
  # row1
  "AUS","HKG","JPN","MYS",
  # row2 (Asia 3 + Europe 1)
  "KOR","TWN","THA","BEL",
  # row3
  "DNK","FRA","DEU","ITA",
  # row4
  "NLD","PRT","ESP","GBR"
)
stopifnot(length(countries_16) == 16)

country_name <- c(
  AUS="Australia", BEL="Belgium", DEU="Germany", DNK="Denmark",
  ESP="Spain", FRA="France", GBR="United Kingdom", HKG="Hong Kong",
  ITA="Italy", JPN="Japan", KOR="South Korea", MYS="Malaysia",
  NLD="Netherlands", PRT="Portugal", THA="Thailand", TWN="Taiwan"
)

# --- top：Asia / bottom：Europe ---
asia_oceania <- c("AUS","HKG","JPN","MYS","KOR","THA","TWN")
europe       <- c("BEL","DNK","FRA","DEU","ITA","NLD","PRT","ESP","GBR")

# ---- y_lim upper ----
get_ylim_country <- function(country_i, continent_i, use_smooth, headroom = 1.05){
  
  df <- case_ts %>%
    dplyr::filter(country == country_i, continent == continent_i, age_group %in% age_levels) %>%
    dplyr::mutate(y_inc_raw = if (use_smooth) inc7_per100k else inc_per100k) %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(
      total_plot = if (continent_i == "Asia–Oceania") {
        sum(trans_asia(y_inc_raw), na.rm = TRUE)
      } else {
        sum(y_inc_raw, na.rm = TRUE)
      },
      .groups = "drop"
    )
  
  max(df$total_plot, na.rm = TRUE) * headroom
}

# ---- common theme ----
base_pt <- 7

theme_fig <- theme(
  text = element_text(family="Arial"),
  panel.grid = element_blank(),
  panel.border = element_rect(color="black", fill=NA, linewidth=0.4),
  axis.ticks.length = unit(0.08, "cm"),

  axis.text.x = element_text(size = base_pt*0.65, margin = margin(t = 1), lineheight = 0.75),
  axis.text.y = element_text(size = base_pt*0.7, margin = margin(r = 1)),
  axis.text.y.right = element_text(size = base_pt*0.0, margin = margin(l = 1)),
  axis.ticks  = element_line(color = "black", linewidth = 0.3),
  axis.text   = element_text(color = "black"),
  axis.title  = element_text(color = "black"),
  panel.spacing.x = unit(0.4, "lines"),
  panel.spacing.y = unit(0.3, "lines"),
  plot.margin = margin(t=0, r=6, b=4, l=7),
  legend.position = "none",
  plot.title = element_text(size=base_pt*0.7, hjust=0.3,
                            margin = margin(t = 0, b = -2))
)
# =========================
# One country fuction （Asia:log10(y+1)）
# =========================
add_linear_right_axis <- function(p_main,
                                  mask_max   = 100,
                                  expand_fac = expansion(mult = c(0.02, 0.02)),
                                  base_pt_   = base_pt) {
  
  p_axis <- ggplot(data.frame(x = 1, y = c(0, mask_max)), aes(x, y)) +
    geom_blank() +
    scale_y_continuous(
      limits   = c(0, mask_max),
      breaks   = seq(0, 100, 25),
      labels   = paste0(seq(0, 100, 25), "%"),
      position = "right",
      expand   = expand_fac
    ) +
    theme_void(base_family = "Arial") +
    theme(
      # ★ 0.80 → 0.6（他パネル axis.text.y と統一）
      axis.text.y.right      = element_text(size  = base_pt_ * 0.6,
                                            hjust = 0,
                                            margin = margin(l = 1)),  # ★ 2 → 1
      axis.ticks.y.right     = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length.y.right = unit(0.08, "cm"),                    # ★ 0.12 → 0.1
      panel.background = element_rect(fill = NA, color = NA),
      plot.background  = element_rect(fill = NA, color = NA)
    )
  
  p_main + patchwork::inset_element(
    p_axis,
    left   = 0,
    right  = 1.21,    
    bottom = 0,
    top    = 1,
    align_to = "panel"
  )
}

#-----BEL------ 

add_left_axis_inset <- function(p_main,
                                ylim_upper,
                                base_pt_ = base_pt) {
  
  y_breaks <- pretty(c(0, ylim_upper), n = 5)
  y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= ylim_upper]
  
  p_axis <- ggplot(data.frame(x = 1, y = c(0, ylim_upper)), aes(x, y)) +
    geom_blank() +
    scale_y_continuous(
      limits   = c(0, ylim_upper),
      breaks   = y_breaks,
      labels   = scales::label_comma()(y_breaks),
      position = "left",
      expand   = expansion(mult = c(0.02, 0.02))
    ) +
    theme_void(base_family = "Arial") +
    theme(
      axis.text.y.left      = element_text(size   = base_pt_ * 0.6,
                                           hjust  = 1,
                                           margin = margin(r = 1)),  
      axis.ticks.y.left     = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length.y.left = unit(0.08, "cm"),                     
      panel.background = element_rect(fill = NA, color = NA),
      plot.background  = element_rect(fill = NA, color = NA)
    )
  
  p_main + patchwork::inset_element(
    p_axis,
    left     = -0.157,  
    right    = 0,
    bottom   = 0,
    top      = 1,
    align_to = "panel"
  )
}

##-----plot-----#####
make_dual_country <- function(country_i, continent_i, ylim_upper, use_smooth,
                              show_right = FALSE, show_legend = FALSE){
  
  is_asia <- (continent_i == "Asia–Oceania")
  force_left <- (country_i == "BEL")   
  
  # case
  case_plot <- case_ts %>%
    dplyr::filter(country == country_i, continent == continent_i,
                  age_group %in% age_levels) %>%
    dplyr::mutate(
      y_inc_raw = if (use_smooth) inc7_per100k else inc_per100k,
      age_group = factor(age_group, levels = age_levels)
    ) %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(
      total_raw = sum(y_inc_raw, na.rm = TRUE),
      y_inc = if (is_asia) {
        log1p(total_raw) * (y_inc_raw / total_raw)  
      } else {
        y_inc_raw
      }
    ) %>%
    dplyr::ungroup()
  
  # mask
  mask_plot <- mask_ts %>%
    dplyr::filter(country == country_i, continent == continent_i, age_group %in% age_levels) %>%
    dplyr::mutate(
      mask_pct = 100 * mask,
      y_mask_on_inc = mask_pct * (ylim_upper / 100),  
      age_group = factor(age_group, levels = age_levels)
    )
  
  policy_sub <- policy_df %>% dplyr::filter(country == country_i)
  
  lim_age <- rev(age_levels)
  lab_age <- unname(age_labels[lim_age])
  
  sc_fill <- if (show_legend) {
    scale_fill_manual(
      values = age_cols_fill,
      limits = lim_age, breaks = lim_age, labels = lab_age,
      drop = FALSE, na.translate = FALSE,
      name = "Age group",
      guide = guide_legend(reverse = TRUE, nrow = 1)
    )
  } else {
    scale_fill_manual(values = age_cols_fill, limits = lim_age,
                      drop = FALSE, na.translate = FALSE, guide = "none")
  }
  
  sc_col <- scale_color_manual(values = age_cols_line, limits = lim_age,
                               drop = FALSE, na.translate = FALSE, guide = "none")
  
  secax_common <- if (show_right && !is_asia) {     
         sec_axis(
           trans  = ~ . * (100 / ylim_upper),
          name   = NULL,
         　breaks = seq(0, 100, 25),
           labels = function(x) paste0(x, "%")
         )
       } else {
         waiver()
       }
  
  y_scale <- if (is_asia) {
    scale_y_continuous(
      limits = c(0, ylim_upper),
      expand = expansion(mult = c(0.02, 0.02)),
      breaks = log1p(c(0, 1, 5, 10, 20, 50)),
      labels =       c(0, 1, 5, 10, 20, 50),
      sec.axis = secax_common
    )
  } else {
    scale_y_continuous(
      limits   = c(0, ylim_upper),
      expand   = expansion(mult = c(0.02, 0.02)),
      sec.axis = secax_common
    )
  }
  ggplot() +
    geom_col(data = case_plot, aes(x=date, y=y_inc, fill=age_group), width=1, alpha=0.75) +
    geom_line(data = mask_plot, aes(x=date, y=y_mask_on_inc, color=age_group, group=age_group),
              linewidth=0.5, alpha=0.99, na.rm=TRUE) +
    geom_vline(xintercept = who_date, linetype="dashed", color="black", linewidth=0.3) +
    geom_vline(data = policy_sub, aes(xintercept=policy_date),
               inherit.aes=FALSE, linetype="dashed", color="red", linewidth=0.3) +
    sc_fill +
    sc_col +
    y_scale+
    scale_x_date(
      name=NULL, limits=c(start_date, end_date),
      breaks=seq(as.Date("2020-02-01"), as.Date("2021-04-01"), by="4 months"),
      labels=scales::date_format("%b\n%Y")
    ) +
    labs(title = unname(country_name[country_i])) +   
    theme_bw(base_size = base_pt) +
    theme_fig +
    theme(
  
      axis.ticks.y.right = element_line(color="black", linewidth=0.3),
      axis.text.y.right  = if (show_right) element_text(size=base_pt*0.6) else element_blank(),
      
      axis.text.y  = element_text(size = base_pt * 0.6),
      axis.ticks.y = element_line(color = "black", linewidth = 0.3),
      
      axis.title.y       = element_blank(),
      axis.title.y.right = element_blank(),
      
      plot.title = element_text(
        #face = "bold",
        hjust = 0.5,
        vjust = -0.01,                    
        margin = margin(t = -0.01, b = 2),
        size = base_pt * 0.9
      )
    )
}
# =========================
# make grid
# =========================

make_one <- function(cc, i) {
  
  cont  <- if (cc %in% asia_oceania) "Asia–Oceania" else "Europe"
  ylim  <- if (cont == "Asia–Oceania") ylim_asia else ylim_eur
  
  col   <- ((i - 1) %% 4) + 1
  row   <- ((i - 1) %/% 4) + 1
  
  show_right  <- if (cont == "Asia–Oceania") FALSE else col == 4
  show_x      <- (row == 4)
  show_legend <- (i == 1)
  
  p <- make_dual_country(
    country_i   = cc,
    continent_i = cont,
    ylim_upper  = ylim,
    use_smooth  = use_smooth,
    show_right  = show_right,
    show_legend = show_legend
  )
  
  if (!show_x) {
    p <- p + theme(axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank())
  }
  
  if (col != 1) {
    p <- p + theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  # ---- MYS----
  if (cc == "MYS") {
    p <- add_linear_right_axis(p, mask_max = 100, base_pt_ = base_pt)
  }
  
  # ---- BEL ----
  if (cc == "BEL") {
    p <- add_left_axis_inset(p, ylim_upper = ylim_eur, base_pt_ = base_pt)
  }
  
  p
}

make_grid <- function(codes, continent_label, ncol_grid = 4) {
  
  make_one_inner <- function(cc, i) {
    
    cont  <- if (cc %in% asia_oceania) "Asia–Oceania" else "Europe"
    ylim  <- if (cont == "Asia–Oceania") ylim_asia else ylim_eur
    
    col   <- ((i - 1) %% 4) + 1
    row   <- ((i - 1) %/% 4) + 1
    
    show_right  <- if (cont == "Asia–Oceania") FALSE else col == 4
    show_x      <- (row == 4)
    show_legend <- (i == 1)
    
    p <- make_dual_country(
      country_i   = cc,
      continent_i = cont,
      ylim_upper  = ylim,
      use_smooth  = use_smooth,
      show_right  = show_right,
      show_legend = show_legend
    )
    
    if (!show_x) {
      p <- p + theme(axis.text.x  = element_blank(),
                     axis.ticks.x = element_blank())
    }
    
    if (col != 1) {
      p <- p + theme(axis.text.y  = element_blank(),
                     axis.ticks.y = element_blank())
    }
    
    if (cc == "KOR") {
      p <- add_linear_right_axis(p, mask_max = 100, base_pt_ = base_pt)
    }
    
    if (cc == "BEL") {
      p <- add_left_axis_inset(p, ylim_upper = ylim_eur, base_pt_ = base_pt)
    }
    
    p
  }
  
  plist <- Map(make_one_inner, codes, seq_along(codes))
  wrap_plots(plist, ncol = ncol_grid, byrow = TRUE)
}

p_asia <- make_grid(asia_oceania, "Asia–Oceania", ncol_grid=4)
p_eur  <- make_grid(europe,       "Europe",       ncol_grid=4)

plist  <- Map(make_one, countries_16, seq_along(countries_16))
p_grid <- patchwork::wrap_plots(plist, ncol = 4, byrow = TRUE)

p_core <- (p_asia / p_eur) + plot_layout(heights=c(1,1))

p_spacer <- patchwork::plot_spacer()

p_out <- p_grid / p_spacer +
  patchwork::plot_layout(heights = c(1.1, 0.01))

# =========================
# line between Asia and Europe
# =========================
h_leg <- 0.11  
lw <- 0.4
x_cut <- 0.705    
y_mid <- 0.548　
y_top <- 0.766

# =========================
# merge all 
# =========================
p_grid_titled <- p_grid +
  plot_annotation(
    caption = NULL,
    theme = theme(plot.margin = margin(0, 0, 0, 0))
  )
x_label_y <- 0.03  

leg_y       <- 0.020              
leg_box_h   <- 0.022              
leg_box_w   <- 0.014              
leg_text_sz <- base_pt * 0.63     
leg_gap     <- 0.005              

p_final <- cowplot::ggdraw(p_out,
                           xlim = c(-0.01,1.01),   
                           ylim = c(-0.02, 1.02)    
                           ) +  
 
cowplot::draw_line(x = c(0.02, x_cut),   y = c(y_mid, y_mid),
                   linewidth = lw, linetype = "F1", color = "#CC79A7") +　
cowplot::draw_line(x = c(x_cut, 0.965),  y = c(y_top, y_top),
                     linewidth = lw, linetype = "F1", color = "#CC79A7") +
cowplot::draw_line(x = c(x_cut, x_cut),  y = c(y_mid, y_top),
                     linewidth = lw, linetype = "F1", color ="#CC79A7") +

cowplot::draw_label("A",
                    x = 0.03, y = y_top + 0.22,
                    hjust = 0, vjust = 0,
                    size = base_pt , fontface = "bold",
                    fontfamily = "Arial") +
                  
cowplot::draw_label("(Log scale)",
                      x = 0.05, y = y_top + 0.22,  
                      hjust = 0, size = base_pt * 0.6,
                      color = "black", fontfamily = "Arial")+
  
cowplot::draw_label("B",
                      x = 0.03, y = y_mid - 0.023,
                      hjust = 0, vjust = 0,
                      size = base_pt, fontface = "bold",
                      fontfamily = "Arial") +

cowplot::draw_label("Reported cases per 100,000 people",
                    x     = 0.0001,   
                    y     = 0.52,    
                    angle = 90,
                    hjust = 0.5, vjust = 1,
                    size  = base_pt,
                    fontfamily = "Arial") +
 
cowplot::draw_label("Mask-wearing coverage (%)",
                    x     = 0.996,   
                    y     = 0.52,
                    angle = 270,
                    hjust = 0.5, vjust = 1,
                    size  = base_pt ,
                    fontfamily = "Arial") +
  
cowplot::draw_label("Date",
                    x     = 0.5,
                    y     = x_label_y + 0.015,   
                    hjust = 0.5, vjust = 0,
                    size  = base_pt*0.9,
                    fontfamily = "Arial") +

cowplot::draw_grob(
  grid::rectGrob(gp = grid::gpar(fill = as.vector(age_cols_fill["age0019"]),
                                 col = NA, alpha = 0.75)),
  x = 0.030, y = leg_y, width = leg_box_w, height = leg_box_h,
  hjust = 0, vjust = 0.5) +
                  
cowplot::draw_label("0–19 years",
                      x = 0.030 + leg_box_w + leg_gap, y = leg_y,
                      hjust = 0, vjust = 0.5,
                      size = leg_text_sz, fontfamily = "Arial") +
  
cowplot::draw_grob(
    grid::rectGrob(gp = grid::gpar(fill = as.vector(age_cols_fill["age2039"]),
                                   col = NA, alpha = 0.75)),
    x = 0.140, y = leg_y, width = leg_box_w, height = leg_box_h,
    hjust = 0, vjust = 0.5) +
                  
cowplot::draw_label("20–39 years",
                      x = 0.140 + leg_box_w + leg_gap, y = leg_y,
                      hjust = 0, vjust = 0.5,
                      size = leg_text_sz, fontfamily = "Arial") +

cowplot::draw_grob(
    grid::rectGrob(gp = grid::gpar(fill = as.vector(age_cols_fill["age4059"]),
                                   col = NA, alpha = 0.75)),
    x = 0.255, y = leg_y, width = leg_box_w, height = leg_box_h,
    hjust = 0, vjust = 0.5) +
                  
cowplot::draw_label("40–59 years",
                      x = 0.255 + leg_box_w + leg_gap, y = leg_y,
                      hjust = 0, vjust = 0.5,
                      size = leg_text_sz, fontfamily = "Arial") +
  
cowplot::draw_grob(
    grid::rectGrob(gp = grid::gpar(fill = as.vector(age_cols_fill["age60over"]),
                                   col = NA, alpha = 0.75)),
    x = 0.370, y = leg_y, width = leg_box_w, height = leg_box_h,
    hjust = 0, vjust = 0.5) +
                  
cowplot::draw_label("60+ years",
                      x = 0.370 + leg_box_w + leg_gap, y = leg_y,
                      hjust = 0, vjust = 0.5,
                      size = leg_text_sz, fontfamily = "Arial") + 
                  
cowplot::draw_line(x = c(0.475, 0.495), y = c(leg_y, leg_y),
                   linetype = "solid", color = "red", linewidth = 0.8) +
                  
cowplot::draw_label("National mask-wearing policy",
                      x = 0.5, y = leg_y,
                      hjust = 0, vjust = 0.5,
                      size = leg_text_sz, fontfamily = "Arial") +
                  
cowplot::draw_line(x = c(0.730, 0.750), y = c(leg_y, leg_y),
                     linetype = "solid", color = "black", linewidth = 0.8) +
                  
cowplot::draw_label("WHO mask-wearing recommendation",
                      x = 0.755, y = leg_y,
                      hjust = 0, vjust = 0.5,
                      size = leg_text_sz, fontfamily = "Arial")


p_final

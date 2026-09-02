suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(cowplot)
  library(posterior); library(stringr); library(showtext)
})
showtext_auto()

base_pt      <- 8
country_name <- "JPN"


scalar_var        <- "Ct"   
susceptibility_var <- "rel_sus"    
infectivity_var    <- "rel_inf"  

age_3 <- c("1", "2", "3")
age_labels_3 <- c(
  "1" = "0–19 years old",
  "2" = "20–39 years old",
  "3" = "40–59 years old"
)

col_scalar <- "#0072B2"        
col_susc   <- "#D85A30"         
col_infec  <- "#009E73"         

intervention_colors <- c(
  "SOE_1" = "#E69F00",
  "SOE_2" = "#56B4E9",
  "schoolclosure" = "#CC79A7"
)

intervention_labels <- c(
  "SOE_1" = "State of Emergency 1",
  "SOE_2" = "State of Emergency 1",
  "schoolclosure" =  "School closure"
)
intervention_df <- read.csv("../JPN_policy.csv")

# =================
# load fit 
# 1-1 get standata from 2_2_prepare_runstan
# 1-2 run "3_1_read_csv.R" 
#=====================  
csvs <- dir("../JPN", #folder including stan result csv files
            pattern = "2_1_prepare_runstan-XXX-.*\\.csv$", #XXX reflects saved time
            full.names = TRUE)
fit_base <- make_csv_fit(csvs)

date_df <- tibble(
  t    = seq_len(nrow(data)),
  date = as.Date(data$date)
)

# -----------------------------------------------------------------------------
# 2. posterior
# -----------------------------------------------------------------------------

# ---- scalar ----
make_scalar_summary <- function(fit, varname) {
  fit$draws(varname, format = "draws_matrix") %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "var", values_to = "value") %>%
    mutate(t = as.integer(str_match(var, "\\[(\\d+)\\]")[, 2])) %>%
    group_by(t) %>%
    summarise(
      med  = median(value, na.rm = TRUE),
      low  = quantile(value, 0.025, na.rm = TRUE),
      high = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(date_df, by = "t")
}

# ---- 2b.t × age----
step_df <- tibble(
  t       = seq_len(stan_data$T),
  step_id = stan_data$step_id,   
  date    = as.Date(data$date)
)

make_age_summary <- function(fit, varname) {
  fit$draws(varname, format = "draws_matrix") %>%
    as.data.frame() %>%
    mutate(.draw = row_number()) %>%
    pivot_longer(-.draw, names_to = "var", values_to = "value") %>%
    mutate(
      step_id = as.integer(str_match(var, "\\[(\\d+),(\\d+)\\]")[, 2]),
      age     = as.character(str_match(var, "\\[(\\d+),(\\d+)\\]")[, 3])
    ) %>%
    filter(age %in% age_3) %>%
    group_by(step_id, age) %>%
    summarise(
      med  = median(value, na.rm = TRUE),
      low  = quantile(value, 0.025, na.rm = TRUE),
      high = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(step_df, by = "step_id") %>%
    mutate(
      age       = factor(age, levels = age_3),
      age_label = age_labels_3[as.character(age)]
    ) %>%
    arrange(age, t)
}

df_scalar <- make_scalar_summary(fit_base, scalar_var)
df_susc   <- make_age_summary(fit_base, susceptibility_var)
df_infec  <- make_age_summary(fit_base, infectivity_var)

# -----------------------------------------------------------------------------
# 3. visualise
# -----------------------------------------------------------------------------
theme_fig4 <- theme_classic(base_family = "Arial") +
  theme(
    axis.text    = element_text(size = base_pt * 0.6, color = "black"),
    axis.title.y = element_text(size = base_pt * 0.6, color = "black",
                                margin = margin(r = 4)),
    axis.title.x = element_blank(),
    axis.line    = element_line(color = "black", linewidth = 0.4),
    axis.ticks   = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "none",
    plot.margin  = margin(t = 15, r = 4, b = 10, l = 4)
  )

x_scale <- scale_x_date(
  limits      = c(min(date_df$date), max(date_df$date)),
  date_breaks = "4 months",
  labels      = scales::date_format("%b\n%Y")
)


add_interventions <- function(p) {
  for (i in seq_len(nrow(intervention_df))) {
    p <- p + annotate(
      "rect",
      xmin  = as.Date(intervention_df$start_date[i]),  
      xmax  = as.Date(intervention_df$end_date[i]),    
      ymin  = -Inf, ymax = Inf,
      fill  = intervention_colors[intervention_df$policy[i]],
      alpha = 0.15
    )
  }
  p
}

p_A_base <- ggplot(df_scalar, aes(x = date)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "gray50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = low, ymax = high),
              fill = col_scalar, alpha = 0.20) +
  geom_line(aes(y = med),
            color = col_scalar, linewidth = 0.3) +
  x_scale +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(title = "Scalar parameter",
    y = "Scalar parameter") +
  theme_fig4 +
  theme (
    plot.title = element_text(
    size = base_pt * 0.8,
    hjust = 0, family = "Arial",
    margin = margin(b = 15))
  )
  

p_A <- add_interventions(p_A_base)
p_A

make_age_panel <- function(df, col, y_label, shared_ylim = NULL) {
  age_list <- levels(df$age)

  plist <- lapply(age_list, function(a) {
    df_a <- df %>% filter(age == a)
    label <- age_labels_3[a]

    p <- ggplot(df_a, aes(x = date)) +
      geom_hline(yintercept = 1, linetype = "dashed",
                 color = "gray50", linewidth = 0.3) +
      geom_ribbon(aes(ymin = low, ymax = high),
                  fill = col, alpha = 0.20) +
      geom_line(aes(y = med),
                color = col, linewidth = 0.3) +
      x_scale +
      labs(y = label) +
      theme_fig4 +
      theme(
        axis.text.x  = if (a == "3") element_text(size = base_pt * 0.6,
                                                  color = "black", angle = 0,
                                                  hjust = 1, vjust = 0.5,
                                                  margin = margin(r = 4)) else element_blank(),
        axis.title.y = element_text(size = base_pt * 0.6, color = "black",
                                     angle = 90, hjust = 0.5,
                                     margin = margin(r = 4)))

    if (!is.null(shared_ylim)) {
      p <- p + scale_y_continuous(
        limits = shared_ylim,
        expand = expansion(mult = c(0.02, 0.05))
      )
    } else {
      p <- p + scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.05))
      )
    }
    add_interventions(p)
  })


  wrap_plots(plist, ncol = 1) +
    plot_annotation(title = y_label,
                    theme = theme(
                      plot.title = element_text(
                        size = base_pt * 0.8,
                        hjust = 0.3, family = "Arial",
                        margin = margin(b = 0)
                      )
                    ))
}


ylim_susc  <- range(c(df_susc$low,  df_susc$high),  na.rm = TRUE)
ylim_infec <- range(c(df_infec$low, df_infec$high), na.rm = TRUE)
ylim_BC    <- range(c(ylim_susc, ylim_infec))  

p_B <- make_age_panel(df_susc,  col_susc,  "Relative Susceptibility")
p_B
p_C <- make_age_panel(df_infec, col_infec, "Relative Infectiousness") 
p_C

p_BC <- cowplot::plot_grid(
  p_B, p_C,
  ncol       = 2,
  rel_heights = c(1.1, 1),
  labels     = c("B", "C"),
  label_size       = base_pt * 1.1,
  label_fontface   = "bold",
  label_fontfamily = "Arial",
  label_x = c(0.1, 0.1),
  label_y = c(0.99, 0.99)
)
p_BC

legend_intervention_df <- data.frame(
  policy = names(intervention_colors),
  label  = unname(intervention_labels),
  color  = unname(intervention_colors)
)

p_leg_intervention <- ggplot(legend_intervention_df,
                             aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                                 fill = policy)) +
  geom_rect(alpha = 0.4) +
  scale_fill_manual(
    values = intervention_colors,
    labels = intervention_labels,
    name   = "Intervention"
  ) +
  guides(fill = guide_legend(
    title.position = "top",
    override.aes   = list(alpha = 0.4)
  )) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = base_pt * 0.6, face = "bold"),
    legend.text      = element_text(size = base_pt * 0.5),
    legend.key.size  = unit(0.6, "lines"),
    legend.key.width = unit(1.0, "lines"),
    legend.spacing.y = unit(3, "pt")
  )


param_levels <- c("Scalar parameter", "Relative Susceptibility", "Relative Infectiousness") 

param_legend_df <- data.frame(
  param = factor(param_levels, levels = param_levels),   
  color = c(col_scalar, col_susc, col_infec)
)

p_leg_param <- ggplot(param_legend_df,
                      aes(x = 1, y = param, color = param)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = setNames(param_legend_df$color, param_legend_df$param),
    labels = c(Scalar         = "Scalar parameter",
               Susceptibility = "Relative susceptibility",
               Infectiousness = "Relative infectiousness"),
    breaks = param_levels,      
    name   = "Parameter"
  ) + 
  guides(color = guide_legend(
    title.position = "top",
    override.aes   = list(linewidth = 0.5)
  )) +
  theme_void(base_family = "Arial") +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = base_pt * 0.52, face = "bold"),
    legend.text      = element_text(size = base_pt * 0.48),
    legend.key.size  = unit(0.5, "lines"),
    legend.key.width = unit(1.0, "lines"),
    legend.spacing.y = unit(3, "pt")
  )

leg_int   <- cowplot::get_legend(p_leg_intervention)
leg_param <- cowplot::get_legend(p_leg_param)
p_legend  <- cowplot::plot_grid(
  leg_int, leg_param,
  ncol    = 2,
  rel_widths = c(1, 1)
)

p_fig4 <- cowplot::plot_grid(
  p_A, p_BC, p_legend,
  ncol       = 1,
  rel_heights = c(0.8, 1.6, 0.3),
  labels     = c("A", "", ""),
  label_size       = base_pt * 1.1,
  label_fontface   = "bold",
  label_fontfamily = "Arial",
  label_x = c(0.02, 0.02, 0.02),
  label_y = c(0.92, 1.00, 1.00)
)
p_fig4

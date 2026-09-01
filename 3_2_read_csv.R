suppressPackageStartupMessages({
  library(data.table)
  library(posterior)
  library(stringr)
})

# find first row
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
  stop("ヘッダ行を特定できません: ", file)
}

# get header names
.read_header_names <- function(file, header_row){
  ln <- readLines(file, n = header_row)[header_row]
  trimws(strsplit(ln, ",", fixed = TRUE)[[1]])
}

.normalize_names <- function(nms){
  # 2 dimention
  n1 <- sub("^([A-Za-z]\\w*)\\.(\\d+)\\.(\\d+)$", "\\1[\\2,\\3]", nms, perl = TRUE)
  # 1 dimention）
  n2 <- sub("^([A-Za-z]\\w*)\\.(\\d+)$", "\\1[\\2]", n1, perl = TRUE)
  n2
}

make_csv_fit <- function(csvs){
  stopifnot(length(csvs) >= 1, all(file.exists(csvs)))
  
  # get header from first csv file
  header_row <- .find_header_row(csvs[1])
  orig_names <- .read_header_names(csvs[1], header_row)   
  norm_names <- .normalize_names(orig_names)              
  
  # normalize name
  orig_by_norm <- setNames(orig_names, norm_names)
  
  variables <- function() norm_names 
  
  .read_subset <- function(select_norm, format = c("draws_matrix","draws_df")){
    format <- match.arg(format)
    
    sel_norm <- intersect(select_norm, norm_names)
    if (!length(sel_norm)) stop("cannot find: ", paste(select_norm, collapse=", "))
    
    sel_orig <- unname(orig_by_norm[sel_norm])
    name2pos <- setNames(seq_along(orig_names), orig_names)
    sel_idx  <- unname(name2pos[sel_orig])
  
    dts <- lapply(seq_along(csvs), function(i){
      h_i <- .find_header_row(csvs[i])
      dt  <- data.table::fread(csvs[i],
                               skip = h_i, header = FALSE,
                               select = sel_idx,
                               showProgress = FALSE)
      data.table::setnames(dt, sel_norm) 
      dt[, .chain := i]
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

## create fit 
csvs <- dir("../country", #folder including stan result csv files
            pattern = "2_1_prepare_runstan-XXX-.*\\.csv$", #XXX reflects saved time
            full.names = TRUE)
fit <- make_csv_fit(csvs)

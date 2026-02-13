###############################################################################
# µEHT ANALYSIS PIPELINE — STEP 1
# -----------------------------------------------------------------------------
#
# PURPOSE:
#   Interactive GUI tool to process engineered heart tissue (EHT) outputs.
#
#   The script:
#     1. Recursively scans an input directory for "Summary_Data.csv" files
#     2. Reads files robustly (automatic delimiter + header detection)
#     3. Adds experimental metadata (Batch, Donor, Week_Time)
#     4. Extracts pacing step from file names (step_001–step_006)
#     5. Computes expected pacing frequency
#     6. Calculates deviation and relative error
#     7. Determines whether each recording is "Following"
#        based on a user-defined tolerance
#     8. Exports a timestamped combined dataset
#
# OUTPUT:
#   A timestamped folder containing:
#     Combined_Summary_with_Following_<tag>_tolXXpct.csv
#
# NOTES:
#   - Uses ttk widgets for macOS dark-mode compatibility.
#   - Designed for interactive use (folder pickers + GUI).
#   - Intended as Step 1 of a 3-step EHT processing pipeline.
#
# AUTHOR:
#   Michele Buono (2026)
###############################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tcltk)
  library(purrr)
})

# ----------------------------- Optional theme nudge --------------------------
# On macOS this helps ttk respect Aqua rendering.
try(tcl("ttk::style", "theme", "use", "aqua"), silent = TRUE)
# ---------------------------------------------------------------------------

# ----------------------------- Helpers ---------------------------------------

os_open_folder <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  sys <- Sys.info()[["sysname"]]
  if (sys == "Darwin") {
    system2("open", shQuote(path))
  } else if (sys == "Windows") {
    system2("explorer", shQuote(path), wait = FALSE)
  } else {
    system2("xdg-open", shQuote(path), wait = FALSE)
  }
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

show_error <- function(msg) tkmessageBox(title = "Error", message = msg, icon = "error", type = "ok")
show_info  <- function(msg) tkmessageBox(title = "Info",  message = msg, icon = "info",  type = "ok")

detect_delim <- function(first_line) {
  counts <- c(
    comma  = stringr::str_count(first_line, ","),
    semico = stringr::str_count(first_line, ";"),
    tab    = stringr::str_count(first_line, "\t")
  )
  delim_name <- names(which.max(counts))
  switch(delim_name, comma = ",", semico = ";", tab = "\t")
}

looks_like_header <- function(path) {
  lns <- readr::read_lines(path, n_max = 2)
  if (length(lns) == 0) return(FALSE)
  delim <- detect_delim(lns[[1]])
  tok1  <- strsplit(lns[[1]], delim, fixed = TRUE)[[1]]
  has_letters <- any(grepl("[A-Za-z]", tok1))
  if (length(lns) == 1) return(has_letters)
  tok2  <- strsplit(lns[[2]], delim, fixed = TRUE)[[1]]
  second_is_data <- mean(
    grepl("^[-+]?[0-9]*\\.?[0-9]+$", trimws(tok2)) |
      grepl("\\.dat$", trimws(tok2), ignore.case = TRUE)
  ) > 0.5
  has_letters && second_is_data
}

clean_names_safe <- function(x) {
  x <- iconv(x, to = "UTF-8")
  x <- str_replace_all(x, "[\\r\\n\\t]", " ")
  x <- str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace(x, "^_+|_+$", "")
  make.unique(x, sep = "_")
}

read_summary_file <- function(path) {
  first_line <- readr::read_lines(path, n_max = 1)
  delim <- detect_delim(first_line)
  header_present <- looks_like_header(path)
  
  if (header_present) {
    df <- readr::read_delim(
      file = path, delim = delim, col_names = TRUE,
      quote = "\"", na = c("", "NA"), trim_ws = TRUE,
      show_col_types = FALSE, guess_max = 10000
    )
    names(df) <- clean_names_safe(names(df))
  } else {
    df <- readr::read_delim(
      file = path, delim = delim, col_names = FALSE,
      quote = "\"", na = c("", "NA"), trim_ws = TRUE,
      show_col_types = FALSE, guess_max = 10000
    )
    names(df) <- paste0("Col", seq_len(ncol(df)))
  }
  
  dplyr::mutate(df, source_file = path, .before = 1)
}

# ----------------------------- Pipeline logic --------------------------------

run_step1 <- function(input_dir, output_parent_dir, Patient_ID, Week, Batch, FOLLOW_TOL, open_after = TRUE) {
  if (!dir.exists(input_dir)) stop("Input folder does not exist.")
  if (!dir.exists(output_parent_dir)) dir.create(output_parent_dir, recursive = TRUE)
  
  files_EHTs <- list.files(
    path       = input_dir,
    pattern    = "^Summary_Data\\.csv$",
    full.names = TRUE,
    recursive  = TRUE
  )
  if (!length(files_EHTs)) stop("No 'Summary_Data.csv' files found in the input folder (recursive).")
  
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")
  out_dir <- file.path(output_parent_dir, paste0("_EHT_Analysis_", timestamp))
  dir.create(out_dir, recursive = TRUE)
  
  force_data <- purrr::map_dfr(files_EHTs, read_summary_file)
  
  pacing_map <- c(`001` = 0, `002` = 1, `003` = 1.5, `004` = 2, `005` = 2.5, `006` = 3)
  
  req_cols <- c("File_Name", "Frequency_Hz")
  missing <- setdiff(req_cols, names(force_data))
  if (length(missing)) {
    stop(
      "Required column(s) missing: ", paste(missing, collapse = ", "),
      "\nAvailable columns:\n", paste(names(force_data), collapse = ", ")
    )
  }
  
  force_data <- force_data %>%
    mutate(
      Batch = Batch,
      Donor = Patient_ID,
      Week_Time = Week,
      step_code    = stringr::str_match(File_Name, "(?i)[_-]step[_-](\\d{3})")[, 2],
      Expected_Hz  = as.numeric(pacing_map[step_code]),
      Deviation_Hz = Frequency_Hz - Expected_Hz,
      Rel_Error    = dplyr::if_else(Expected_Hz == 0, NA_real_, Deviation_Hz / Expected_Hz)
    ) %>%
    mutate(
      Following = dplyr::case_when(
        is.na(Expected_Hz) | is.na(Frequency_Hz) ~ NA_character_,
        Expected_Hz == 0                         ~ dplyr::if_else(Frequency_Hz == 0, "Yes", "No"),
        abs(Deviation_Hz) <= FOLLOW_TOL * Expected_Hz ~ "Yes",
        TRUE                                      ~ "No"
      )
    ) %>%
    relocate(Batch, .before = 1) %>%
    relocate(Donor, .after = File_Name) %>%
    relocate(Week_Time, .after = Following)
  
  parent_folder <- basename(dirname(out_dir))
  prefix_tag <- substr(parent_folder, 1, 2)
  
  out_file <- file.path(
    out_dir,
    paste0("Combined_Summary_with_Following_", prefix_tag, "_tol", as.integer(FOLLOW_TOL * 100), "pct.csv")
  )
  readr::write_csv(force_data, out_file)
  
  if (open_after) try(os_open_folder(out_dir), silent = TRUE)
  
  list(
    out_dir = out_dir,
    out_file = out_file,
    n_files = length(files_EHTs),
    n_rows = nrow(force_data),
    following_table = table(force_data$Following, useNA = "ifany")
  )
}

# ----------------------------- GUI (ttk) -------------------------------------

build_step1_gui <- function() {
  tt <- tktoplevel()
  tkwm.title(tt, "EHT Step 1 — Build Combined + Following (ttk)")
  
  # vars
  v_input  <- tclVar("")
  v_output <- tclVar("")
  v_patient <- tclVar("P2")
  v_week    <- tclVar("2")
  v_batch   <- tclVar(format(Sys.Date(), "%y%m%d"))
  v_tol_mode <- tclVar("10%")
  v_tol_custom <- tclVar("0.10")
  v_open_after <- tclVar("1")
  
  row <- 0
  
  # Input folder
  tkgrid(ttklabel(tt, text = "Input folder (contains Summary_Data.csv recursively):"),
         row = row, column = 0, sticky = "w", padx = 8, pady = 6)
  tkgrid(ttkentry(tt, textvariable = v_input, width = 60),
         row = row, column = 1, sticky = "we", padx = 8, pady = 6)
  tkgrid(ttkbutton(tt, text = "Browse…", command = function() {
    p <- tk_choose.dir(caption = "Select input folder (contains Summary_Data.csv files)")
    if (!is.na(p) && nchar(p) > 0) tclvalue(v_input) <- p
  }), row = row, column = 2, padx = 8, pady = 6)
  row <- row + 1
  
  # Output folder
  tkgrid(ttklabel(tt, text = "Output parent folder (timestamped subfolder will be created):"),
         row = row, column = 0, sticky = "w", padx = 8, pady = 6)
  tkgrid(ttkentry(tt, textvariable = v_output, width = 60),
         row = row, column = 1, sticky = "we", padx = 8, pady = 6)
  tkgrid(ttkbutton(tt, text = "Browse…", command = function() {
    p <- tk_choose.dir(caption = "Select output folder")
    if (!is.na(p) && nchar(p) > 0) tclvalue(v_output) <- p
  }), row = row, column = 2, padx = 8, pady = 6)
  row <- row + 1
  
  # Patient / Week / Batch
  tkgrid(ttklabel(tt, text = "Patient_ID:"), row = row, column = 0, sticky = "w", padx = 8, pady = 6)
  tkgrid(ttkentry(tt, textvariable = v_patient, width = 20),
         row = row, column = 1, sticky = "w", padx = 8, pady = 6)
  row <- row + 1
  
  tkgrid(ttklabel(tt, text = "Week (Week_Time):"), row = row, column = 0, sticky = "w", padx = 8, pady = 6)
  tkgrid(ttkentry(tt, textvariable = v_week, width = 20),
         row = row, column = 1, sticky = "w", padx = 8, pady = 6)
  row <- row + 1
  
  tkgrid(ttklabel(tt, text = "Batch (yymmdd):"), row = row, column = 0, sticky = "w", padx = 8, pady = 6)
  tkgrid(ttkentry(tt, textvariable = v_batch, width = 20),
         row = row, column = 1, sticky = "w", padx = 8, pady = 6)
  row <- row + 1
  
  # Tolerance (combobox + custom)
  tkgrid(ttklabel(tt, text = "Following tolerance:"), row = row, column = 0, sticky = "w", padx = 8, pady = 6)
  tol_frame <- ttkframe(tt)
  tkgrid(tol_frame, row = row, column = 1, sticky = "w", padx = 8, pady = 6)
  
  tol_choices <- c("5%", "10%", "15%", "Custom")
  cb <- ttkcombobox(tol_frame, values = tol_choices, textvariable = v_tol_mode, width = 8, state = "readonly")
  tcl(cb, "set", "10%")
  tkgrid(cb, row = 0, column = 0, padx = 4)
  
  tkgrid(ttklabel(tol_frame, text = "Custom (e.g. 0.07):"), row = 0, column = 1, padx = 8)
  tkgrid(ttkentry(tol_frame, textvariable = v_tol_custom, width = 10),
         row = 0, column = 2, padx = 4)
  row <- row + 1
  
  # Open after
  tkgrid(ttkcheckbutton(tt, text = "Open output folder after run", variable = v_open_after),
         row = row, column = 1, sticky = "w", padx = 8, pady = 6)
  row <- row + 1
  
  # Log (tktext is fine; if it still looks bad we can make it separate or remove)
  tkgrid(ttklabel(tt, text = "Log:"), row = row, column = 0, sticky = "nw", padx = 8, pady = 6)
  status <- tktext(tt, height = 8, width = 85)
  tkgrid(status, row = row, column = 1, columnspan = 2, sticky = "we", padx = 8, pady = 6)
  row <- row + 1
  
  log_line <- function(x) {
    tkinsert(status, "end", paste0(x, "\n"))
    tksee(status, "end")
    tcl("update")
  }
  
  # Buttons
  btn_frame <- ttkframe(tt)
  tkgrid(btn_frame, row = row, column = 1, sticky = "w", padx = 8, pady = 10)
  
  tkgrid(ttkbutton(btn_frame, text = "Run Step 1", command = function() {
    in_dir  <- tclvalue(v_input)
    out_dir <- tclvalue(v_output)
    patient <- tclvalue(v_patient)
    week    <- tclvalue(v_week)
    batch   <- tclvalue(v_batch)
    
    if (nchar(in_dir) == 0 || !dir.exists(in_dir)) return(show_error("Please select a valid input folder."))
    if (nchar(out_dir) == 0) return(show_error("Please select an output folder."))
    if (nchar(patient) == 0) return(show_error("Patient_ID cannot be empty."))
    if (nchar(week) == 0) return(show_error("Week cannot be empty."))
    if (!grepl("^\\d{6}$", batch)) return(show_error("Batch should be yymmdd (6 digits), e.g. 251124."))
    
    tol_mode <- tclvalue(v_tol_mode)
    tol <- switch(tol_mode,
                  "5%" = 0.05,
                  "10%" = 0.10,
                  "15%" = 0.15,
                  "Custom" = safe_num(tclvalue(v_tol_custom))
    )
    if (is.na(tol) || tol <= 0 || tol >= 1) return(show_error("FOLLOW_TOL must be a number between 0 and 1."))
    
    open_after <- (tclvalue(v_open_after) == "1")
    
    log_line("------------------------------------------------------------")
    log_line(paste0("Running Step 1 with: Patient_ID=", patient,
                    " | Week=", week, " | Batch=", batch, " | FOLLOW_TOL=", tol))
    log_line(paste0("Input:  ", in_dir))
    log_line(paste0("Output: ", out_dir))
    
    res <- try(run_step1(in_dir, out_dir, patient, week, batch, tol, open_after), silent = TRUE)
    if (inherits(res, "try-error")) {
      show_error(as.character(res))
      log_line(paste0("❌ Error: ", as.character(res)))
      return(invisible(NULL))
    }
    
    log_line(paste0("✅ Done. Files parsed: ", res$n_files))
    log_line(paste0("✅ Rows written: ", res$n_rows))
    log_line(paste0("✅ Output file: ", res$out_file))
    log_line("Following table:")
    log_line(paste(capture.output(print(res$following_table)), collapse = "\n"))
    
    show_info(paste0("Step 1 completed.\n\nOutput:\n", res$out_file))
  }), row = 0, column = 0, padx = 6)
  
  tkgrid(ttkbutton(btn_frame, text = "Close", command = function() tkdestroy(tt)),
         row = 0, column = 1, padx = 6)
  
  tkgrid.columnconfigure(tt, 1, weight = 1)
  invisible(NULL)
}

build_step1_gui()

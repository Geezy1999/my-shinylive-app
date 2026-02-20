library(shiny)
library(shinyjs)

# ============================================================
# Shared INPUTS (must exist in your environment)
# ============================================================
W <- as.matrix(read.csv("data/wts.csv"))[, -1]
colnames(W) <- c("Momentum", "Quality", "Value")

SAA <- c(0.4, 0.25, 0.35)
names(SAA) <- c("Momentum", "Quality", "Value")

# For App 2
M <- as.matrix(read.csv("data/mat_idxs.csv"))[, -1]
colnames(M) <- c("P(Alpha>AT)", "Information Ratio", "P(Alpha<0|12M)", "Batting Avg", "P(Alpha<-5|12M)")
stopifnot(is.matrix(M), nrow(M) == nrow(W))
stopifnot(!is.null(colnames(M)))
stopifnot(all(is.finite(M)))

metric_dir <- c(1, 1, -1, 1, -1)  # <-- must match ncol(M)
stopifnot(is.numeric(metric_dir), length(metric_dir) == ncol(M))
stopifnot(all(metric_dir %in% c(-1, 1)))
names(metric_dir) <- colnames(M)

# Shared checks
stopifnot(is.matrix(W), ncol(W) == 3)
stopifnot(all(is.finite(W)))
stopifnot(all(abs(rowSums(W) - 1) < 1e-8))
stopifnot(is.numeric(SAA), length(SAA) == 3)

# Shared helpers
labs  <- colnames(W)
W_min <- apply(W, 2, min, na.rm = TRUE)
W_max <- apply(W, 2, max, na.rm = TRUE)

# Display-only rounding helpers
round_up_2dp   <- function(x) ceiling(x * 100) / 100
round_down_2dp <- function(x) floor(x * 100) / 100
W_min_disp <- round_up_2dp(W_min)
W_max_disp <- round_down_2dp(W_max)

STEP_2DP <- 0.01

# ============================================================
# Shared UI styling (robust lock + red input)
# ============================================================
APP_CSS <- "
.locked-wrap { opacity: 0.45 !important; pointer-events: none !important; }
input.out-of-range { border: 2px solid #d11 !important; background-color: #ffecec !important; }
.small-error { color: #d11; font-size: 12px; margin-top: -8px; margin-bottom: 10px; }
"

# ============================================================
# App 1 UI
# ============================================================
app1UI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    tags$style(HTML(APP_CSS)),
    titlePanel("Factor Portfolio Positioning Explorer"),
    tags$p(
      "Explore feasible factor tilts given long-term portfolio return objectives and risk constraints",
      style = "margin-top:-10px; margin-bottom:18px; color:#666; font-size:14px;"
    ),
    sidebarLayout(
      sidebarPanel(
        tags$div(
          id = ns("mom_wrap"),
          numericInput(ns("mom"), "Momentum", value = NULL,
                       min = W_min["Momentum"], max = W_max["Momentum"], step = STEP_2DP),
          uiOutput(ns("mom_err")),
          helpText(sprintf("SAA: %.2f | Allowed range: [%.2f, %.2f]",
                           SAA["Momentum"], W_min_disp["Momentum"], W_max_disp["Momentum"]))
        ),
        tags$div(
          id = ns("qua_wrap"),
          numericInput(ns("qua"), "Quality", value = NULL,
                       min = W_min["Quality"], max = W_max["Quality"], step = STEP_2DP),
          uiOutput(ns("qua_err")),
          helpText(sprintf("SAA: %.2f | Allowed range: [%.2f, %.2f]",
                           SAA["Quality"], W_min_disp["Quality"], W_max_disp["Quality"]))
        ),
        tags$div(
          id = ns("val_wrap"),
          numericInput(ns("val"), "Value", value = NULL,
                       min = W_min["Value"], max = W_max["Value"], step = STEP_2DP),
          uiOutput(ns("val_err")),
          helpText(sprintf("SAA: %.2f | Allowed range: [%.2f, %.2f]",
                           SAA["Value"], W_min_disp["Value"], W_max_disp["Value"]))
        ),
        tags$hr(),
        actionButton(ns("reset"), "Reset (clear all)"),
        tags$hr(),
        htmlOutput(ns("status"))
      ),
      mainPanel(plotOutput(ns("bxp_plot"), height = "520px"))
    )
  )
}

# ============================================================
# App 1 Server
# ============================================================
app1Server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Static boxplot object
    stats_full <- sapply(seq_len(ncol(W)), function(j) {
      quantile(W[, j], probs = c(0.001, 0.25, 0.50, 0.75, 0.999), na.rm = TRUE, names = FALSE)
    })
    bxp_obj_full <- list(stats = stats_full, n = colSums(is.finite(W)), names = labs)
    
    K_NEIGH <- 20
    
    neighbour_idx_1d <- function(W, driver, target, k = 20) {
      j <- match(driver, colnames(W))
      ord <- order(abs(W[, j] - target))
      ord[seq_len(min(k, nrow(W)))]
    }
    
    mini_bxp_from_neighbours <- function(W, driver, target, k = 20, probs = c(0, 0.25, 0.50, 0.75, 1)) {
      nn <- neighbour_idx_1d(W, driver = driver, target = target, k = k)
      other <- setdiff(colnames(W), driver)
      stats_mini <- sapply(other, function(nm) quantile(W[nn, nm], probs = probs, na.rm = TRUE, names = FALSE))
      stats_mini <- as.matrix(stats_mini)
      list(
        nn = nn,
        other = other,
        bxp = list(stats = stats_mini, n = rep(length(nn), length(other)), names = other)
      )
    }
    
    draw_mini_box <- function(stats, x, halfwidth = 0.14,
                              boxfill = rgb(198/255, 226/255, 255/255, alpha = 0.35),
                              border = "navy", whiskcol = "navy", staplecol = "navy", lwd = 2) {
      stats <- as.numeric(stats)
      if (length(stats) != 5 || any(!is.finite(stats))) return(invisible(NULL))
      lo <- stats[1]; q1 <- stats[2]; med <- stats[3]; q3 <- stats[4]; hi <- stats[5]
      segments(x, lo, x, hi, col = whiskcol, lwd = lwd)
      segments(x - halfwidth * 0.6, lo, x + halfwidth * 0.6, lo, col = staplecol, lwd = lwd)
      segments(x - halfwidth * 0.6, hi, x + halfwidth * 0.6, hi, col = staplecol, lwd = lwd)
      rect(x - halfwidth, q1, x + halfwidth, q3, col = boxfill, border = border, lwd = lwd)
      segments(x - halfwidth, med, x + halfwidth, med, col = border, lwd = lwd + 0.5)
      invisible(NULL)
    }
    
    # State
    active_driver <- reactiveVal(NULL)
    mini_info <- reactiveVal(NULL)
    err_msg <- reactiveValues(mom = NULL, qua = NULL, val = NULL)
    
    # Locking that actually works with numericInput
    set_input_locked <- function(which, locked) {
      input_id <- session$ns(which)
      wrap_id  <- session$ns(paste0(which, "_wrap"))
      if (locked) {
        shinyjs::runjs(sprintf("$('#%s').prop('disabled', true);", input_id))
        shinyjs::addClass(wrap_id, "locked-wrap")
      } else {
        shinyjs::runjs(sprintf("$('#%s').prop('disabled', false);", input_id))
        shinyjs::removeClass(wrap_id, "locked-wrap")
      }
    }
    
    set_out_of_range <- function(which, is_bad) {
      input_id <- session$ns(which)
      if (is_bad) shinyjs::addClass(input_id, "out-of-range") else shinyjs::removeClass(input_id, "out-of-range")
    }
    
    validate_target <- function(drv, target) {
      if (is.null(target) || is.na(target)) return(list(ok = FALSE, msg = NULL))
      if (target < W_min[drv] || target > W_max[drv]) {
        return(list(ok = FALSE, msg = sprintf("Out of range: must be in [%.2f, %.2f]", W_min_disp[drv], W_max_disp[drv])))
      }
      list(ok = TRUE, msg = NULL)
    }
    
    # HARD clear that actually blanks the browser field
    hard_clear_input <- function(which) {
      id <- session$ns(which)
      shinyjs::runjs(sprintf("$('#%s').val(''); $('#%s').trigger('input'); $('#%s').trigger('change');", id, id, id))
    }
    
    observeEvent(input$reset, {
      # 1) reset state first
      active_driver(NULL)
      mini_info(NULL)
      
      # 2) clear errors
      err_msg$mom <- NULL; err_msg$qua <- NULL; err_msg$val <- NULL
      
      # 3) unlock first (so browser accepts the clearing)
      set_input_locked("mom", FALSE)
      set_input_locked("qua", FALSE)
      set_input_locked("val", FALSE)
      
      # 4) remove red highlight
      set_out_of_range("mom", FALSE)
      set_out_of_range("qua", FALSE)
      set_out_of_range("val", FALSE)
      
      # 5) freeze to avoid observers re-firing during clear
      freezeReactiveValue(input, "mom")
      freezeReactiveValue(input, "qua")
      freezeReactiveValue(input, "val")
      
      # 6) clear BOTH: Shiny value + browser UI value
      updateNumericInput(session, "mom", value = NULL)
      updateNumericInput(session, "qua", value = NULL)
      updateNumericInput(session, "val", value = NULL)
      
      hard_clear_input("mom")
      hard_clear_input("qua")
      hard_clear_input("val")
      
    }, ignoreInit = TRUE)
    
    # Lock after first selection
    observeEvent(input$mom, {
      if (is.null(active_driver()) && !is.null(input$mom) && !is.na(input$mom)) {
        active_driver("Momentum")
        set_input_locked("qua", TRUE)
        set_input_locked("val", TRUE)
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$qua, {
      if (is.null(active_driver()) && !is.null(input$qua) && !is.na(input$qua)) {
        active_driver("Quality")
        set_input_locked("mom", TRUE)
        set_input_locked("val", TRUE)
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$val, {
      if (is.null(active_driver()) && !is.null(input$val) && !is.na(input$val)) {
        active_driver("Value")
        set_input_locked("mom", TRUE)
        set_input_locked("qua", TRUE)
      }
    }, ignoreInit = TRUE)
    
    # Compute mini info if valid
    observeEvent(list(input$mom, input$qua, input$val), {
      drv <- active_driver()
      if (is.null(drv)) { mini_info(NULL); return() }
      
      target <- switch(drv,
                       "Momentum" = input$mom,
                       "Quality"  = input$qua,
                       "Value"    = input$val)
      
      v <- validate_target(drv, target)
      
      err_msg$mom <- NULL; err_msg$qua <- NULL; err_msg$val <- NULL
      set_out_of_range("mom", FALSE); set_out_of_range("qua", FALSE); set_out_of_range("val", FALSE)
      
      active_input_id <- switch(drv, "Momentum" = "mom", "Quality" = "qua", "Value" = "val")
      
      if (!isTRUE(v$ok)) {
        if (!is.null(v$msg)) {
          err_msg[[active_input_id]] <- v$msg
          set_out_of_range(active_input_id, TRUE)
        }
        mini_info(NULL)
        return()
      }
      
      info <- mini_bxp_from_neighbours(W, driver = drv, target = target, k = K_NEIGH)
      info$target <- target
      mini_info(info)
    }, ignoreInit = TRUE)
    
    output$mom_err <- renderUI(if (is.null(err_msg$mom)) NULL else tags$div(class = "small-error", err_msg$mom))
    output$qua_err <- renderUI(if (is.null(err_msg$qua)) NULL else tags$div(class = "small-error", err_msg$qua))
    output$val_err <- renderUI(if (is.null(err_msg$val)) NULL else tags$div(class = "small-error", err_msg$val))
    
    output$status <- renderUI({
      drv <- active_driver()
      info <- mini_info()
      
      lines <- c()
      lines <- c(lines, sprintf("Active driver: %s", if (is.null(drv)) "(none yet)" else drv))
      lines <- c(lines, "")
      lines <- c(lines, "Inputs (what you typed; other boxes remain unchanged):")
      lines <- c(lines, sprintf("  Momentum: %s", ifelse(is.null(input$mom), "NULL", as.character(input$mom))))
      lines <- c(lines, sprintf("  Quality : %s", ifelse(is.null(input$qua), "NULL", as.character(input$qua))))
      lines <- c(lines, sprintf("  Value   : %s", ifelse(is.null(input$val), "NULL", as.character(input$val))))
      lines <- c(lines, "")
      lines <- c(lines, sprintf("Mini boxplot based on locally similar portfolios (N = %d)", K_NEIGH))
      
      if (is.null(drv) || is.null(info)) {
        if (!is.null(drv)) {
          active_input_id <- switch(drv, "Momentum" = "mom", "Quality" = "qua", "Value" = "val")
          if (!is.null(err_msg[[active_input_id]])) {
            lines <- c(lines, "")
            lines <- c(lines, sprintf("NOTE: %s", err_msg[[active_input_id]]))
          } else {
            lines <- c(lines, "(none yet — enter a value in one box)")
          }
        } else {
          lines <- c(lines, "(none yet — enter a value in one box)")
        }
        return(tags$pre(HTML(paste(lines, collapse = "\n"))))
      }
      
      lines <- c(lines, "")
      lines <- c(lines, sprintf("Driver target: %s = %.4f", drv, info$target))
      lines <- c(lines, sprintf("Other factors: %s", paste(info$other, collapse = ", ")))
      
      # ---- ADD: mini-boxplot min/median/max for the unspecified factors ----
      stats_m <- as.matrix(info$bxp$stats)  # 5 x length(info$other)
      if (nrow(stats_m) == 5 && ncol(stats_m) == length(info$other)) {
        
        mins <- round(stats_m[1, ], 2)
        meds <- round(stats_m[3, ], 2)
        maxs <- round(stats_m[5, ], 2)
        
        names(mins) <- info$other
        names(meds) <- info$other
        names(maxs) <- info$other
        
        lines <- c(lines, "")
        lines <- c(lines, "Dynamic (local) tilt ranges from mini boxplots:")
        lines <- c(lines, paste0("  Min   : ", paste(sprintf("%s=%.2f", names(mins), mins), collapse = " | ")))
        lines <- c(lines, paste0("  Median: ", paste(sprintf("%s=%.2f", names(meds), meds), collapse = " | ")))
        lines <- c(lines, paste0("  Max   : ", paste(sprintf("%s=%.2f", names(maxs), maxs), collapse = " | ")))
      } else {
        lines <- c(lines, "")
        lines <- c(lines, "Dynamic (local) tilt ranges: (unavailable)")
      }
      
      tags$pre(HTML(paste(lines, collapse = "\n")))
    })
  
  output$bxp_plot <- renderPlot({
    bxp(
      bxp_obj_full,
      main = "Boxplot of TAA Ranges",
      boxfill   = "grey95",
      border    = "grey80",
      whiskcol  = "grey80",
      staplecol = "grey80",
      lwd       = 2.5,
      outline   = FALSE
    )
    
    points(1:3, SAA[1:3], pch = 21, bg = "red",  col = "red",   cex = 1.65, lwd = 2)
    points(1:3, SAA[1:3], pch = 21, bg = "red",  col = "white", cex = 1.30, lwd = 2)
    
    info <- mini_info()
    if (!is.null(info)) {
      drv_pos <- match(active_driver(), labs)
      if (is.finite(drv_pos)) points(drv_pos, info$target, pch = 23, bg = "blue", col = "white", cex = 1.2)
      
      at_pos <- match(info$other, labs)
      stats_m <- as.matrix(info$bxp$stats)
      for (i in seq_along(at_pos)) draw_mini_box(stats = stats_m[, i], x = at_pos[i], halfwidth = 0.14)
      
      legend("topright",
             legend = c("Strategic factor allocation", "Dynamic (local) tilt ranges"),
             pch = c(16, 15), col = c("red", "blue"),
             pt.cex = c(1.4, 1.2), bty = "n")
    } else {
      legend("topright", legend = c("Strategic factor allocation"),
             pch = 16, col = "red", pt.cex = 1.4, bty = "n")
    }
  })
  })
}

# ============================================================
# App 2 UI  (REPLACE YOUR WHOLE app2UI WITH THIS)
# ============================================================
app2UI <- function(id) {
  ns <- NS(id)
  
  metric_choices <- setNames(
    colnames(M),
    paste0(colnames(M), ifelse(metric_dir[colnames(M)] == 1, " (max)", " (min)"))
  )
  
  fluidPage(
    tags$style(HTML(APP_CSS)),
    titlePanel("Factor Portfolio Positioning Explorer"),
    tags$p(
      "Suggesting complimentary factor tilts based on risk/return objectives",
      style = "margin-top:-10px; margin-bottom:18px; color:#666; font-size:14px;"
    ),
    sidebarLayout(
      sidebarPanel(
        radioButtons(
          ns("method"), "How should the portfolio be selected?",
          choices = c(
            "Average factor tilts for nearby portfolios" = "nn_mean",
            "Best portfolio for chosen objective (green point)" = "metric_opt"
          ),
          selected = "nn_mean"
        ),
        selectInput(ns("metric"), "Which objective should be prioritised?",
                    choices = metric_choices, selected = colnames(M)[1]),
        tags$hr(),
        
        tags$div(
          id = ns("mom_wrap"),
          numericInput(ns("mom"), "Momentum", value = NULL,
                       min = W_min["Momentum"], max = W_max["Momentum"], step = STEP_2DP),
          uiOutput(ns("mom_err")),
          helpText(sprintf("SAA: %.2f | Allowed range: [%.2f, %.2f]",
                           SAA["Momentum"], W_min_disp["Momentum"], W_max_disp["Momentum"]))
        ),
        tags$div(
          id = ns("qua_wrap"),
          numericInput(ns("qua"), "Quality", value = NULL,
                       min = W_min["Quality"], max = W_max["Quality"], step = STEP_2DP),
          uiOutput(ns("qua_err")),
          helpText(sprintf("SAA: %.2f | Allowed range: [%.2f, %.2f]",
                           SAA["Quality"], W_min_disp["Quality"], W_max_disp["Quality"]))
        ),
        tags$div(
          id = ns("val_wrap"),
          numericInput(ns("val"), "Value", value = NULL,
                       min = W_min["Value"], max = W_max["Value"], step = STEP_2DP),
          uiOutput(ns("val_err")),
          helpText(sprintf("SAA: %.2f | Allowed range: [%.2f, %.2f]",
                           SAA["Value"], W_min_disp["Value"], W_max_disp["Value"]))
        ),
        
        tags$hr(),
        actionButton(ns("reset"), "Reset (clear all)"),
        tags$hr(),
        
        # CHANGED: status now sits inside a box with a dynamic title
        tags$div(
          style = "border:1px solid #ddd; border-radius:6px; padding:10px; background:#fafafa;",
          tags$div(
            style = "font-weight:600; margin-bottom:8px;",
            uiOutput(ns("console_title"))
          ),
          htmlOutput(ns("status"))
        )
      ),
      mainPanel(plotOutput(ns("bxp_plot"), height = "520px"))
    )
  )
}

# ============================================================
# App 2 Server  (REPLACE YOUR WHOLE app2Server WITH THIS)
# ============================================================
# ============================================================
# App 2 Server  (FULL MODULE - COPY/PASTE)
# ============================================================
# ============================================================
# App 2 Server  (FULL MODULE - COPY/PASTE)
#   - Hides ALL objective-value text when nn_mean is selected
#   - Disables the metric dropdown when nn_mean is selected
# ============================================================
app2Server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    stats_full <- sapply(seq_len(ncol(W)), function(j) {
      quantile(W[, j], probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE, names = FALSE)
    })
    bxp_obj_full <- list(stats = stats_full, n = colSums(is.finite(W)), names = labs)
    
    K_NEIGH <- 20
    
    neighbour_idx_1d <- function(W, driver, target, k = 20) {
      j <- match(driver, colnames(W))
      ord <- order(abs(W[, j] - target))
      ord[seq_len(min(k, nrow(W)))]
    }
    
    adjust_from_neighbours_no_sum1 <- function(W, driver, target, k = 20) {
      nn <- neighbour_idx_1d(W, driver = driver, target = target, k = k)
      out <- numeric(3); names(out) <- colnames(W)
      out[driver] <- target
      other <- setdiff(colnames(W), driver)
      out[other] <- colMeans(W[nn, other, drop = FALSE], na.rm = TRUE)
      out
    }
    
    best_idx_by_metric <- function(M, metric, metric_dir, idx = NULL) {
      if (is.null(idx)) idx <- seq_len(nrow(M))
      vals <- M[idx, metric]
      if (metric_dir[metric] == 1) idx[which.max(vals)] else idx[which.min(vals)]
    }
    
    nearest_idx_in_W <- function(W, p) {
      p <- as.numeric(p)
      d2 <- rowSums((W - matrix(p, nrow = nrow(W), ncol = ncol(W), byrow = TRUE))^2)
      which.min(d2)
    }
    
    active_driver <- reactiveVal(NULL)
    adjusted <- reactiveVal(NULL)
    best_idx <- reactiveVal(NULL)
    err_msg <- reactiveValues(mom = NULL, qua = NULL, val = NULL)
    
    set_input_locked <- function(which, locked) {
      input_id <- session$ns(which)
      wrap_id  <- session$ns(paste0(which, "_wrap"))
      if (locked) {
        shinyjs::runjs(sprintf("$('#%s').prop('disabled', true);", input_id))
        shinyjs::addClass(wrap_id, "locked-wrap")
      } else {
        shinyjs::runjs(sprintf("$('#%s').prop('disabled', false);", input_id))
        shinyjs::removeClass(wrap_id, "locked-wrap")
      }
    }
    
    set_out_of_range <- function(which, is_bad) {
      input_id <- session$ns(which)
      if (is_bad) shinyjs::addClass(input_id, "out-of-range") else shinyjs::removeClass(input_id, "out-of-range")
    }
    
    validate_target <- function(drv, target) {
      if (is.null(target) || is.na(target)) return(list(ok = FALSE, msg = NULL))
      if (target < W_min[drv] || target > W_max[drv]) {
        return(list(ok = FALSE, msg = sprintf("Out of range: must be in [%.2f, %.2f]",
                                              W_min_disp[drv], W_max_disp[drv])))
      }
      list(ok = TRUE, msg = NULL)
    }
    
    hard_clear_input <- function(which) {
      id <- session$ns(which)
      shinyjs::runjs(sprintf("$('#%s').val(''); $('#%s').trigger('input'); $('#%s').trigger('change');", id, id, id))
    }
    
    observeEvent(input$reset, {
      active_driver(NULL)
      adjusted(NULL)
      best_idx(NULL)
      
      err_msg$mom <- NULL; err_msg$qua <- NULL; err_msg$val <- NULL
      
      set_input_locked("mom", FALSE)
      set_input_locked("qua", FALSE)
      set_input_locked("val", FALSE)
      
      set_out_of_range("mom", FALSE)
      set_out_of_range("qua", FALSE)
      set_out_of_range("val", FALSE)
      
      freezeReactiveValue(input, "mom")
      freezeReactiveValue(input, "qua")
      freezeReactiveValue(input, "val")
      
      updateNumericInput(session, "mom", value = NULL)
      updateNumericInput(session, "qua", value = NULL)
      updateNumericInput(session, "val", value = NULL)
      
      hard_clear_input("mom")
      hard_clear_input("qua")
      hard_clear_input("val")
    }, ignoreInit = TRUE)
    
    observeEvent(input$mom, {
      if (is.null(active_driver()) && !is.null(input$mom) && !is.na(input$mom)) {
        active_driver("Momentum")
        set_input_locked("qua", TRUE)
        set_input_locked("val", TRUE)
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$qua, {
      if (is.null(active_driver()) && !is.null(input$qua) && !is.na(input$qua)) {
        active_driver("Quality")
        set_input_locked("mom", TRUE)
        set_input_locked("val", TRUE)
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$val, {
      if (is.null(active_driver()) && !is.null(input$val) && !is.na(input$val)) {
        active_driver("Value")
        set_input_locked("mom", TRUE)
        set_input_locked("qua", TRUE)
      }
    }, ignoreInit = TRUE)
    
    # ------------------------------------------------------------
    # NEW: Disable metric dropdown when averaging (nn_mean)
    # ------------------------------------------------------------
    observeEvent(input$method, {
      shinyjs::toggleState(id = session$ns("metric"), condition = identical(input$method, "metric_opt"))
    }, ignoreInit = FALSE)
    
    observeEvent(list(input$mom, input$qua, input$val, input$method, input$metric), {
      drv <- active_driver()
      
      if (is.null(drv)) {
        adjusted(NULL)
        if (identical(input$method, "metric_opt")) {
          best_idx(best_idx_by_metric(M, metric = input$metric, metric_dir = metric_dir, idx = NULL))
        } else {
          best_idx(NULL)
        }
        return()
      }
      
      target <- switch(drv, "Momentum" = input$mom, "Quality" = input$qua, "Value" = input$val)
      v <- validate_target(drv, target)
      
      err_msg$mom <- NULL; err_msg$qua <- NULL; err_msg$val <- NULL
      set_out_of_range("mom", FALSE); set_out_of_range("qua", FALSE); set_out_of_range("val", FALSE)
      
      active_input_id <- switch(drv, "Momentum" = "mom", "Quality" = "qua", "Value" = "val")
      
      if (!isTRUE(v$ok)) {
        if (!is.null(v$msg)) {
          err_msg[[active_input_id]] <- v$msg
          set_out_of_range(active_input_id, TRUE)
        }
        adjusted(NULL); best_idx(NULL)
        return()
      }
      
      nn <- neighbour_idx_1d(W, driver = drv, target = target, k = K_NEIGH)
      
      if (identical(input$method, "nn_mean")) {
        adjusted(adjust_from_neighbours_no_sum1(W, driver = drv, target = target, k = K_NEIGH))
        best_idx(NULL)
      } else {
        adjusted(NULL)
        best_idx(best_idx_by_metric(M, metric = input$metric, metric_dir = metric_dir, idx = nn))
      }
    }, ignoreInit = TRUE)
    
    output$mom_err <- renderUI(if (is.null(err_msg$mom)) NULL else tags$div(class = "small-error", err_msg$mom))
    output$qua_err <- renderUI(if (is.null(err_msg$qua)) NULL else tags$div(class = "small-error", err_msg$qua))
    output$val_err <- renderUI(if (is.null(err_msg$val)) NULL else tags$div(class = "small-error", err_msg$val))
    
    output$console_title <- renderUI({
      if (identical(input$method, "nn_mean")) return(tags$span("Console — Average Portfolio"))
      metric_nm <- input$metric
      if (is.null(metric_nm) || is.na(metric_nm) || !nzchar(metric_nm)) return(tags$span("Console"))
      dir_short <- if (metric_dir[metric_nm] == 1) "max" else "min"
      tags$span(paste0("Console — ", metric_nm, " (", dir_short, ")"))
    })
    
    output$status <- renderUI({
      drv <- active_driver()
      
      is_avg <- identical(input$method, "nn_mean")
      method_txt <- if (is_avg) "Average portfolio (blue)" else "Best portfolio (green)"
      
      metric_nm <- input$metric
      dir_txt   <- if (metric_dir[metric_nm] == 1) "Maximise" else "Minimise"
      
      p  <- adjusted()
      bi <- best_idx()
      
      scenario_wts <- NULL
      scenario_obj <- NA
      
      if (!is.null(p)) {
        scenario_wts <- p
      } else if (!is.null(bi)) {
        scenario_wts <- as.numeric(W[bi, ])
        scenario_obj <- M[bi, metric_nm]
      }
      
      lines <- c()
      lines <- c(lines, sprintf("Active driver: %s", if (is.null(drv)) "(none yet)" else drv))
      lines <- c(lines, sprintf("Selection: %s", method_txt))
      lines <- c(lines, sprintf("Chosen objective: <b>%s</b>%s",
                                if (is_avg) "Average Portfolio" else metric_nm,
                                if (is_avg) "" else paste0("  [", dir_txt, "]")))
      lines <- c(lines, "")
      lines <- c(lines, sprintf("<b>Recommended portfolio weights</b> (%s):", method_txt))
      
      if (is.null(scenario_wts)) {
        lines <- c(lines, "  (none)")
      } else {
        wtxt <- paste0(colnames(W), "=", sprintf("%.4f", as.numeric(scenario_wts)), collapse = ", ")
        lines <- c(lines, paste0("  <b>", wtxt, "</b>"))
      }
      
      # IMPORTANT: Only show objective value lines when metric_opt is selected
      if (!is_avg) {
        saa_i   <- nearest_idx_in_W(W, SAA)
        saa_obj <- M[saa_i, metric_nm]
        
        lines <- c(lines, "")
        lines <- c(lines, sprintf("<b>%s</b> (%s) for <b>Strategic factor allocation</b>: <b>%.6f</b>",
                                  metric_nm, dir_txt, saa_obj))
        
        if (is.null(bi) || is.null(scenario_wts)) {
          lines <- c(lines, "")
          lines <- c(lines, sprintf("<b>%s</b> (%s) for scenario: <b>(none)</b>", metric_nm, dir_txt))
        } else {
          lines <- c(lines, "")
          lines <- c(lines, sprintf("<b>%s</b> (%s) for scenario: <b>%.6f</b>", metric_nm, dir_txt, scenario_obj))
        }
      }
      
      tags$pre(HTML(paste(lines, collapse = "\n")))
    })
    
    output$bxp_plot <- renderPlot({
      bxp(bxp_obj_full, main = "Boxplot of Factor Ranges", boxfill = "grey85", border = "black")
      
      points(1:3, SAA[1:3], pch = 21, bg = "red",  col = "red",   cex = 1.65, lwd = 2)
      points(1:3, SAA[1:3], pch = 21, bg = "red",  col = "white", cex = 1.30, lwd = 2)
      
      leg_labels <- c("Strategic factor allocation")
      leg_cols   <- c("red")
      leg_pch    <- c(16)
      leg_cex    <- c(1.4)
      
      p <- adjusted()
      if (!is.null(p)) {
        points(1:3, p[1:3], pch = 21, bg = "blue", col = "blue",  cex = 1.35, lwd = 2)
        points(1:3, p[1:3], pch = 21, bg = "blue", col = "white", cex = 1.05, lwd = 2)
        leg_labels <- c(leg_labels, "Average Portfolio")
        leg_cols   <- c(leg_cols, "blue")
        leg_pch    <- c(leg_pch, 16)
        leg_cex    <- c(leg_cex, 1.2)
      }
      
      bi <- best_idx()
      if (!is.null(bi)) {
        pts <- as.numeric(W[bi, ])
        points(1:3, pts, pch = 21, bg = "darkgreen", col = "darkgreen", cex = 1.35, lwd = 2)
        points(1:3, pts, pch = 21, bg = "darkgreen", col = "white",     cex = 1.05, lwd = 2)
        leg_labels <- c(leg_labels, paste0("Best portfolio for: ", input$metric))
        leg_cols   <- c(leg_cols, "darkgreen")
        leg_pch    <- c(leg_pch, 16)
        leg_cex    <- c(leg_cex, 1.2)
      }
      
      legend("topright", legend = leg_labels, pch = leg_pch, col = leg_cols, pt.cex = leg_cex, bty = "n")
    })
  })
}

# ============================================================
# Main wrapper
# ============================================================
ui <- navbarPage(
  title = "Factor Portfolio Tools",
  useShinyjs(),
  tabPanel("Tool 1: Local Tilt Ranges", app1UI("tool1")),
  tabPanel("Tool 2: Objective Optimiser", app2UI("tool2"))
)

server <- function(input, output, session) {
  app1Server("tool1")
  app2Server("tool2")
}

shinyApp(ui, server)

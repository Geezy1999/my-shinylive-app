library(shiny)
library(shinyjs)

# ============================================================
# Shared INPUTS (must exist in your environment)
# ============================================================
W <- as.matrix(read.csv("data/wts.csv"))[,-1]
colnames(W) <- c("Momentum", "Quality", "Value")

SAA <- c(0.4,0.25,0.35)
names(SAA) <- c("Momentum", "Quality", "Value")

# For App 2
M <- as.matrix(read.csv("data/mat_idxs.csv"))[,-1]
colnames(M) = c("P(Alpha>AT)","Information Ratio","P(Alpha<0|12M)","Batting Avg","P(Alpha<-5|12M)")
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

# Display-only rounding helpers (shared style)
round_up_2dp   <- function(x) ceiling(x * 100) / 100
round_down_2dp <- function(x) floor(x * 100) / 100
W_min_disp <- round_up_2dp(W_min)
W_max_disp <- round_down_2dp(W_max)

# Use 2dp spinner step in both apps (does NOT force rounding of typed values)
STEP_2DP <- 0.01

# ============================================================
# App 1 (Module): Mini boxplots from neighbours
# ============================================================
app1UI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    titlePanel("Factor Portfolio Positioning Explorer"),
    
    tags$p(
      "Explore feasible factor tilts given long-term portfolio return objectives and risk constraints",
      style = "margin-top:-10px; margin-bottom:18px; color:#666; font-size:14px;"
    ),
    
    sidebarLayout(
      sidebarPanel(
        numericInput(
          ns("mom"), "Momentum",
          value = NULL,
          min = W_min["Momentum"], max = W_max["Momentum"],
          step = STEP_2DP
        ),
        helpText(sprintf(
          "SAA: %.2f | Allowed range: [%.2f, %.2f]",
          SAA["Momentum"], W_min_disp["Momentum"], W_max_disp["Momentum"]
        )),
        
        numericInput(
          ns("qua"), "Quality",
          value = NULL,
          min = W_min["Quality"], max = W_max["Quality"],
          step = STEP_2DP
        ),
        helpText(sprintf(
          "SAA: %.2f | Allowed range: [%.2f, %.2f]",
          SAA["Quality"], W_min_disp["Quality"], W_max_disp["Quality"]
        )),
        
        numericInput(
          ns("val"), "Value",
          value = NULL,
          min = W_min["Value"], max = W_max["Value"],
          step = STEP_2DP
        ),
        helpText(sprintf(
          "SAA: %.2f | Allowed range: [%.2f, %.2f]",
          SAA["Value"], W_min_disp["Value"], W_max_disp["Value"]
        )),
        
        tags$hr(),
        actionButton(ns("reset"), "Reset (clear all)"),
        tags$hr(),
        verbatimTextOutput(ns("status"))
      ),
      mainPanel(
        plotOutput(ns("bxp_plot"), height = "520px")
      )
    )
  )
}

app1Server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # ---- Static boxplot object (App 1: 0.001/0.25/0.50/0.75/0.999)
    stats_full <- sapply(seq_len(ncol(W)), function(j) {
      quantile(
        W[, j],
        probs = c(0.001, 0.25, 0.50, 0.75, 0.999),
        na.rm = TRUE,
        names = FALSE
      )
    })
    
    bxp_obj_full <- list(
      stats = stats_full,
      n     = colSums(is.finite(W)),
      names = labs
    )
    
    # ---- Helpers (unchanged)
    K_NEIGH <- 20
    
    neighbour_idx_1d <- function(W, driver, target, k = 20) {
      j <- match(driver, colnames(W))
      ord <- order(abs(W[, j] - target))
      ord[seq_len(min(k, nrow(W)))]
    }
    
    mini_bxp_from_neighbours <- function(
    W, driver, target, k = 20,
    probs = c(0, 0.25, 0.50, 0.75, 1)
    ) {
      nn <- neighbour_idx_1d(W, driver = driver, target = target, k = k)
      other <- setdiff(colnames(W), driver)
      
      stats_mini <- sapply(other, function(nm) {
        quantile(W[nn, nm], probs = probs, na.rm = TRUE, names = FALSE)
      })
      
      stats_mini <- as.matrix(stats_mini)
      
      list(
        nn = nn,
        other = other,
        bxp = list(
          stats = stats_mini,
          n     = rep(length(nn), length(other)),
          names = other
        )
      )
    }
    
    draw_mini_box <- function(
    stats, x,
    halfwidth = 0.14,
    boxfill = rgb(198/255, 226/255, 255/255, alpha = 0.35),
    border = "navy",
    whiskcol = "navy",
    staplecol = "navy",
    lwd = 2
    ) {
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
    
    # ---- State (unchanged)
    active_driver <- reactiveVal(NULL)
    mini_info <- reactiveVal(NULL)
    
    observeEvent(input$reset, {
      active_driver(NULL)
      mini_info(NULL)
      
      updateNumericInput(session, "mom", value = NULL)
      updateNumericInput(session, "qua", value = NULL)
      updateNumericInput(session, "val", value = NULL)
      
      shinyjs::enable(session$ns("mom"))
      shinyjs::enable(session$ns("qua"))
      shinyjs::enable(session$ns("val"))
    }, ignoreInit = TRUE)
    
    observeEvent(input$mom, {
      if (is.null(active_driver()) && !is.null(input$mom) && !is.na(input$mom)) {
        active_driver("Momentum")
        shinyjs::disable(session$ns("qua"))
        shinyjs::disable(session$ns("val"))
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$qua, {
      if (is.null(active_driver()) && !is.null(input$qua) && !is.na(input$qua)) {
        active_driver("Quality")
        shinyjs::disable(session$ns("mom"))
        shinyjs::disable(session$ns("val"))
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$val, {
      if (is.null(active_driver()) && !is.null(input$val) && !is.na(input$val)) {
        active_driver("Value")
        shinyjs::disable(session$ns("mom"))
        shinyjs::disable(session$ns("qua"))
      }
    }, ignoreInit = TRUE)
    
    observeEvent(list(input$mom, input$qua, input$val), {
      
      drv <- active_driver()
      if (is.null(drv)) {
        mini_info(NULL)
        return()
      }
      
      target <- switch(
        drv,
        "Momentum" = input$mom,
        "Quality"  = input$qua,
        "Value"    = input$val
      )
      
      req(!is.null(target), !is.na(target))
      req(target >= W_min[drv], target <= W_max[drv])
      
      info <- mini_bxp_from_neighbours(W, driver = drv, target = target, k = K_NEIGH)
      info$target <- target
      mini_info(info)
      
    }, ignoreInit = TRUE)
    
    output$status <- renderPrint({
      drv <- active_driver()
      cat("Active driver:", if (is.null(drv)) "(none yet)" else drv, "\n\n")
      
      cat("Inputs (what you typed; other boxes remain unchanged):\n")
      print(c(Momentum = input$mom, Quality = input$qua, Value = input$val))
      
      cat("\nMini boxplot based on K nearest neighbours (K =", K_NEIGH, ")\n")
      info <- mini_info()
      
      if (is.null(info)) {
        cat("(none yet — enter a value in one box)\n")
        return()
      }
      
      cat("Driver target:\n")
      print(setNames(info$target, drv))
      
      cat("\nNeighbour count:", length(info$nn), "\n")
      cat("Other factors:", paste(info$other, collapse = ", "), "\n\n")
      
      stats_m <- as.matrix(info$bxp$stats)
      if (nrow(stats_m) == 5 && ncol(stats_m) == length(info$other)) {
        
        mins <- round(stats_m[1, ], 2)
        meds <- round(stats_m[3, ], 2)
        maxs <- round(stats_m[5, ], 2)
        
        names(mins) <- info$other
        names(meds) <- info$other
        names(maxs) <- info$other
        
        cat("Dynamic (local) tilt ranges from mini boxplots:\n")
        cat("  Min   :", paste(sprintf("%s=%.2f", names(mins), mins), collapse = " | "), "\n")
        cat("  Median:", paste(sprintf("%s=%.2f", names(meds), meds), collapse = " | "), "\n")
        cat("  Max   :", paste(sprintf("%s=%.2f", names(maxs), maxs), collapse = " | "), "\n")
        
      } else {
        cat("Dynamic (local) tilt ranges: (unavailable)\n")
      }
    })
    
    output$bxp_plot <- renderPlot({
      
      bxp(
        bxp_obj_full,
        main = paste0("Boxplot of TAA Ranges; N=", nrow(W)),
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
        if (is.finite(drv_pos)) {
          points(drv_pos, info$target, pch = 23, bg = "blue", col = "white", cex = 1.2)
        }
        
        at_pos <- match(info$other, labs)
        if (any(is.na(at_pos)) || any(!is.finite(at_pos))) return()
        
        stats_m <- as.matrix(info$bxp$stats)
        if (ncol(stats_m) != length(at_pos)) return()
        
        fixed_halfwidth <- 0.14
        for (i in seq_along(at_pos)) {
          draw_mini_box(stats = stats_m[, i], x = at_pos[i], halfwidth = fixed_halfwidth)
        }
        
        legend(
          "topright",
          legend = c("Current factor SAA","Dynamic (local) tilt ranges"),
          pch = c(16, 15),
          col = c("red", "blue"),
          pt.cex = c(1.4, 1.2),
          bty = "n"
        )
        
      } else {
        
        legend(
          "topright",
          legend = c("Current factor SAA"),
          pch = 16,
          col = "red",
          pt.cex = 1.4,
          bty = "n"
        )
      }
    })
  })
}

# ============================================================
# App 2 (Module): Metric objective / green point
# ============================================================
app2UI <- function(id) {
  ns <- NS(id)
  
  metric_choices <- setNames(
    colnames(M),
    paste0(colnames(M), ifelse(metric_dir[colnames(M)] == 1, " (max)", " (min)"))
  )
  
  fluidPage(
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
        
        selectInput(
          ns("metric"), "Which objective should be prioritised?",
          choices = metric_choices,
          selected = colnames(M)[1]
        ),
        
        tags$hr(),
        
        numericInput(
          ns("mom"), "Momentum",
          value = NULL,
          min = W_min["Momentum"], max = W_max["Momentum"],
          step = STEP_2DP
        ),
        helpText(sprintf(
          "SAA: %.2f | Allowed range: [%.2f, %.2f]",
          SAA["Momentum"], W_min_disp["Momentum"], W_max_disp["Momentum"]
        )),
        
        numericInput(
          ns("qua"), "Quality",
          value = NULL,
          min = W_min["Quality"], max = W_max["Quality"],
          step = STEP_2DP
        ),
        helpText(sprintf(
          "SAA: %.2f | Allowed range: [%.2f, %.2f]",
          SAA["Quality"], W_min_disp["Quality"], W_max_disp["Quality"]
        )),
        
        numericInput(
          ns("val"), "Value",
          value = NULL,
          min = W_min["Value"], max = W_max["Value"],
          step = STEP_2DP
        ),
        helpText(sprintf(
          "SAA: %.2f | Allowed range: [%.2f, %.2f]",
          SAA["Value"], W_min_disp["Value"], W_max_disp["Value"]
        )),
        
        tags$hr(),
        actionButton(ns("reset"), "Reset (clear all)"),
        tags$hr(),
        verbatimTextOutput(ns("status"))
      ),
      mainPanel(
        plotOutput(ns("bxp_plot"), height = "520px")
      )
    )
  )
}

app2Server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # ---- Static boxplot object (App 2: min/0.25/0.5/0.75/max)
    stats_full <- sapply(seq_len(ncol(W)), function(j) {
      quantile(W[, j], probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE, names = FALSE)
    })
    
    bxp_obj_full <- list(
      stats = stats_full,
      n     = colSums(is.finite(W)),
      names = labs
    )
    
    # ---- Helpers (unchanged)
    K_NEIGH <- 20
    
    neighbour_idx_1d <- function(W, driver, target, k = 20) {
      j <- match(driver, colnames(W))
      ord <- order(abs(W[, j] - target))
      ord[seq_len(min(k, nrow(W)))]
    }
    
    adjust_from_neighbours_no_sum1 <- function(W, driver, target, k = 20) {
      nn <- neighbour_idx_1d(W, driver = driver, target = target, k = k)
      out <- numeric(3)
      names(out) <- colnames(W)
      out[driver] <- target
      other <- setdiff(colnames(W), driver)
      out[other] <- colMeans(W[nn, other, drop = FALSE], na.rm = TRUE)
      out
    }
    
    best_idx_by_metric <- function(M, metric, metric_dir, idx = NULL) {
      stopifnot(metric %in% colnames(M))
      stopifnot(all(colnames(M) %in% names(metric_dir)))
      
      if (is.null(idx)) idx <- seq_len(nrow(M))
      vals <- M[idx, metric]
      
      if (metric_dir[metric] == 1) {
        idx[which.max(vals)]
      } else {
        idx[which.min(vals)]
      }
    }
    
    # ---- State (unchanged)
    active_driver <- reactiveVal(NULL)
    adjusted <- reactiveVal(NULL)
    best_idx <- reactiveVal(NULL)
    
    observeEvent(input$reset, {
      active_driver(NULL)
      adjusted(NULL)
      best_idx(NULL)
      
      updateNumericInput(session, "mom", value = NULL)
      updateNumericInput(session, "qua", value = NULL)
      updateNumericInput(session, "val", value = NULL)
      
      shinyjs::enable(session$ns("mom"))
      shinyjs::enable(session$ns("qua"))
      shinyjs::enable(session$ns("val"))
    }, ignoreInit = TRUE)
    
    observeEvent(input$mom, {
      if (is.null(active_driver()) && !is.null(input$mom) && !is.na(input$mom)) {
        active_driver("Momentum")
        shinyjs::disable(session$ns("qua"))
        shinyjs::disable(session$ns("val"))
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$qua, {
      if (is.null(active_driver()) && !is.null(input$qua) && !is.na(input$qua)) {
        active_driver("Quality")
        shinyjs::disable(session$ns("mom"))
        shinyjs::disable(session$ns("val"))
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$val, {
      if (is.null(active_driver()) && !is.null(input$val) && !is.na(input$val)) {
        active_driver("Value")
        shinyjs::disable(session$ns("mom"))
        shinyjs::disable(session$ns("qua"))
      }
    }, ignoreInit = TRUE)
    
    observeEvent(list(input$mom, input$qua, input$val, input$method, input$metric), {
      
      drv <- active_driver()
      
      if (is.null(drv)) {
        adjusted(NULL)
        
        if (identical(input$method, "metric_opt")) {
          bi <- best_idx_by_metric(M, metric = input$metric, metric_dir = metric_dir, idx = NULL)
          best_idx(bi)
        } else {
          best_idx(NULL)
        }
        return()
      }
      
      target <- switch(
        drv,
        "Momentum" = input$mom,
        "Quality"  = input$qua,
        "Value"    = input$val
      )
      
      req(!is.null(target), !is.na(target))
      req(target >= W_min[drv], target <= W_max[drv])
      
      nn <- neighbour_idx_1d(W, driver = drv, target = target, k = K_NEIGH)
      
      if (identical(input$method, "nn_mean")) {
        p <- adjust_from_neighbours_no_sum1(W, driver = drv, target = target, k = K_NEIGH)
        adjusted(p)
        best_idx(NULL)
      } else {
        adjusted(NULL)
        bi <- best_idx_by_metric(M, metric = input$metric, metric_dir = metric_dir, idx = nn)
        best_idx(bi)
      }
      
    }, ignoreInit = TRUE)
    
    output$status <- renderPrint({
      drv <- active_driver()
      cat("Active driver:", if (is.null(drv)) "(none yet)" else drv, "\n")
      cat("Selection:", if (identical(input$method, "nn_mean")) "Nearby feasible tilt (blue)" else "Best portfolio (green)", "\n")
      
      dir_txt <- if (metric_dir[input$metric] == 1) "Maximise" else "Minimise"
      cat("Chosen objective:", input$metric, " [", dir_txt, "]\n\n", sep = "")
      
      cat("Inputs (what you typed; other boxes remain unchanged):\n")
      print(c(Momentum = input$mom, Quality = input$qua, Value = input$val))
      
      p <- adjusted()
      bi <- best_idx()
      
      cat("\nNearby feasible tilt (blue point):\n")
      if (is.null(p)) cat("(none)\n") else print(round(p, 2))
      
      cat("\nBest portfolio for chosen objective (green point):\n")
      if (is.null(bi)) {
        cat("(none)\n")
      } else {
        cat("Portfolio weights:", paste0(colnames(W), "=", round(W[bi, ], 2), collapse = ", "), "\n")
        cat("Objective value:", round(M[bi, input$metric], 6), "\n")
      }
    })
    
    output$bxp_plot <- renderPlot({
      
      bxp(
        bxp_obj_full,
        main = paste0("Boxplot of Factor Ranges; N=", nrow(W)),
        boxfill = "grey85",
        border  = "black"
      )
      
      points(1:3, SAA[1:3], pch = 21, bg = "red",  col = "red",   cex = 1.65, lwd = 2)
      points(1:3, SAA[1:3], pch = 21, bg = "red",  col = "white", cex = 1.30, lwd = 2)
      
      leg_labels <- c("Current portfolio position")
      leg_cols   <- c("red")
      leg_pch    <- c(16)
      leg_cex    <- c(1.4)
      
      p <- adjusted()
      if (!is.null(p)) {
        points(1:3, p[1:3], pch = 21, bg = "blue", col = "blue",  cex = 1.35, lwd = 2)
        points(1:3, p[1:3], pch = 21, bg = "blue", col = "white", cex = 1.05, lwd = 2)
        
        leg_labels <- c(leg_labels, "Nearby feasible tilt")
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
      
      legend(
        "topright",
        legend = leg_labels,
        pch = leg_pch,
        col = leg_cols,
        pt.cex = leg_cex,
        bty = "n"
      )
    })
  })
}

# ============================================================
# Main "two-tab website" wrapper (no behavioural changes)
# ============================================================
ui <- navbarPage(
  title = "Factor Portfolio Tools",
  useShinyjs(),  # one global shinyjs init is enough for both tabs/modules
  
  tabPanel("Tool 1: Local Tilt Ranges", app1UI("tool1")),
  tabPanel("Tool 2: Objective Optimiser", app2UI("tool2"))
)

server <- function(input, output, session) {
  app1Server("tool1")
  app2Server("tool2")
}

shinyApp(ui, server)
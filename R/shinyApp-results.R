.results_plot_ui <- function() {
  renderUI({
    tabBox(width = NULL, id = "input_tabs",
      tabPanel("Feature Explorer",
        plotlyOutput("plot"),
        shinyjs::disabled(
          actionBttn('update_annot', "Use selected molecule", style="bordered")
        ),
        br(),
        DTOutput('feature_tbl')
      ),
      tabPanel("Results table",
        downloadButton("tbldld","Download all results"),
        DTOutput('tbl')
      )
    )
  })
}
.extend_dt_btn <- function(name) {
  list(extend = name, filename = 'lipID_annotated',
    title = "Exported from lipID", exportOptions = list(
      modifier = list(page = "all")
    ))
}
.render_results_dt <- function(tbl, ...) {
  tbl %>% select(-n_and, -n_or, -n_and_true, -n_or_true, -and_cols, -or_cols) %>%
    datatable(extensions = 'Buttons', rownames = FALSE, ...) %>%
    formatSignif(sapply(., is.numeric), digits = 3)
}

.get_feature_tbl <- function(tbl, datum) {
  if (is.null(datum)) {
    datum <- list(x=0, y=0)
  }

  tbl %>% rename(mz=1, rt=2) %>%
    filter(mz == as.numeric(datum$y), rt == as.numeric(datum$x)) %>%
    ungroup() %>%
    select(-mz, -rt) %>%
    arrange(-best_match, -partial_match)
}

.select_row_as_best <- function(tbl, row) {
  tbl$best_match = ifelse(
    tbl$mz == row$mz & tbl$rt == row$rt,
    # For features, modify best_match to selection, leave others untouched.
    tbl$ms2_file == row$ms2_file & tbl$precursor == row$precursor &
      tbl$ms2_rt == row$ms2_rt & tbl$precursor == row$precursor &
      tbl$name == row$name,
    tbl$best_match
  )
  tbl
}
.update_tbl_data <- function(data, id) {
  data = data %>%
    select(-n_and, -n_or, -n_and_true, -n_or_true, -and_cols, -or_cols) %>%
    as.data.frame()
  replaceData(dataTableProxy(id), data, rownames = FALSE)
  return(data)
}

.update_feature_tbl_data <- function(data, id) {
  data <- .update_tbl_data(data, id)
  selectRows(dataTableProxy(id), which(data$best_match))
  return(data)
}

.render_results <- function(input, output, results) {
  showNotification("Rendering results", duration = 10)
  output$result <- .results_plot_ui()

  # Initial render using non-reactive results
  output$tbl <- renderDT(.render_results_dt(
    results,
    callback = JS("$('div.dwnld').append($('#tbldld'));"),
    options = list(dom = 'B<"dwnld">frtip', buttons=list())))

  output$feature_tbl <- renderDT({
    .get_feature_tbl(results, NULL) %>% .render_results_dt(
      selection = list(mode = "single", selected = which(.$best_match)),
      options = list(dom = 'Bfrtip', buttons=list('copy'))
    )
  }, server = TRUE)

  # Create reactives for results and feature tables
  ## These will be used to update / replace the data
  ## for rendered components (using JS, rather than re-rendering).
  tbl <- reactiveVal(results)
  output$plot <- renderPlotly(plot_features(tbl()) %>% event_register("plotly_click"))
  output$tbldld <- downloadHandler(
    filename = function() {
      paste("lipID_annotated-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write_csv(isolate(tbl()), file)
    }
  )

  selected_feature <-  reactiveVal()
  observe({
    .get_feature_tbl(tbl(), event_data("plotly_click")) %>%
      .update_feature_tbl_data('feature_tbl') %>%
      selected_feature()
  })

  observeEvent(input$feature_tbl_rows_selected, {
    req(input$feature_tbl_rows_selected)
    if (!is.null(selected_feature()) &&
        length(input$feature_tbl_rows_selected) == 1 &&
        !selected_feature()$best_match[[input$feature_tbl_rows_selected]]
        ) {
      shinyjs::enable('update_annot')
      shinyjs::addCssClass('update_annot', 'bttn-primary')
    } else {
      shinyjs::disable('update_annot')
      shinyjs::removeCssClass('update_annot', 'bttn-primary')
    }
  })

  observeEvent(input$update_annot, {
    selected_row <- selected_feature()[input$feature_tbl_rows_selected, ]
    mzrt <- isolate(event_data("plotly_click"))
    selected_row$rt = mzrt$x
    selected_row$mz = mzrt$y
    tbl(.select_row_as_best(tbl(), selected_row))
    .update_tbl_data(tbl(), 'tbl')
  })
}

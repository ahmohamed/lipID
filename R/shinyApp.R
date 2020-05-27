options(DT.options = list(dom = 'Bfrtip', buttons=c("copy", "csv", "excel")))
gobttn <- function(id, label) {
  actionBttn(id, loadingLabel(id, label),
    style="bordered", color = "primary", block = TRUE
  )
}

loadingLabel <- function(id, label) {
  span(
    label,
    hidden(span(
      icon("spinner", class = "btn-loading-indicator fa-spin"),
      id = paste0(id, "label")
    ))
  )
}

withLoading <- function(session, id, expr) {
  ns <- session$ns
  shinyjs::disable(id)
  shinyjs::show(paste0(id, "label"))
  # updateActionButton(session, id, icon = icon("spinner", class = "btn-loading-indicator fa-spin"))
  on.exit({
    shinyjs::hide(paste0(id, "label"))
    shinyjs::enable(id)
    # updateActionButton(session, 'dsbttn', icon = icon(FALSE))
  })
  tryCatch({
    value <- expr
    # Sys.sleep(5)
    return(value)
  })
}

#' @importFrom purrr quietly safely
quietly <- function(.f) {
  fun <- .f %>% purrr::quietly() %>% purrr::safely()
  function(...) {
    res <- fun(...)

    if(!is.null(res$error)) {  # safely output
      print("erros")
      if( !"shiny.silent.error" %in% class(res$error) ) {
        shinyalert(text=res$error$message, type="error")
        return(res$result)
      }
    }
    res <- res$result # quiely output
    if(!is.null(res$warnings) && length(res$warnings) > 0) {
      lapply(unique(res$warnings), showNotification, duration = 10, type="warning")
    }
    return(res$result)
  }
}

get_ui <- function(init_libs) {
  init_libs <- get_libs()
  ui <- tagList(
    useShinyjs(), useShinyalert(),
    singleton(tags$head(
      tags$style(
        ".datatables.html-widget{overflow-x:auto;}",
        ".pg-loading-screen {background-color: #17a2b8 !important;}")
    )),
    fluidPage(

      # Application title
      titlePanel("lipID: Fast Lipid Identification from MS/MS spectra"),

      sidebarLayout(
        sidebarPanel(
          fileInput('ms2_file', "Upload MS2 files", multiple = TRUE,
            accept = c(".ms2")
          ),
          numericInput('ppm_tol', 'MS2 matching tolerance (PPM)', 30, 0, 1000),
          numericInput('intensity_cutoff', 'MS2 fragment intensity cutoff', 100, 0, 100000),
          radioGroupButtons(
            inputId = "libmode",
            label = "Polarity",
            choices = c("Pos", "Neg"),
            checkIcon = list(yes = icon("ok", lib = "glyphicon"))
          ),
          pickerInput("liblist", "Lipid Classes",
            choices = init_libs$file, multiple = TRUE, selected = init_libs$file,
            options = list(`actions-box` = TRUE)),
          uiOutput('match_ms2_opts'),

          fileInput('features_file', "Upload features table",
            accept = c("txt/csv", "text/comma-separated-values,text/plain", ".csv")
          ),
          numericInput('mz_window', 'M/Z widnow (Daltons) for merging MS2 annotations with feature table', 1, 0, 10, 0.01),
          numericInput('rt_window', 'RT widnow (minutes) for merging MS2 annotations with feature table', 1, 0, 10, 0.1),
          gobttn('go', 'GIMME IDs!')
        ),
        # Show a plot of the generated distribution
        mainPanel(
          DTOutput('tbl')
        )
      )
    ))
}

server <- function(input, output, session) {
  libs <- get_libs()
  observeEvent(input$libmode, {
    print(input$libmode)
    libs <- get_libs(input$libmode)
    updatePickerInput(session, 'liblist',
      choices = libs$file, selected = libs$file)
  })

  observeEvent(input$go, {
    req(input$go > 0)
    print("go")
    if(is.null(input$ms2_file)){
      shinyalert(text="No input MS2 files", type="error")
      return()
    }
    if(is.null(input$liblist)){
      shinyalert(text="Please select one or more lipid classes", type="error")
      return()
    }
    if(!is.null(input$features_file)){
      features <- quietly(readr::read_csv)(input$features_file$datapath)
      f_copy <- quietly(.check_features_df)(features)
      req(f_copy)
    }

    withLoading(session, "go", {
      showNotification("Reading MS2 files", duration = 20)
      ms2_files <- setNames(input$ms2_file$datapath, input$ms2_file$name)
      print("ms2_files")
      print(ms2_files)
      ms2_data <- quietly(read_ms2)(ms2_files)

      showNotification("Matching against lipid libraries", duration = 30)
      selected_libs <- libs %>% filter(file %in% input$liblist)
      ms2_annotated <- match_ms2(ms2_data, selected_libs, input$ppm_tol, input$intensity_cutoff) %>%
        filter(partial_match >= 1)#input$partial_match_cutoff)

      if(!is.null(input$features_file)){
        showNotification("Merging MS2 annotations with feature table", duration = 10)
        tbl <- merge_ms2(features, ms2_annotated, input$mz_window, input$rt_window)
      } else {
        tbl <- ms2_annotated
      }

      showNotification("Rendering results", duration = 10)
      output$tbl <- renderDT({
        req(tbl)
        tbl %>% select(-n_and, -n_or, -n_and_true, -n_or_true, -and_cols, -or_cols) %>%
          datatable(extensions = 'Buttons', rownames = FALSE) %>%
          formatSignif(sapply(., is.numeric), digits = 3)
      }, server = FALSE)
    })
  })
}



#' Launch Shiny Web App  interface
#'
#' @return None
#' @import shiny
#' @importFrom shinyjs useShinyjs hidden
#' @importFrom shinyalert useShinyalert shinyalert
#' @importFrom shinyWidgets actionBttn radioGroupButtons pickerInput updatePickerInput
#' @importFrom DT datatable formatSignif renderDT DTOutput
#' @export
lipIDApp <- function(){
  shinyApp(ui = get_ui(), server = server)
}



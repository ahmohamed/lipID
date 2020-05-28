# options(DT.options = list(dom = 'Bfrtip', buttons=c("copy", "csv", "excel")))
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
how_to_use <- function() {

  tagList(
    tags$h1('How to use:'),
    tags$h3('1. Prepare you MS2 files:'),
    tags$p('Use MSConvert to format your raw MS2 data as *.ms2 files. ',
      'Upload all your ms2 files here.'),
    tags$h3('2. Select your libraries:'),
    tags$p('Choose the appropriate polarity and lipid classes that you want to search.'),
    tags$h3('3. Prepare your feature table:'),
    tags$p('Feature table should be a CSV file, with the first and second columns',
      ' corresponding to m/z and RT, similar to the example below. This is optional.'),
    tags$img(
      id = "ftable",
      src = "feature_table.png"
    ),
    tags$h3('4. Start matching!')
  )
}

get_ui <- function(init_libs) {
  init_libs <- get_libs()
  ui <- tagList(
    useShinyjs(), useShinyalert(),
    singleton(tags$head(
      tags$style(
        ".datatables.html-widget{overflow-x:auto;}",
        ".pg-loading-screen {background-color: #17a2b8 !important;}",
        "#ftable {max-width: 500px; width: 100%; height: auto}"
      ),
      tags$head(tags$title("lipID: Fast lipid ID from MS/MS spectra"))
    )),
    dashboardPage(skin = "black",

      # Application title
      dashboardHeader(title=tags$span(
        tags$img(
          src = "logo.png",
          height='40'
        ),
        "lipID: Fast lipid ID from MS/MS spectra"
        ),
        titleWidth = '100%'
      ),
      dashboardSidebar(disable = TRUE),
      dashboardBody(fluidRow(
        box(status = "primary", width = 4,
          fileInput('ms2_file', "Upload MS2 files", multiple = TRUE,
            accept = c(".ms2")
          ),
          numericInput('ppm_tol', 'MS2 matching tolerance (PPM)', 30, 0, 1000),
          numericInput('intensity_cutoff', 'MS2 fragment intensity cutoff', 100, 0, 100000),
          sliderInput('partial_match_cutoff', 'Partial matching cutoff', 0, 100, 100, 1, post = "%"),

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
          sliderInput(
            'mz_window',
            'Quadropole isolation window (Daltons)',
            0, 10, 1, 0.01, post=" Da"),
          sliderInput('rt_window',
            'RT window (minutes) for merging MS2 annotations with feature table',
            0, 10, 1, 0.1, post=" min"),
          gobttn('go', 'GIMME IDs!')
        ), #end box side
        box(status = "danger", width = 8,
          uiOutput('result')
        )
      )) #dashboard body
    ) #page
  )# taglist
}

server <- function(input, output, session) {
  output$result <- renderUI({
    how_to_use()
  })

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
      output$result <- renderUI(DTOutput('tbl'))

      showNotification("Reading MS2 files", duration = 20)
      ms2_files <- setNames(input$ms2_file$datapath, input$ms2_file$name)
      print("ms2_files")
      print(ms2_files)
      ms2_data <- quietly(read_ms2)(ms2_files)

      showNotification("Matching against lipid libraries", duration = 30)
      selected_libs <- libs %>% filter(file %in% input$liblist)
      ms2_annotated <- match_ms2(ms2_data, selected_libs, input$ppm_tol, input$intensity_cutoff) %>%
        filter(partial_match >= (input$partial_match_cutoff/100))

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
          datatable(extensions = 'Buttons', rownames = FALSE,
            options = list(
              dom = 'Bfrtip',
              buttons=list('copy', list(
                extend = 'collection',
                buttons = list(
                  list(extend = 'csv', filename = 'lipID_annotated'),
                  list(extend = 'excel', filename = 'lipID_annotated', title = "Exported from lipID"),
                  list(extend = 'pdf', filename = 'lipID_annotated', title = "Exported from lipID")
                ),
                text = 'Download'
              ))
          )) %>%
          formatSignif(sapply(., is.numeric), digits = 3)
      }, server = FALSE)
    })
  })
}

.launchApp <- function(){
  shinyApp(ui = get_ui(), server = server)
}

#' Launch Shiny Web App  interface
#'
#' @return None
#' @import shiny
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar dashboardBody box
#' @importFrom shinyjs useShinyjs hidden
#' @importFrom shinyalert useShinyalert shinyalert
#' @importFrom shinyWidgets actionBttn radioGroupButtons pickerInput updatePickerInput
#' @importFrom DT datatable formatSignif renderDT DTOutput
#' @export
lipIDApp <- function(){
  appDir <- system.file("shinyApp", package = "lipID")
  runApp(appDir = appDir)
}



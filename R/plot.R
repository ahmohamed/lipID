#' Interactively plot identified features
#'
#' @param results Annotated features from `[lipID]` or `[merge_ms2]`.
#' @param color_by Color points by whether they "annotated",
#' how they are assigned: "assigement", or their "class_name".
#'
#' @import plotly
#'
#' @return interactive plot
#' @export
plot_features <- function(results,
  color_by = c("class_name", "annotated", "assigement")) {
  color_by = match.arg(color_by)
  colors <- .get_colors(results[[color_by]])

  if (colnames(results)[[1]] == "ms2_file") { # MS2 only matching
    results <- results %>% mutate(mz=precursor, rt=ms2_rt) %>%
      select(mz, rt, everything())
  }

  results = results %>%
    group_by_at(vars(1,2)) %>% arrange(-best_match) %>%
    summarise_all(first) %>%
    mutate(
      annotated = ifelse(is.na(name) | !best_match, "Not annotated", "Annotated"),
      class_name = ifelse(
        is.na(class_name),
        "Not annotated", sub("(.{40}).*$", "\\1...", class_name))
    )

  if (!"assigement" %in% colnames(results)) {
    results$assigement = ifelse(results$annotated == "Not annotated", "Not Assigned", "Auto Assigned")
  }
  # TODO: Allow MS2_file res
  results <- rename(results, mz=1, rt=2)
  results$color_by <- results[[color_by]]

  plot_ly(size=I(35), hoverinfo="text", colors = colors) %>%
    add_markers(data=highlight_key(results %>% filter(annotated == "Not annotated"), ~mz+rt),
      x= ~rt, y= ~mz, color=I("grey"), text = "Not Annotated", name="Not Annotated", mode="markers", visible="legendonly") %>%
    add_markers(data = highlight_key(results %>% filter(annotated != "Not annotated"), ~mz+rt),
      x= ~rt, y= ~mz, color= ~color_by, text= ~name, mode="markers") %>%
    highlight(
      off='plotly_doubleclick', color="red",
      selected = attrs_selected(showlegend=FALSE)) %>%
    layout(dragmode = "zoom",
      xaxis= list(title = "Retention Time (min)"),
      yaxis= list(title = "m/z")) %>%
    config(doubleClick=FALSE, displayModeBar=TRUE, displaylogo=FALSE,
      showTips=FALSE,
      modeBarButtonsToRemove = c("lasso2d", "select2d", "resetScale2d",
        "hoverCompareCartesian", "hoverClosestCartesian", "toggleSpikelines"))
}


.get_colors <- function(values) {
  values <- na.omit(unique(values))
  cols <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(values))
  names(cols) <- sub("(.{40}).*$", "\\1...", values)
  cols
}


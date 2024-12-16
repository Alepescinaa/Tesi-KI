power_categorical <- function(significant_covs, data, scheme){
  data <- data[[scheme-1]]

  df_merged <- data %>%
    left_join(significant_covs, by = c( "transition"))
  df_merged <- df_merged %>%
    arrange(transition)
  
  cell_color <- function(value, condition) {
    if (condition) {
      if (value > 80) return("green")
      if (value > 50) return("yellow")
      return("red")
    }
    return("white")
  }
  
  styled_table <- df_merged %>%
    rowwise() %>%
    mutate(
      style_cov1.x = cell_color(cov1.x, cov1.y == 1),
      style_cov2.x = cell_color(cov2.x, cov2.y == 1),
      style_cov3.x = cell_color(cov3.x, cov3.y == 1)
    ) %>%
    ungroup()
  
  
  styled_table <- styled_table %>%
    select(model, transition, style_cov1.x, style_cov2.x, style_cov3.x) %>%
    mutate(cov1=style_cov1.x, cov2=style_cov2.x, cov3= style_cov3.x, style_cov1.x=NULL, style_cov2.x=NULL, style_cov3.x=NULL )
  
  styled_table <- styled_table %>%
    mutate(
      cov1 = cell_spec("",  background = cov1),
      cov2 = cell_spec("",  background = cov2),
      cov3 = cell_spec("",  background = cov3)
    )
  

  kbl(styled_table, escape = FALSE) %>%
    kable_styling("striped", full_width = FALSE)%>%
    add_header_above(c(" "=2, "Power" = 3)) %>%
    footnote(
      general = "Red : Power < 50%, Yellow : 50% < Power < 80%, Green : Power > 80%",
      general_title = "Legend:",
      footnote_as_chunk = TRUE
    )
  

  
  
  
}
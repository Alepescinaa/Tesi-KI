type_1_error <- function(significant_covs, data, scheme){
  data <- data[[scheme-1]]
  
  df_merged <- data %>%
    left_join(significant_covs, by = c( "transition"))
  df_merged <- df_merged %>%
    arrange(transition)
  
  cell_color <- function(value, condition) {
    if (condition) {
      if (value > 5) return("yellow")
      if (value > 10) return("red")
      return("green")
    }
    return("white")
  }
  
  styled_table <- df_merged %>%
    rowwise() %>%
    mutate(
      style_cov1.x = cell_color(cov1.x, cov1.y == 0),
      style_cov2.x = cell_color(cov2.x, cov2.y == 0),
      style_cov3.x = cell_color(cov3.x, cov3.y == 0)
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
    add_header_above(c(" "=2, "Error type I" = 3)) %>%
    footnote(
      general = "Green : Error_type_I < 5%, Yellow : 5% < Error_type_I < 10%, Green : Error_type_I > 10%",
      general_title = "Legend:",
      footnote_as_chunk = TRUE
    )
  
  
  
  
  
}
library(stringr)

# A bootstrap panel 
bsPanel <- function(..., heading = NULL, footer = NULL, type = "default") {

  if (is.null(heading)) {
    heading = ""
  }
  else {
    heading = div(class="panel-heading", heading)
  }
  
  if (is.null(footer)) {
    footer = ""
  }
  else {
    footer = div(class = "panel-footer", footer)
  }

  div(class = str_c("panel", str_c("panel-", type, sep = ""), sep = " "),
      heading,
      div(class = "panel-body", ...),
      footer
  )
}


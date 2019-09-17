library(shinythemes)
library(shinyjs)
library(shinycssloaders)
library(GenomicFeatures)
library(Gviz)
library(tidyverse)
library(stringr)
library(shiny)
library(DT)
library(here)
library(stringdist)

## options for knitting/rendering rmarkdown chunks
knitr::opts_chunk$set(echo = FALSE, comment = NA, cache = FALSE,
                      message = FALSE, warning = FALSE)

## function to render .md files to html
inclMD <- function(path)
  markdown::markdownToHTML(path, fragment.only = TRUE, options = "", stylesheet = "")

## function to render .Rmd files to html - does not embed image or add css
inclRmd <- function(path, r_env = parent.frame()) {
  paste(readLines(path, warn = FALSE), collapse = '\n') %>%
    knitr::knit2html(text = ., fragment.only = TRUE, envir = r_env,  options = "",
                     stylesheet = "") %>%
    gsub("&lt;!--/html_preserve--&gt;","",.) %>%  ## knitr adds this
    gsub("&lt;!--html_preserve--&gt;","",.) %>%   ## knitr adds this
    HTML
}


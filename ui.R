source("helpPage.R")
source("wcWidgets.R")

tags <- shiny::tags

tagList(

useShinyjs(),
  
navbarPage(
  id = "navbar_page",
  title = "WheatCrispr",
  theme = shinytheme("flatly"),
  tabPanel(
    title = "Home", 
    div(
      class = "container",
      fluidRow(
        column(
          width = 12,
          bsPanel(
            actionButton( "do_goto_grna", list("Search for gRNAs", HTML("&#8608;"))),
            h2("Introduction"),
            p(
              "WheatCrispr is an interactive tool for selecting CRISPR single
               guide RNAs (sgRNAs) in bread wheat. sgRNAs are scored according
               to predicted on-target activity and off-target potential using
               models devised by", em("Doench et. al., Nat. Biotech. 2016")
            ),
            p(
              "WheatCrispr uses the Doench on-target and off-target scores,
               along with information about the location of off-target hits
               (coding, intronic, intergenic, etc.) to produce a single
               overall score for each sgRNA.  By default, detailed information
               is displayed for the ten highest scoring sgRNAs to facilitate
               rapid identification of the most likely candidate sequences.
               An interactive interface allows the user to browse all other
               sgRNAs if desired. "
            ),
            p(
              "Polyploidy presents a unique challenge to sgRNA selection as
              the high similarity between homoeologues results in a greatly
               increased chance of sgRNAs targetting multiple homoeologues.
               Depending on the nature of the experiment this may be desired
               or it may not.  WheatCrispr supports two overall scoring
               schemes, one that favors sgRNAs that target only the queried
               gene by treating off-target hits in homoeologues as it would
               any other gene, and a second that targets all homoeologous
               genes by rewarding off-target activity in homoeologues."
            ),
            h2("Contact"),
            p(
              "For questions or comments contact us at Sateesh.Kagale(at)nrc.gc.ca and Dustin.Cram(at)nrc.gc.ca"
            )
          )
        )
      )
    )
  ),
  tabPanel(
    title = "Find gRNAs",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        radioButtons(
          inputId = "input_method",
          label = "Select an input method",
          choices = c("Gene Name" = "gene", "Paste a Sequence" = "sequence")
        ),
        conditionalPanel(
          "input.input_method == 'gene'",
	  textInput( inputId = "gene_id", label = "Gene Name"
          ),
          p(
            class = "form-text text-muted",
            "Enter a IWGSC v1 gene name.  For example: TraesCS6B02G093900"
          ),
          radioButtons(
            inputId = "ontarget_class",
            label = "Ontarget Set",
            choices = c("Coding" = "coding", "Promoters" = "promoter")
          ),
          p(
            class = "form-text text-muted", 
            "Select the grna target: the coding regions of the gene, or the
	     promoter region, defined as 2kbp upstream of the gene."
          ),
          checkboxGroupInput(
            inputId = "do_target_homoeologues",
            label = "Target Homoeologues?",
            choices = ""
          ),
          p(
            class = "form-text text-muted",
            "By default, offtarget hits to homoeologues penalize the overall
	     score, as with any other gene.  Select this box to adjust the
	     scoring scheme to reward offtarget hits in homoeologues."
          )
        ),
        conditionalPanel( 
          'input.input_method == "sequence"',
          textAreaInput(
            inputId = "user_sequence",
            label = "Sequence",
            placeholder = "paste a raw DNA sequence.",
            rows = 10,
            width = "100%"
          ) 
        )
      ),
      mainPanel(
        width = 9,
        fluidRow(
          column(
            width = 12,
            uiOutput("matched_gene")
          )
        ),
        fluidRow(
          column(
            width = 12,
            bsPanel(
              heading = "gRNAs Table",
              DT::dataTableOutput("grnas_table") %>% withSpinner()
            )
          )
        ),
        fluidRow(
          column(
            width = 6,
            bsPanel(
              heading = "gRNAs Plot",
              plotOutput("grnas_plot", height = 540) %>% withSpinner()
            )
          ),
          column(
            width = 6,
            bsPanel(
              heading = "Gene Plot",
              plotOutput("gene_plot") %>% withSpinner()
	    )
          )
        )
      )
    )
  ),
  tabPanel(
    title = "Help",
    fluidRow(
      column(
        width = 12,
        div(
          class="container", 
          #style="padding-left: 0px; padding-right: 0px;",
          wcHelp()
        )
      )
    )
  )
),

tags$head(
  tags$link(rel = "stylesheet", type = "text/css", href = "WheatCrispr.css"),
  
  tags$script(type = "text/javascript",
              async = "",
              src = "https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML")

  )
)



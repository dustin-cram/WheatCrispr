
tags <- shiny::tags

wcHelp <- function() {
  tagList(
    
    fluidRow(h2("Inputs")),
    fluidRow(h3("Select an Input Method")),
    fluidRow(
      column(
        width = 6,
	wellPanel(
	  radioButtons(
	    "help_input_method",
	    label = "Input Method",
	    choices = list("Gene Name" = "gene_name", "Paste a Sequence" = "user_sequence")
	  )
	)
      ),
      column(
        width = 6,
	bsPanel(
	  type = "info",
	  "Choose to target a gene by providing its name, or to target
           arbitrary sequence by manually providing that sequence.
           Providing a gene name is preferred whenever possible.  Details
           are described below"
	)
      )
    ),
    fluidRow(h3("Select a Gene")),
    fluidRow(
      column(
        width = 6,
        wellPanel(
          selectizeInput(
            inputId = "help_gene_id",
            label = NULL,
            choices = c("TraesCS3D02G273600"),
            selected = "TraesCS3D02G273600"
          )
        )
      ),
      column(
        width = 6,
        bsPanel(
          type = "info",
          "Start by entering the name of any gene here.  This must be a
          IWGSC annotation v1.0 or v1.1 identifier, for example \"TraesCS3D02G273600\""
        )
      )
    ),
    
    fluidRow(h3("Paste a Sequence")),
    fluidRow(
      column(
        width = 6,
	wellPanel(
	  textAreaInput(
	    inputId = "help_user_sequence",
	    label = "Sequence",
	    placeholder = "Paste a raw DNA sequence",
	    rows = 10
	  )
	)
      ),
      column(
        width = 6,
	bsPanel(
	  type = "info",
	  "Paste a DNA sequence.  This must be just raw sequence, no fasta headers or other characters."
	)
      )
    ),


    fluidRow(h3("Select Ontarget Set")),
    fluidRow(
      column(
        width = 6,
        wellPanel(
          radioButtons(
            "help_ontarget_set",
            label = "Ontarget Set",
            choices = list("Coding" = "coding", "Promoters" = "promoter"),
            inline = TRUE
          )
        )
      ),
      column(
        width = 6,
        bsPanel(
          type = "info",
          "Two sets of gRNAs are available for each gene: those located in
          the coding regions of the gene, and those located in the promoter
          regions (defined here as the 2kbp immediately upstream of the gene.)"
        )
      )
    ),
    
    fluidRow(h3("Select homoeologue scoring system")),
    fluidRow(
      column(
        width = 6,
        wellPanel(
          checkboxInput(
            "do_target_homoeologues",
            label = "Target Homoeologues? ",
            value = FALSE
          )
        )
      ),
      column(
        width = 6,
        bsPanel(
          type = "info",
          p("If this box is left blank then the overall score will be
            calculated using a method that treats homoeologues like any other
            gene, preferring gRNAs that are unique to the selected
            gene only."
          ),
          p("If this box is checked then the scoring scheme is adjusted to
            reward offtarget hits in homoeologues while still penalizing all
            other offtarget hits. "
          )
        )
      )
    ),
    
    fluidRow(
      column(
        width = 6,
        hr(),
             
        h2("Activity scores"),
             
        p("A brief explanation of the on-target and off-target activity
               scores is provided here.  For more information about the RS2
               and CFD scores, please refer to", em("Doench et. al., Nat. Biotech. 2016")),
             
        h3("RS2 (on-target) score"),
               
        p("The RS2 (ruleset 2) scores measures the predicted cutting
          efficiency of the gRNA.  This is a function only of the gRNA
          sequence itself plus a small flanking region on either side.  The
          score ranges from 0 (no predicted activity) to 1.0 (maximum
          activity)."),
        
             
        h3("CFD (off-target) scores"),
               
        p("The CFD (cutting frequency determination) score measures the
          predicted cutting efficiency of an off-target sequence relative to
          the on-target sequence.  This score ranges from 0 (no predicted
          activity), to 1.0 (full activity, ie. equal activity to the
          on-target site).  An off-target site with identical sequence to the
          on-target site will therefore always have a score of 1.0."),
             
        h3("Overall score"),
               
        p("The overall score is a weighted average of the RS2 and maximum CFD
          scores.  This score is not described by Doench _et. al._, it is
          specific to WheatCrispr.  __Note__: this scoring function is not
          based on any empirical evidence that suggests it provides the best
          possible balance between on-target and off-target activity.  It is
          simply an intuitive estimate designed to help accelerate the process
          of finding effective gRNAs.  Users are strongly encouraged to
          consider the individual RS2 and CFD scores, and other factors,
          before selecting a gRNA.  The exact function used when not
          targetting all homoeologues (the default mode) is:"
        ),
        
        p(" `0.5(text{rs2}) + 1 - (0.5( 0.7(max(text{cfd_coding}, text{cfd_promoter})) + 0.2(max(text{cfd_other_genic})) + 0.1(max(text{cfd_intergenic})) )`"),
        
        p("and when targetting homoeologues is enabled:"),
         
        p("`0.33(text{rs2}) + (1 - (0.33( 0.7(max(text{cfd_coding}, text{cfd_promoter})) + 0.2(max(text{cfd_other_genic})) + 0.1(max(text{cfd_intergenic})) ))) + 0.34(mean(text{cfd_hmlgs}))`")
        
      )
    ),

    fluidRow(hr()),
    fluidRow(h2("Outputs")),
    fluidRow(h3("gRNAs Table")),
    
    fluidRow(
      column(
        width = 12,
        bsPanel(
          img(src="img/sgrna_table.png")
        ),
        bsPanel(
          type = "info",
          p("This table summarizes all gRNAs for the selected gene.  The
            columns are: "),
          tags$dl(
            tags$dt("Sequence"),
            tags$dd("The protospacer sequence"),
            tags$dt("Overall Score"),
            tags$dd("The WheatCrispr overall score, as defined above"),
            tags$dt("RS2"),
            tags$dd("The Doench RS2 score, as defined above"),
            tags$dt("coding"),
            tags$dd("The maximum offtarget CFD score in a coding region"),
            tags$dt("promoter"),
            tags$dd("The maximum offtarget CFD score in a promoter region (defined
               as 2kbp upstream of a gene)"),
            tags$dt("other genic"),
            tags$dd("The maximum offtarget CFD score in other genic regions
               (introns and UTR)"),
            tags$dt("intergenic"),
            tags$dd("The maxmimum offtarget CFD score in the intergenic regions")
          ),
          p("Note that when the \"Target Homoeologue\" option is enabled that
            the maximum CFD scores will be adjusted to show the maximum CFD
            score outside of a homoeologue.  That is, the table always shows
            the \"worst\" CFD score."
          ),
          p("Clicking on a row while open a new table below that displays the
            set of all offtarget hits for the selected gRNA."),
          p("The table can be sorted by any column")
        )
      )
     ),

    fluidRow(hr()),
    fluidRow(h3("gRNAs Plot")),
    fluidRow(
      column(
        width = 12,
        bsPanel(
          img(src="img/sgrna_plot.png")
        ),
        bsPanel(
          type = "info",
          p("This plot displays a visualization of the scores found in the
            table, plus the score for any offtargets hits to homoeologues.
            The blue-gray bars shows the RS2 (ontarget) score, the black
            points indicate the worst CFD scores for each region (coding,
            promoter, other genic, intergenic), and the green points, if any,
            show the CFD scores for homoeologues."
          )
        )
      )
    ),
    
    fluidRow(hr()),
    fluidRow(h3("Offtargets Table")),
    fluidRow(
      column(
        width = 12,
        bsPanel(
          img(src="img/offtarget_table.png")
        ),
        bsPanel(
          type = "info",
          p("The offtarget table displays all offtarget hits for the selected gRNA.  The columns are:"),
          tags$dl(
            tags$dt("protospacer"),
            tags$dd("Sequence of the 20bp protospacer region"),
            tags$dt("pam"),
            tags$dd("Sequence of the 3bp PAM region"),
            tags$dt("partition"),
            tags$dd("The genomic region, one of \"coding\", \"promoter\", \"other genic\" (UTRs and introns), and \"intergenic\""),
            tags$dt("mismatches"),
            tags$dd("Number of mismatches to the gRNA"),
            tags$dt("cfd"),
            tags$dd("CFD score of the gRNA - offtarget site pair"),
            tags$dt("gene"),
            tags$dd("The gene in which this offtarget hits occurs (applicable only to coding and promoter regions")
          )
        )
      )
    ),
    
    fluidRow(hr()),
    fluidRow(h3("Gene Plot")),
    
    fluidRow(
      width = 12,
      bsPanel(
        img(src="img/gene_plot.png")
      ),
      bsPanel(
        type = "info",
        p("The gene plot displays the physical location of the gRNAs against
          the gene models.  Each row in the gene model represents an isoform
          of the gene.  The gray lines are introns and yellow bars exons.  The
          thinner bars indicate UTRs."
        )
      )
    )
  )
}

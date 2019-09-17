]84;0;0coptions(ucscChromosomeNames = FALSE)
options(shiny.host = "0.0.0.0")

# TODO: these should be read in from a config file
datadir <- "data/"
conda_source_script <- "/usr/local/anaconda/etc/profile.d/conda.sh"
conda_ss_env <- str_c("CONDA_SOURCE_SCRIPT", conda_source_script, sep = "=")


# print a tbl summary to stderr
tbl_perr <- function(tbl, n = 10) {
  message(paste(capture.output({print(tbl, n = n, width = 480)}), collapse = "\n"))
}

# Get the TxDB object for the genes from GTF (or the cached sqliteDB)
get_tx <- function() {

  sqlite_db <- here(datadir, "ref/txdb.sqlite")
  if (file.exists(sqlite_db)) {
    chr_tx <- loadDb(sqlite_db)
  }
  else {
    gff_file <- here(datadir, "ref/iwgsc.tidy.gff3")
    chrom_sizes_file <- here(datadir, "ref/scaffolds.sizes")
    chrominfo <- get_chrominfo(chrom_sizes_file)
    chr_tx <- makeTxDbFromGFF(gff_file, chrominfo = chrominfo)
    cat(file = stderr(), "made the TxDB!")
    utrs5 <- fiveUTRsByTranscript(chr_tx)
    saveDb(chr_tx, file=sqlite_db)
  }  
  
  return(chr_tx)
}

# Get the chrominfo for our genome
get_chrominfo <- function(file) {
    col_types <- 'ci'
    col_names <- c("chrom", "length")
    chrom_df <- read_tsv(file, col_types = col_types, col_names = col_names)
    return(chrom_df)
}

# Get a mapping of each gene to its homoeologues
get_homoeologues <- function() {
  filename <- here(datadir, "ref/homoeologues.tsv")
  col_types <- 'cc--'
  col_names <- c("gene_name", "hom_name")
  raw_df <- read_tsv(filename, col_types = col_types, col_names = col_names)
  
  df <- raw_df %>% 
    filter(!is.na(hom_name)) %>%
    arrange(gene_name)
  
  return(df)
}

# Compute the overall score.
grna_overall_score <- function(rs2, cfd_coding, cfd_promoter, cfd_other_genic,
                                cfd_intergenic) {
  
  # This is a single value between 0 and 1 that attempts to describe the overall
  # goodness of the grna, considering the rs2 and cfds.  Currently this
  # function is something I've basically just made up. I have no real evidence
  # that this is a particularly good function although my subjective feeling
  # browsing my data is that it does a reasonable job.  I view this score as a
  # way to help the user find interesting grnas faster, not an authoritative
  # value of the "best" grna.

  # Values may be NA if there are no offtargets at all in a partition (at least
  # none within the upstream limits).  Treat these as 0.
  
  cfd_coding <- if_else(is.na(cfd_coding), 0, cfd_coding)
  cfd_promoter <- if_else(is.na(cfd_promoter), 0, cfd_promoter)
  cfd_other_genic <- if_else(is.na(cfd_other_genic), 0, cfd_other_genic)
  cfd_intergenic <- if_else(is.na(cfd_intergenic), 0, cfd_intergenic)
  
  cfd_coding_promoter <- pmax(cfd_coding, cfd_promoter)
  
  cfd_overall <- cfd_coding_promoter * 0.7 +
    cfd_other_genic * 0.2 +
    cfd_intergenic * 0.1
  score_overall <- (rs2 * 0.5) + ((1 - cfd_overall) * 0.5)
    
  return(score_overall)
}

# An alternative overall_score function that rewards high CFD scores against
# homoeologues, while still penalizing high-CFDs against non-homoeologues.
hmlg_grna_overall_score <- function(rs2, cfd_coding, cfd_promoter,
                                     cfd_other_genic, cfd_intergenic,
                                     hmlgs_cfd_mean) {
  # As with the non-homoeologue overall score, the actual function implemented
  # here is just something a made up.
  
  cfd_coding <- if_else(is.na(cfd_coding), 0, cfd_coding)
  cfd_promoter <- if_else(is.na(cfd_promoter), 0, cfd_promoter)
  cfd_other_genic <- if_else(is.na(cfd_other_genic), 0, cfd_other_genic)
  cfd_intergenic <- if_else(is.na(cfd_intergenic), 0, cfd_intergenic)
  
  cfd_coding_promoter <- pmax(cfd_coding, cfd_promoter)
  cfd_no_hmlg_overall <- cfd_coding_promoter * 0.7 +
    cfd_other_genic * 0.2 +
    cfd_intergenic * 0.1
  
  hmlgs_cfd_mean <- if_else( is.nan(hmlgs_cfd_mean), 0, hmlgs_cfd_mean)
  
  score_overall <- (rs2 * 0.33) + ((1 - cfd_no_hmlg_overall) * 0.33) + (hmlgs_cfd_mean * 0.34)
  
  return(score_overall)  
}

# The full table of all homoeologues
message("Get the homeologues...")
all_homoeologues <- get_homoeologues()

# A list of genes/transcripts as TxDB objects
message("Get the transcript database")
gene_tx <- get_tx()

# Define server logic
shinyServer(function(input, output, session) {
 
  # a shinyjs snippet to change some navbar styling
  addClass(class = "container", selector = "nav.navbar > .container-fluid")
  removeClass(class = "container-fluid", selector = "nav.navbar > .container")
  
  genes_dir <- reactive({
      here(str_c(datadir, req(input$ontarget_class), "/"))
  })
  
  # Read the input file for this query (if "gene" then the precomputed result
  # for this gene, or if "sequence", then the result of an on-the-fly search)
  grnas_df <- reactive({
    if ( input$input_method == "gene" ) {
      return( grnas_df_gene_input() )
    }
    else if ( input$input_method == "sequence" ) {
      return( grnas_df_sequence_input() )
    }
  })
  

  # A reactive expression to validate/tweak the gene name provided.
  # Specifically, it will make sure the value is a valid IWGSC identifier, and
  # convert IWGSC annotation v1.0 gene names to v1.1 names
  valid_gene_id <- reactive({
   req(input$gene_id)
   validate(
      need(str_detect(input$gene_id, "^\\s*TraesCS[1-7][ABD]0[12]G[0-9]{6}"), "The provided gene ID does not look like a valid IWGSC name")
    )

    # Convert annotation version 1.0 names (01G) to annotation version 1.1 names (02G), and
    # get rid of leading/trailing whitespace
    clean_input <-
        str_replace(input$gene_id, "(TraesCS[1-7][ABD]0)[12](G[0-9]{6})", "\\12\\2") %>%
	str_trim()
    return(clean_input)
  })



  # Read in the pre-computed result for a gene and then tweak the format to make
  # some downstream work easier. This tibble (both input and final output) have
  # one row per onmer,off_20,pam_type tuple.  The counts for the four partitions
  # are listed in four distinct columns.
  grnas_df_gene_input <- reactive({
    
    req(valid_gene_id())
    
    # Read in the file for this gene from disk
    col_names <- c("on_ref", "on_start", "on_strand", "on_30mer",
                   "on_rs2_score", "off_20mer", "pam_type", "cfd_score", 
                   "count_coding", "count_promoter", "count_other_genic", 
                   "count_intergenic", "genes_coding", "genes_promoter",
                   "genes_other_genic")
    col_types <- c('ciccdccdiiiiccc')
    filename <- str_c(genes_dir(), valid_gene_id(), ".tsv.gz")
    
    tbl <- read_tsv(filename, col_names = col_names, col_types = col_types)
    
    # Get the homoeologues for this gene
    hmlgs <- homoeologues()$hom_name
    
    if ( req(input$ontarget_class) == "coding") {
      count_column_name = "count_coding"
      count_column = quo(count_coding)
    }
    else {
      count_column_name = "count_promoter"
      count_column = quo(count_promoter)
    }

    tbl %>%
      # Add an on_20mer column to make comparisons to the off_20mer easier
      mutate(on_20mer = str_sub(on_30mer, 5, 24)) %>%
      # Decrement the coding count when the offmer is == onmer and the pam is
      # canonical because one of those must be self
      mutate(!! count_column_name := if_else(on_20mer == off_20mer & pam_type == "GG",
                                             !! count_column - 1L,
                                             !! count_column)) %>%
      # Collapse the comma-separated genes lists into list-columns
      mutate(genes_coding = str_split(genes_coding, ","),
             genes_promoter = str_split(genes_promoter, ","),
             genes_other_genic = str_split(genes_other_genic, ",")) %>%
      # Get the list of homoeologues
      mutate(hmlgs_coding = map(genes_coding, ~ .[. %in% hmlgs])) %>%
      mutate(hmlgs_promoter = map(genes_promoter, ~ .[. %in% hmlgs])) %>%
      # Now summarize those gene lists to indicate the number of hits to a
      # homoeologue
      mutate(num_hmlgs_coding = genes_coding %>%
               map(~ . %in% hmlgs) %>% map_int(sum)) %>%
      mutate(num_hmlgs_promoter = genes_promoter %>% 
               map(~ . %in% hmlgs) %>% map_int(sum)) %>%
      # Add columns for cfd per partition
      mutate(cfd_coding = if_else(count_coding > 0, cfd_score, 0.0)) %>%
      mutate(cfd_promoter = if_else(count_promoter > 0, cfd_score, 0.0)) %>%
      mutate(cfd_other_genic = if_else(count_other_genic > 0,
                                       cfd_score, 0.0)) %>%
      mutate(cfd_intergenic = if_else(count_intergenic > 0, cfd_score, 0.0)) %>%
      select(-cfd_score)
  })
  
  # Do an on-the-fly search for the user-provided sequence and read in the
  # output file of that search.  Then tweak that tibble to make some downstream
  # work easier.  The format of the tibble (both the input and output) is one
  # row per (onmer, off_20mer, pam_type) tuple.  The four count_{partition}
  # columns have a different meaning here than they do in the
  # grnas_df_gene_input() function above.  In that function they are actual
  # counts, here they just indicate presence/absence, so that the value is never
  # higher than 1.  This is super annoying but became necessary due to the
  # compute requirements for doing the on-the-fly search.
  grnas_df_sequence_input <- reactive({

    user_sequence <- req(input$user_sequence)
    #TODO: safety and sanity checking on this sequence!
    
    # Do the onthefly search   
    message("Do onthefly search...")
    onfly_results_tmpfile <- tempfile(pattern = "onfly")
    onfly_results <- system2(here("scripts/onthefly_wrapper.sh"), args = "4",
                           stdout = onfly_results_tmpfile, input = user_sequence,
                           env = conda_ss_env)
    message("Reading in onthefly results...")
    
    # Read in the results of that search from the tempfile it produces
    col_names <- c("on_ref", "on_start", "on_strand", "on_30mer",
                      "on_rs2_score", "off_20mer", "pam_type", "cfd_score",
                      "count_coding", "count_promoter", "count_other_genic",
                      "count_intergenic")
    col_types <- 'ciccdccdiiii'
    tbl <- read_tsv(onfly_results_tmpfile, col_names = col_names,
                    col_types = col_types)
    
    message(str_glue("num search results: ", nrow(tbl)))

    message("Adding on_20mer column")
    # Add on on_20mer column for convenience
    tbl <- tbl %>%
      # Add on on_20mer column to make comparisons to the off_20mer easier
      mutate(on_20mer = str_sub(on_30mer, 5, 24))
    
    # Now blast the user provided sequence against the wheat genes, to
    # guess what the "self" gene is.
    message("Do the blast search to guess self...")
    predicted_gene <- get_matching_gene_for_sequence()
    
    # The offtarget info that we read in includes a hit for self, ie. for the
    # actual onmer that we want to target.  If we don't exclude self then we
    # will always have a hit with CFD = 1.0 in the ontarget partition which will
    # have a large negative effect our scoring function.  For the gene input
    # version of this function it's relatively easy to exclude the self hit by
    # decrementing the count by 1 in the appropriate partition.  But here, we
    # don't actually have a true count, just a presence/absence indicator.  If
    # we were to simply switch the exact match count from one to zero then we
    # could have the opposite problem where we would not drop the scoring
    # function in the case where there do exist other, true, exact match
    # offmers.  We will get around this limitation by looking up the actual
    # count for the exact match offmer in a sqlite database.
    message("Do the offtarget exact search...")
    offtarget_db <- here("data/ref/otd_all.db")
    ot_tmpfile <- tempfile(pattern = "ot")
    input_30mers <- str_c(unique(tbl$on_30mer), collapse = "\n")
    system2(here("scripts/query_otd_db_exact_wrapper.sh"),
            args = offtarget_db, stdout = ot_tmpfile, input = input_30mers,
            env = conda_ss_env)
    
    message("Reading in exact matches search...")
    ot_col_names <- c('on_20mer', 'gene_name')
    ot_col_types <- 'cc'
    exact_matches_tbl <- read_tsv(ot_tmpfile, col_names = ot_col_names,
                                  col_types = ot_col_types)
    
    message("Exact matches results: ")
    
    # TODO: Here I'm assuming that a single exact match is self.  That isn't
    # necessarily true if the user has pasted a sequence that does not match the
    # reference sequence. 
    message("Summarize exact matches to counts")
    exact_counts_tbl <- exact_matches_tbl %>%
      group_by(on_20mer) %>%
      summarize(exact_count = n()) %>%
      ungroup()
    
    message("Join the exact counts tbl to the main tbl")
    # Drop the count from one to zero in cases where there is only one exact match
    tbl <- tbl %>%
      left_join(exact_counts_tbl, by = c("on_20mer")) %>%
      mutate(count_coding = if_else((on_20mer == off_20mer & exact_count == 1),
                                    0L,
                                    count_coding))

    message("Add columns for cfd per partition")
    tbl <- tbl %>%
      # Add columns for cfd per partition
      group_by(on_ref, on_start, on_strand) %>%
      mutate(cfd_coding = if_else(count_coding > 0, cfd_score, 0.0)) %>%
      mutate(cfd_promoter = if_else(count_promoter > 0, cfd_score, 0.0)) %>%
      mutate(cfd_other_genic = if_else(count_other_genic > 0,
                                       cfd_score, 0.0)) %>%
      mutate(cfd_intergenic = if_else(count_intergenic > 0, cfd_score, 0.0)) %>%
      select(-cfd_score) %>% 
      # Add empty genes_coding and genes_promoter columns. This is so that the tbl
      # has the same shape as the gene_input method
      add_column(genes_coding = list(character(0)),
                 genes_promoter = list(character(0)),
                 hmlgs_coding = list(character(0)),
                 hmlgs_promoter = list(character(0)))
    
    message("join table results: ")
    
    return(tbl)
    
  })
  
  # Take the grnas_df above, and summarize it to have one row per grna with
  # interesting offtarget matches summarized into new columns
  grnas_spread <- reactive({
    
    if (input$input_method == "gene" && length(input$do_target_homoeologues) && 
        nrow(homoeologues())) {
	message("Doing grnas_spread with homoeologues")
      return( grnas_spread_with_homoeologues() )
    }
    else {
      message("Doing grnas_spread with no homoeologues")
      return( grnas_spread_no_homoeologues() )
    }
  })

  grnas_spread_no_homoeologues <- reactive({
    
    grnas_df() %>%
      group_by(on_ref, on_start, on_strand) %>%
      # Summarize by grabbing the max cfd score
      summarize(on_rs2_score = first(on_rs2_score), on_30mer = first(on_30mer),
                coding = max(cfd_coding), promoter = max(cfd_promoter),
                other_genic = max(cfd_other_genic), intergenic = max(cfd_intergenic)) %>%
      ungroup() %>%
      # Add a column with an overall_score which is a weighted average of the
      # rs2 and the max_cfds
      mutate(overall_score = grna_overall_score(on_rs2_score, coding, promoter,
                                                other_genic, intergenic)) %>%
      # Add an explicit rank by descending overall_score
      mutate(rank = row_number(desc(overall_score))) %>%
      # Extract the protospacer sequence from the full 30mer sequence
      mutate(seq = str_sub(on_30mer, 5, 24),
             pam = str_sub(on_30mer, 25, 27)) %>%
      # Select the columns we care about
      select(on_ref, on_start, on_strand, seq, pam, overall_score, rank,
             on_rs2_score, coding, promoter, other_genic, intergenic) %>%
      # Arrange from best to worst (descending overall_score)
      arrange(desc(overall_score))
  })
  
  grnas_spread_with_homoeologues <- reactive({
    
    hmlgs <- homoeologues()$hom_name
    
    grnas_df() %>%
      group_by(on_ref, on_start, on_strand) %>%
      summarize(on_30mer = first(on_30mer),
                on_rs2_score = first(on_rs2_score),
                coding =
                  if_else(req(input$ontarget_class) == "coding",
                          max( cfd_coding[num_hmlgs_coding < count_coding] ),
                          max( cfd_coding )),
                promoter = 
                  if_else(req(input$ontarget_class) == "promoter",
                          max( cfd_promoter[num_hmlgs_promoter < count_promoter] ),
                          max( cfd_promoter )),
                other_genic = max(cfd_other_genic),
                intergenic = max(cfd_intergenic),
                hmlgs_cfd_mean =
                  if_else(input$ontarget_class == "coding",
                          map_dbl(hmlgs, function(h)   cfd_coding[map_lgl(  genes_coding, function(g) h %in% g)] %>% max(., 0) ) %>% mean(),
                          map_dbl(hmlgs, function(h) cfd_promoter[map_lgl(genes_promoter, function(g) h %in% g)] %>% max(., 0) ) %>% mean()
                  )
      ) %>%
      ungroup() %>%
      mutate(overall_score = 
               hmlg_grna_overall_score(on_rs2_score, coding, promoter,
                                       other_genic, intergenic,
                                       hmlgs_cfd_mean)) %>%
      mutate(rank = rank(dplyr::desc(overall_score))) %>%
      mutate(seq = str_sub(on_30mer, 5, 24),
             pam = str_sub(on_30mer, 25, 27)) %>%
      dplyr::select(on_ref, on_start, on_strand, seq, pam, overall_score, rank, 
             on_rs2_score, coding, promoter, other_genic, intergenic) %>%
      arrange(desc(overall_score))
  })
 
  # Reduce the grnas_df down to "interesting" hits.  That's the max cfd per
  # partition, plus the homoeologues.  Unlike the grnas_spread* functions, this
  # retains each hit as a separate row
  grnas_to_plot <- reactive({
    
    # Gather the cfds per partition, so we now have four rows per [on_20mer,
    # off_20mer, pam_type]
    gathered_df <- grnas_df() %>%
      dplyr::rename(coding = cfd_coding, promoter = cfd_promoter,
             other_genic = cfd_other_genic, intergenic = cfd_intergenic) %>%
      gather(key = "off_partition", value = "cfd_score", coding, promoter,
             other_genic, intergenic) %>%
      mutate(off_partition = factor(off_partition,
                                    levels = c("coding", "promoter",
                                               "other_genic", "intergenic")))
    
    # Add a num_hmlgs column.  If our input type is sequence (ie. user pasted
    # sequence), then we don't support any of the homoeologue-aware features, so
    # just set it to zero.  If the input type is gene, then this just normalizes
    # the ontarget classes to give a consistent column name
    if (input$input_method == "gene") {
      gathered_df <- gathered_df %>%
        mutate(num_hmlgs = if(input$ontarget_class == "coding")
          num_hmlgs_coding else num_hmlgs_promoter) %>%
        select(-num_hmlgs_coding, -num_hmlgs_promoter)
    }
    else { # input$input_method == "sequence"
      gathered_df <- gathered_df %>%
        mutate(num_hmlgs = 0)
    }
    
    # Get a new tibble that filters the input down to just the max cfds for each
    # onmer/partition.  
    max_cfds_df <- gathered_df %>%
      group_by(on_ref, on_start, on_strand, off_partition) %>%
      filter(row_number(desc(cfd_score)) == 1)
      
    # Get a new tibble that includes only hits to homoeologues
    hmlgs_df <- gathered_df %>%
      filter(num_hmlgs > 0) %>%
      filter(off_partition == input$ontarget_class)
    
    # Rowbind the max cfds and hmlgs.  Together these two tibbles contain our
    # "interesting" hits
    to_plot_df <- bind_rows(max_cfds_df, hmlgs_df) %>%
      # Join on the spread table to get the overall score
      left_join(grnas_spread(), by = c("on_ref", "on_start", "on_strand")) %>%
      mutate(is_homoeologue = num_hmlgs > 0) %>%
      select(on_ref, on_start, on_strand, seq, on_rs2_score = on_rs2_score.x,
             off_partition, cfd_score, is_homoeologue, overall_score, rank)
  })
  
  # Get the current "page" of hits to plot, based on the paging of the table
  grnas_to_plot_current <- reactive({
    current_rows <- req(input$grnas_table_rows_current)
    grnas_to_plot() %>%
      filter(rank %in% current_rows)
  })
  
  # Given the user-pasted sequence, returns the most similar gene, or NA if the
  # most similar doesn't meet some threshold
  get_matching_gene_for_sequence <- reactive({
    
    user_sequence <- input$user_sequence
    
    # Now we will try to figure out which gene, if any, the pasted sequence was
    # derived from.  This is important because if we don't know where each grna
    # comes from, then we are forced to interpret all (exact) occurrences of
    # that grna as an offtarget, which results in a large bogus penalty to the
    # overall score for each grna as the true ontarget site is counted as an
    # offtarget.
    
    # submit pasted sequence to an external program that will align
    # against the genes and report the best match
    blastdb = "data/ref/IWGSC_v1.1_cds.fasta"
    predicted_gene <- system2(here("scripts/blast_against_gene_wrapper.sh"),
                              args = blastdb, stdout = TRUE,
                              input = user_sequence,
                              env = conda_ss_env)

  })

  # The offtargets table for the selected grna
  offtargets <- reactive({

    selected_rows <- req(input$grnas_table_rows_selected)
    grna <- grnas_spread()[selected_rows,]
    
    offtargets_df <- grnas_df() %>%
      # Filter to include only the offtargets for this ontarget grna
      filter(on_ref == grna$on_ref, on_start == grna$on_start,
             on_strand == grna$on_strand) %>%
      ungroup()

    # Calculate mimatches between on_20mer and off_20mer
    offtargets_df <- offtargets_df %>%
       mutate(mismatches = as.integer(stringdist(on_20mer, off_20mer, method = "hamming")))
    
    # Get the protospacer and PAM sequences
    offtargets_df <- offtargets_df %>%
      dplyr::rename(protospacer = off_20mer, pam = pam_type)
    
    if (input$input_method == "gene") {
      offtargets_df <- offtargets_df %>%
        mutate(num_hmlgs = if (input$ontarget_class == "coding")
          num_hmlgs_coding else num_hmlgs_promoter)
    }
    else {
      offtargets_df <- offtargets_df %>%
        mutate(num_hmlgs = 0)
    }
    
    # smash the four cfd columns into one
    offtargets_df <- offtargets_df %>%
      mutate(cfd_score = pmax(cfd_coding, cfd_promoter, cfd_other_genic, cfd_intergenic)) %>%
      select(protospacer, pam, starts_with("count_"), num_hmlgs, cfd_score, genes_coding,
             genes_promoter, hmlgs_coding, hmlgs_promoter, mismatches) %>%
      arrange(desc(cfd_score))
  })
  
  # Returns a (single column) tibble listing the homoeologues for this gene
  homoeologues <- reactive({

    if (str_length(valid_gene_id())) {
      gene_id <- valid_gene_id()
    }
    else {
      gene_id <- ""
    }
    all_homoeologues %>%
      filter(gene_name == gene_id) %>%
      select(hom_name)
  })
  
  # Get the Gviz track for the grnas.
  grnas_track <- reactive({

    # Get the txstart because we want coordinates relative to the gene start, not chromosome
    # gene_tbl <- AnnotationDbi::select(gene_tx, 
    #                       columns = c("TXSTRAND", "TXSTART"), 
    #                       keys = req(valid_gene_id()),
    #                       keytype = "GENEID")
    # 
    # gene_start <- min(gene_tbl$TXSTART)
    
    #Get the grnas info (all grnas, or just those displayed?)
    current_rows <- req(input$grnas_table_rows_current)
    
    grnas <- req(grnas_spread()) %>%
      filter(rank %in% current_rows)
    track <- AnnotationTrack(chromosome = grnas$on_ref[1], 
                             start = grnas$on_start, width = 20,
                             strand = grnas$on_strand, group = grnas$rank, 
                             name = "gRNAs", groupAnnotation = "group",
                             just.group = "left", cex.group = 1.0)
    return(track)
  })
  
  # Get a simple gviz track for the case when the user pastes in an arbitrary
  # sequence.  This will just be a solid bar since have no exon/intron/utr, etc.
  # info.  But still somewhat useful because we can show where the grnas occur
  # relative to the user sequence.
  dummy_gene_model_track <- reactive({
    req(input$user_sequence)
    
    seqlen <- str_length(input$user_sequence)
    
    dummy_gene_df <- tribble(
      ~chromosome, ~start, ~end, ~width, ~strand, ~feature, ~gene, ~exon, ~transcript, ~symbol,
      "user_sequence", 1, seqlen, seqlen, "+", "gene", 1, 1, 1, "user_sequence"
    )
    
    gene_model_track <- GeneRegionTrack(dummy_gene_df, name="Gene Model")
  })
  
  # Get the Gviz track for the currently selected gene.
  gene_model_track <- reactive({
    
    if (input$input_method == "sequence") {
      return(dummy_gene_model_track())
    }
    
    # Note that while there are convenience methods in GViz to directly create
    # tracks from a TxDB object, I've found these to be extremely slow (~75sec.)
    # There are also helper methods to create Granges representing gene models
    # from the TxDB, but these aren't much faster.  So I'm doing the query with
    # AnnotationDbi::select manually, and wrangling the output into a dataframe
    # format accepted by GeneTrackRegion manually.  This is super annoying
    # because TxDB does not explicitly record UTR info.
    
    model_raw_df <-
      AnnotationDbi::select(gene_tx, 
                            columns = c("GENEID", "CDSID", "CDSSTART", "CDSEND",
                                        "TXSTRAND", "EXONID", "EXONSTART",
                                        "EXONEND", "TXID", "TXNAME", "TXCHROM",
                                        "TXSTART", "TXEND"), 
                            keys = req(valid_gene_id()),
                            keytype = "GENEID")

    model_df <- model_raw_df %>% 
      # Add 5' UTR (or NAs if there isn't any)      
      mutate(UTR5ID = ifelse(CDSSTART == EXONSTART, NA, EXONID),
             UTR5START = ifelse(CDSSTART == EXONSTART, NA, EXONSTART),
             UTR5END = ifelse(CDSSTART == EXONSTART, NA, CDSSTART - 1)) %>%
      # Add 3' UTR (or NAs if there isn't any)
      mutate(UTR3ID = ifelse(CDSEND == EXONEND, NA, EXONID), 
             UTR3START = ifelse(CDSEND == EXONEND, NA, CDSEND + 1), 
             UTR3END = ifelse(CDSEND == EXONEND, NA, EXONEND)) %>%
      # If this is a UTR only exon it's tough to tell if this is 5' or 3' UTR
      # (yes we can check if it aligns to transcript start/end but UTR only
      # exons aren't necessarily terminal), so I'm just arbitrarily calling it
      # 3'UTR.  For the purpsoses of plotting it doesn't matter if this is
      # correct
      mutate(UTR3ID = ifelse(is.na(CDSID), EXONID, UTR3ID),
             UTR3START = ifelse(is.na(CDSSTART), EXONSTART, UTR3START),
             UTR3END = ifelse(is.na(CDSEND), EXONEND, UTR3END))
    
    # Extract separate dataframes for each of CDS, UTR5, and UTR3 then rbind()
    # them together
    
    model_cds_df <- model_df %>%
      filter(!is.na(CDSID)) %>%
      mutate(feature = "protein_coding", width = CDSEND - CDSSTART + 1) %>%
      select(chromosome = TXCHROM, start = CDSSTART, end = CDSEND,
             width, strand = TXSTRAND, feature, gene = GENEID, exon = EXONID,
             transcript = TXID, symbol = TXNAME, txstart = TXSTART)
    
    model_utr5_df <- model_df %>%
      filter(!is.na(UTR5ID)) %>%
      mutate(feature = "utr", width = UTR5END - UTR5START + 1) %>%
      select(chromosome = TXCHROM, start = UTR5START, end = UTR5END,
             width, strand = TXSTRAND, feature, gene = GENEID, exon = EXONID,
             transcript = TXID, symbol = TXNAME, txstart = TXSTART)
    
    model_utr3_df <- model_df %>%
      filter(!is.na(UTR3ID)) %>%
      mutate(feature = "utr", width = UTR3END - UTR3START + 1) %>%
      select(chromosome = TXCHROM, start = UTR3START, end = UTR3END,
             width, strand = TXSTRAND, feature, gene = GENEID, exon = EXONID,
             transcript = TXID, symbol = TXNAME, txstart = TXSTART)
    
    model_tidy_df <- rbind(model_cds_df, model_utr5_df, model_utr3_df)
    
    # # Now let's convert the coordinates from chromosome space to gene space
    # # First, find the minimum TXSTART as the gene start
    # gene_start <- min(model_tidy_df$txstart)
    # model_tidy_df <- model_tidy_df %>% 
    #   mutate(start = start - txstart + 1,
    #          end = end - txstart + 1)
    
    gene_model_track <- GeneRegionTrack(model_tidy_df, name="Gene Model")
    

    return(gene_model_track)
  })
  
  #### Output Elements
  
  observeEvent(
    input$do_goto_grna,
    {
      updateNavbarPage(
        session,
        inputId = "navbar_page",
        selected = "Find gRNAs"
      )
    },
    priority = 1
  )
  
    observeEvent(
    input$grnas_table_rows_selected,
    {
      showModal(
        modalDialog(
            DT::dataTableOutput("offtargets_table") %>% withSpinner(),
	          size = 'l'
        )
      )
      
    }
  )
  

  # Output the grnas as a DT table
  output$grnas_table <- DT::renderDataTable({
    
    DT::datatable(
      grnas_spread() %>%
        select(on_ref, on_start, on_strand, seq, pam, overall_score,
               on_rs2_score, coding, promoter, other_genic, intergenic),
      extensions = 'Buttons',
      options = list(dom = 'Btip',
                     buttons = c('copy', 'csv', 'excel', 'print'),
                     columnDefs = list(
                       list(visible = FALSE, targets = c(1,2,3,5)),
                       list(defaultContent = "None", targets = c(8,9,10,11))
                     )),
      caption = "Click on any row to view all off-targets for that gRNA",
      selection = "single", 
      colnames = c("", "on_ref", "on_start", "on_strand", "sequence", "pam",
                   "overall score", "rs2", "coding", "promoter", "other genic", 
                   "intergenic")
    ) %>% 
      formatSignif( c("overall_score","on_rs2_score","coding","promoter",
                      "other_genic","intergenic"), 2) %>%
      formatStyle( "seq", fontFamily = "monospace" )
  }, server = FALSE)
  
  # Output a bar/scatter type plot showing a summary of the grnas
  output$grnas_plot <- renderPlot({

    bar_df <- grnas_to_plot_current() %>% group_by(rank) %>% dplyr::slice(1)
    ggplot(data = bar_df) +
      geom_bar(
        aes(
          x = factor(match(rank, req(input$grnas_table_rows_current))),
          y = on_rs2_score, fill = "blue"),
        color = "black", alpha = 0.05, stat = "identity", width = 0.5) +
      geom_point(
        data = grnas_to_plot_current(),
        aes(
          x = factor(match(rank, input$grnas_table_rows_current)),
          y = cfd_score, shape = off_partition, colour = is_homoeologue), 
        size = 4, alpha = 0.6, stroke = 1) +
      scale_shape_manual( values = c(16,17,0,3) ) + 
      scale_colour_manual( values = c("black", "#00CC00") ) +
      scale_fill_identity( guide = 'legend', labels = c("bars") ) +
      scale_x_discrete(labels = c(input$grnas_table_rows_current)) +
      scale_y_continuous(breaks = seq(0,1,0.2), minor_breaks = seq(0,1,0.1),
                         limits = c(0,1)) +
      labs(shape = "Offtarget Partition", colour = "Is homoeologue?",
           fill = "Ontarget (RS2) value",
           x = "Overall Score Rank", y = "RS2 / CFD Scores") +
      theme_minimal() +
      theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
            axis.text = element_text(size = 18),
            axis.title = element_text(size = 18),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 16),
            legend.position = "bottom",
            legend.direction = "vertical",
            legend.box = "horizontal"
            )
  })
  
#  Output the offtargets for a single grna
  output$offtargets_table <- DT::renderDataTable({

  ot_tbl <-  cbind('expand' = '&oplus;', offtargets());
    
    DT::datatable(
      ot_tbl %>%
        select(expand, protospacer, pam, cfd_score, mismatches, num_hmlgs, starts_with("count_"),
               genes_coding, genes_promoter, hmlgs_coding, hmlgs_promoter),
      extensions = 'Buttons',
      options = list(dom = 'Btip',
                     buttons = c('copy', 'csv', 'excel', 'print'),
                     columnDefs = list(
                       list(visible = FALSE, targets = c(0,11,12,13,14)),
                       list(orderable = FALSE, className = "genes_expander", targets = c(1))
                       )
                     ),
      selection = "none",
      filter = 'top',
      escape = c(1),
      colnames = c("", "protospacer", "pam", "cfd", "mismatches", "homoeologue hits",
                   "coding region hits", "promoter region hits", "other genic region hits",
                   "intergenic region hits", "coding genes", "promoter genes",
                   "coding homoeologues", "promoter homoeologues"),


      callback = JS("
var fmt = function(d) {
  var result = '<div style=\"background-color:#fff; padding: .5em;\">';

  /* For some reason, R vectors of length 1 are converted to just 
    strings instead of arrays of length 1 */

  genes_coding = d[11];
  genes_promoter = d[12];
  hmlgs_coding = d[13];
  hmlgs_promoter = d[14];

  if ( typeof(genes_coding) == 'string' ) { genes_coding = [genes_coding] }
  if ( typeof(hmlgs_coding) == 'string' ) { hmlgs_coding = [hmlgs_coding] }
  if ( typeof(genes_promoter) == 'string' ) { genes_promoter = [genes_promoter] }
  if ( typeof(hmlgs_promoter) == 'string' ) { hmlgs_promoter = [hmlgs_promoter] }

  if ( genes_coding != null && genes_coding.length > 0 ) {
    result += 'Coding Genes: ' + genes_coding.join(', ') + '<br />';
  }
  if ( hmlgs_coding != null && hmlgs_coding.length > 0 ) {
    result += 'Coding Homoeologues: ' + hmlgs_coding.join(', ') + '<br />';
  }
  if ( genes_promoter != null && genes_promoter.length > 0 ) {
    result += 'Promoter Genes: ' + genes_promoter.join(', ') + '<br />';
  }
  if ( hmlgs_promoter != null && hmlgs_promoter.length > 0 ) {
    result += 'Promoter Homoeologues: ' + hmlgs_promoter.join(', ');
  }
  result + '</div>';
  return result;
};
table.on('click', 'td.genes_expander', function() {
  var td = $(this), row = table.row(td.closest('tr'));
  if (row.child.isShown()) {
    row.child.hide();
  } else {
    row.child(fmt(row.data())).show();
  }
});
")

      
      
      
      ) %>%
      formatSignif( c("cfd_score"), 2) %>%
#      formatStyle( "off_gene", "is_homoeologue",
#                   backgroundColor = styleEqual( c(1,0), c("#00CC00","") ),
#                   color = styleEqual( c(1,0), c("white",""))) %>%
      formatStyle( "cfd_score",
                   backgroundColor = styleInterval(c(0.05,0.2,0.9999),
                                                   c("#FFFFFF", "#FFD0D0",
                                                     "#FF8080","#FF0000"))

      ) %>%
      formatStyle( "protospacer", fontFamily = "monospace" ) %>%
      formatStyle( "pam", fontFamily = "monospace" )
  })
  

  # Output the gene plot
  output$gene_plot <- renderPlot({

    plotTracks(c(GenomeAxisTrack(cex = 1.2), gene_model_track(), grnas_track()),
               stackHeight = 0.6, col.title = "black", cex.title = 1.2, shape = "arrow")
  }
  )
  
  output$matched_gene <- renderUI({
    req(input$input_method)
    if ( input$input_method == "sequence" ) {
      req(input$user_sequence)

      tagList(
        "The pasted sequence has highest similarity to: ",
        strong(get_matching_gene_for_sequence()),
        ".  Hits against this gene will not count as offtargets."
      )
    }
  })
  
  # Output the help page Rmd
#  output$help <- renderUI({inclRmd("WheatCrispr_help.Rmd")})
})

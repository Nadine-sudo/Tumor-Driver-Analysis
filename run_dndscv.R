##############################################################################
# Universal dNdScv Analysis Pipeline for Multiple Species
# Workflow: 1. Install Dependencies 2. Preprocess MAF Data 3. Build Species-specific RefCDS (if needed)
#           4. Run dNdScv Analysis 5. Filter Candidate Genes 6. Generate Visualizations
##############################################################################

# ==============================
# Step 1: Install & Load Dependencies
# ==============================
install_dependencies <- function() {
  # Core dependencies
  required_pkgs <- c("dplyr", "devtools", "BiocManager")
  for (pkg in required_pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if (pkg == "BiocManager") {
        install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
      } else {
        install.packages(pkg, dependencies = TRUE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
      }
      library(pkg, character.only = TRUE, quietly = TRUE)
    }
  }
  
  # Bioinformatics dependencies (required for dndscv)
  bioc_pkgs <- c("GenomicRanges", "IRanges", "BSgenome", "Biostrings")
  for (pkg in bioc_pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE)
      library(pkg, character.only = TRUE, quietly = TRUE)
    }
  }
  
  # Verify dndscv installation
  if (!require("dndscv", character.only = TRUE, quietly = TRUE)) {
    stop("Error: dndscv package not found. Please install it first.")
  } else {
    library(dndscv)
    cat(paste("✅ dndscv Version: ", as.character(packageVersion("dndscv")), "\n", sep = ""))
    
    # Check if newer version is needed for buildref
    if (packageVersion("dndscv") < "1.1.0") {
      cat("⚠️  Updating dndscv to support buildref function...\n")
      devtools::install_github("im3sanger/dndscv")
      library(dndscv)
    }
  }
}

# Execute dependency installation
install_dependencies()


# ==============================
# Step 2: Configuration Setup
# ==============================
# Analysis configuration
setup_analysis_config <- function(species = "human", genome_build = "hg38") {
  config <- list(
    species = species,
    genome_build = genome_build,
    use_prebuilt_ref = ifelse(species == "human", TRUE, FALSE),
    refdb = ifelse(species == "human", genome_build, NULL),
    numcode = ifelse(species == "human", 1, 1),  # Standard genetic code
    excludechrs = ifelse(species == "human", NULL, "MT")
  )
  
  cat(paste("🔧 Analysis Configuration:\n"))
  cat(paste("  Species: ", config$species, "\n"))
  cat(paste("  Genome Build: ", config$genome_build, "\n"))
  cat(paste("  Use Prebuilt Reference: ", ifelse(config$use_prebuilt_ref, "Yes", "No"), "\n"))
  
  return(config)
}


# ==============================
# Step 3: File Path Configuration
# ==============================
setup_file_paths <- function(base_dir, species, genome_build) {
  paths <- list(
    base_dir = base_dir,
    maf_file = file.path(base_dir, paste0(species, "_tumor.maf")),
    results_dir = file.path(base_dir, "dndscv_results"),
    # For custom species only
    cds_table = ifelse(species == "human", NULL, file.path(base_dir, paste0(species, "_", genome_build, "_cds.tsv"))),
    genome_fasta = ifelse(species == "human", NULL, file.path(base_dir, paste0(genome_build, "_genome.fa"))),
    refcds_output = ifelse(species == "human", NULL, file.path(base_dir, paste0(genome_build, "_refcds.rda")))
  )
  
  # Create directories
  if (!dir.exists(paths$results_dir)) {
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  cat(paste("📁 File Paths Setup:\n"))
  cat(paste("  Base Directory: ", paths$base_dir, "\n"))
  cat(paste("  Results Directory: ", paths$results_dir, "\n"))
  
  return(paths)
}


# ==============================
# Step 4: Preprocess CDS Table (For Non-human Species)
# ==============================
preprocess_cds_table <- function(cds_path) {
  if (is.null(cds_path) || !file.exists(cds_path)) {
    stop("CDS table file not found or not provided.")
  }
  
  # Read and validate CDS table
  cat(paste("📥 Reading CDS table: ", cds_path, "\n"))
  cds <- read.table(cds_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  required_cols <- c("gene.id", "gene.name", "cds.id", "chr", "chr.coding.start", 
                     "chr.coding.end", "cds.start", "cds.end", "length", "strand")
  if (!all(required_cols %in% colnames(cds))) {
    stop("CDS table missing required columns!")
  }
  
  # Clean CDS data
  cds_clean <- cds %>%
    dplyr::filter(chr != "MT") %>%  # Exclude mitochondria
    dplyr::mutate(chr = gsub("^chr", "", chr)) %>%  # Unify chromosome names
    dplyr::filter(!grepl("N|n", cds.id))  # Remove ambiguous transcripts
  
  # Save clean CDS table
  clean_cds_path <- gsub("\\.tsv", "_clean.tsv", cds_path)
  write.table(cds_clean, clean_cds_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("✅ Clean CDS saved to: ", clean_cds_path, "\n"))
  
  return(clean_cds_path)
}


# ==============================
# Step 5: Build Species-specific RefCDS
# ==============================
build_species_refcds <- function(cds_path, genome_path, out_path, excludechrs = "MT", numcode = 1) {
  if (is.null(cds_path) || is.null(genome_path)) {
    stop("CDS or genome file path is missing for RefCDS building.")
  }
  
  cat(paste("🔨 Building ", out_path, " RefCDS (this may take 30-60 mins)...\n"))
  
  buildref(
    cdsfile = cds_path,
    genomefile = genome_path,
    outfile = out_path,
    excludechrs = excludechrs,
    numcode = numcode,
    verbose = TRUE
  )
  
  if (file.exists(out_path)) {
    load(out_path)
    cat(paste("✅ RefCDS built! Genes included: ", length(RefCDS), "\n"))
    return(out_path)
  } else {
    stop("❌ RefCDS build failed! Check Fasta/CDS compatibility.")
  }
}


# ==============================
# Step 6: MAF Data Preprocessing
# ==============================
prepare_maf_data <- function(maf_path, refcds_path = NULL, use_prebuilt_ref = TRUE) {
  # Check if MAF file exists
  if (!file.exists(maf_path)) {
    stop(paste("Error: MAF file not found at path: ", maf_path, sep = ""))
  }
  
  # Read MAF
  cat(paste("📥 Reading MAF file: ", maf_path, "\n"))
  maf_data <- read.delim(
    maf_path,
    comment.char = "#",
    stringsAsFactors = FALSE,
    sep = "\t",
    quote = ""
  )
  
  # Check required MAF columns
  required_cols <- c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", 
                     "Reference_Allele", "Tumor_Seq_Allele2")
  missing_cols <- setdiff(required_cols, colnames(maf_data))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Missing required MAF columns: ", paste(missing_cols, collapse = ", "), sep = ""))
  }
  
  # If using custom RefCDS, load it to get valid chromosomes
  valid_chrs <- NULL
  if (!use_prebuilt_ref && !is.null(refcds_path) && file.exists(refcds_path)) {
    load(refcds_path)
    valid_chrs <- unique(sapply(RefCDS, function(g) g$chr))
  }
  
  # Convert to dndscv standard format
  mutations <- maf_data %>%
    dplyr::select(
      sampleID = Tumor_Sample_Barcode,
      chr = Chromosome,
      pos = Start_Position,
      ref = Reference_Allele,
      mut = Tumor_Seq_Allele2
    ) %>%
    dplyr::mutate(
      chr = gsub("^chr", "", chr),  # Remove "chr" prefix
      ref = toupper(ref),
      mut = toupper(mut)
    ) %>%
    # Filter for single-nucleotide variants
    dplyr::filter(
      nchar(ref) == 1,
      nchar(mut) == 1,
      ref %in% c("A", "T", "C", "G"),
      mut %in% c("A", "T", "C", "G")
    )
  
  # If using custom RefCDS, filter by valid chromosomes
  if (!is.null(valid_chrs)) {
    mutations <- mutations %>%
      dplyr::mutate(chr = ifelse(chr %in% valid_chrs, chr, NA)) %>%
      dplyr::filter(!is.na(chr))
  }
  
  # Remove duplicates
  mutations <- mutations %>%
    dplyr::distinct(sampleID, chr, pos, ref, mut, .keep_all = TRUE)
  
  # Print data summary
  cat(paste("📊 Total Mutations in MAF: ", nrow(maf_data), "\n"))
  cat(paste("📊 Valid Single-Nucleotide Variants: ", nrow(mutations), "\n"))
  cat(paste("📊 Number of Samples: ", length(unique(mutations$sampleID)), "\n"))
  
  return(mutations)
}


# ==============================
# Step 7: Run dNdScv Analysis
# ==============================
run_dndscv_analysis <- function(mutations, refdb = "hg38", refcds_path = NULL, 
                               use_prebuilt_ref = TRUE, numcode = 1) {
  cat(paste("🔬 Starting dNdScv Analysis\n"))
  
  # Determine reference database
  ref_database <- ifelse(use_prebuilt_ref, refdb, refcds_path)
  cat(paste("  Using reference database: ", ref_database, "\n"))
  
  # Run dndscv (compatible with both old/new version parameters)
  if ("mutations" %in% names(formals(dndscv))) {
    # New version: parameter name = "mutations"
    dndsout <- dndscv(
      mutations = mutations,
      refdb = ref_database,
      outmats = TRUE,
      cv = ifelse(use_prebuilt_ref, NULL, NULL),  # Disable human covariates for non-human species
      max_muts_per_gene_per_sample = 5,
      max_coding_muts_per_sample = 1000,
      numcode = numcode
    )
  } else {
    # Old version: parameter name = "x"
    dndsout <- dndscv(
      x = mutations,
      refdb = ref_database,
      outmats = TRUE,
      cv = ifelse(use_prebuilt_ref, NULL, NULL),
      max_muts_per_gene_per_sample = 5,
      max_coding_muts_per_sample = 1000,
      numcode = numcode
    )
  }
  
  # Verify analysis results
  if (is.null(dndsout$sel_cv)) {
    stop("Error: dNdScv analysis failed - 'sel_cv' result not generated.")
  } else {
    cat(paste("✅ Analysis Completed! Total Genes Analyzed: ", nrow(dndsout$sel_cv), "\n"))
  }
  
  return(dndsout)
}


# ==============================
# Step 8: Filter Candidate Genes
# ==============================
filter_candidate_genes <- function(dndsout, results_dir, p_threshold = 0.01) {
  # Handle output directory
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    cat(paste("📁 Created output directory: ", results_dir, "\n"))
  }
  
  # Extract sel_cv results
  sel_cv <- dndsout$sel_cv
  if (!"pallsubs_cv" %in% colnames(sel_cv)) {
    stop("Error: 'pallsubs_cv' column not found in sel_cv. Check data structure.")
  }
  
  # Filter candidate genes by p-value
  candidate_genes <- sel_cv %>%
    dplyr::filter(pallsubs_cv < p_threshold) %>%
    dplyr::arrange(pallsubs_cv) %>%
    dplyr::select(
      gene_name,
      pallsubs_cv,
      wmis_cv, wnon_cv,
      n_syn, n_mis, n_non,
      pmis_cv, ptrunc_cv
    ) %>%
    dplyr::mutate(
      Selection_Strength = dplyr::case_when(
        pallsubs_cv < 1e-04 ~ "Extreme Significance (p<0.0001)",
        pallsubs_cv < 1e-03 ~ "High Significance (p<0.001)",
        pallsubs_cv < p_threshold ~ "Significance (p<0.01)"
      )
    )
  
  # Save candidate gene list
  save_table <- function(data, filename) {
    file_path <- file.path(results_dir, filename)
    tryCatch({
      write.table(data, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
      cat(paste("✅ Saved: ", file_path, "\n"))
    }, error = function(e) {
      alt_path <- file.path(getwd(), filename)
      write.table(data, alt_path, sep = "\t", row.names = FALSE, quote = FALSE)
      cat(paste("⚠️  Saved to fallback path: ", alt_path, "\n"))
    })
  }
  
  save_table(candidate_genes, "candidate_genes_pvalue_based.txt")
  
  # Print result summary
  cat(paste("\n📈 Result Summary (p <", p_threshold, "):\n"))
  cat(paste("Total Genes Analyzed: ", nrow(sel_cv), "\n"))
  cat(paste("Number of Candidate Genes: ", nrow(candidate_genes), "\n"))
  
  if (nrow(candidate_genes) > 0) {
    cat("\n🏆 Top 10 Most Significant Candidate Genes:\n")
    print(head(candidate_genes[, c("gene_name", "pallsubs_cv", "wmis_cv", "Selection_Strength")], 10), digits = 4)
  } else {
    cat(paste("\n⚠️  No candidate genes found with p <", p_threshold, ". Try relaxing the threshold (e.g., p<0.05).\n"))
  }
  
  # Save global dN/dS results
  if (!is.null(dndsout$globaldnds)) {
    save_table(dndsout$globaldnds, "global_dnds_estimation.txt")
    cat("\n🌐 Global dN/dS Estimation:\n")
    print(dndsout$globaldnds, digits = 3)
  }
  
  # Save full analysis results
  tryCatch({
    saveRDS(dndsout, file.path(results_dir, "dndscv_full_results.rds"))
    cat(paste("✅ Full results saved to: ", file.path(results_dir, "dndscv_full_results.rds"), "\n"))
  }, error = function(e) {
    alt_rds <- file.path(getwd(), "dndscv_full_results.rds")
    saveRDS(dndsout, alt_rds)
    cat(paste("⚠️  Full results saved to fallback path: ", alt_rds, "\n"))
  }
  )
  
  # Return result summary
  return(list(
    candidate_genes = candidate_genes,
    global_dnds = if (!is.null(dndsout$globaldnds)) dndsout$globaldnds else NULL,
    output_dir = results_dir,
    p_threshold_used = p_threshold
  ))
}


# ==============================
# Step 9: Visualization Functions
# ==============================
# Safe PDF saving function
safe_pdf <- function(plot_code, filename, output_dir, width = 14, height = 8) {
  file_path <- file.path(output_dir, filename)
  tryCatch({
    pdf(file_path, width = width, height = height)
    eval(plot_code)
    dev.off()
    cat(paste("✅ Plot saved: ", file_path, "\n"))
    return(TRUE)
  }, error = function(e) {
    cat(paste("❌ Plot save failed: ", filename, " | Error: ", e$message, "\n"))
    return(FALSE)
  })
}

# Optimized visualization function
visualize_results <- function(dndsout, results_summary) {
  output_dir <- results_summary$output_dir
  if (is.null(output_dir) || !dir.exists(output_dir)) {
    output_dir <- getwd()
    cat(paste("⚠️  Using fallback output directory: ", output_dir, "\n"))
  }
  
  cat(paste("\n🎨 Starting Visualization\n"))
  
  # 1. Candidate Genes dN/dS Bar Plot
  if (nrow(results_summary$candidate_genes) > 0) {
    top_genes <- head(results_summary$candidate_genes, 20)
    p_threshold <- results_summary$p_threshold_used
    max_dnds <- max(top_genes$wmis_cv)
    
    safe_pdf(
      plot_code = quote({
        # Increase bottom margin for X-axis label
        par(mar = c(18, 10, 6, 2) + 0.1)
        
        bar_heights <- top_genes$wmis_cv
        gene_names <- top_genes$gene_name
        
        bar_plot <- barplot(
          height = bar_heights,
          names.arg = gene_names,
          las = 2,
          col = adjustcolor("steelblue", alpha.f = 0.8),
          border = "white",
          main = paste("Top 20 Candidate Genes (p <", p_threshold, ") - Missense dN/dS Ratio"),
          main.cex = 1.2,
          xlab = "",
          ylab = "Missense dN/dS Ratio (wmis_cv)",
          ylab.cex = 1.1,
          ylim = c(0, max_dnds * 1.2),
          yaxt = "n"
        )
        
        # Add X-axis label lower
        mtext(
          text = "Candidate Genes",
          side = 1,
          line = 12,
          cex = 1.1,
          font = 2
        )
        
        # Y-axis ticks
        y_ticks <- seq(0, ceiling(max_dnds * 1.2), by = ceiling(max_dnds * 1.2 / 6))
        axis(2, at = y_ticks, labels = as.character(y_ticks), cex.axis = 1.0, las = 1)
        
        # Significance markers
        marker_height <- max_dnds * 0.03
        for (i in 1:length(bar_heights)) {
          p_val <- top_genes$pallsubs_cv[i]
          y_pos <- bar_heights[i] + marker_height
          text(x = bar_plot[i], y = y_pos, 
               labels = ifelse(p_val < 1e-04, "***", ifelse(p_val < 1e-03, "**", "*")),
               cex = 1.8, col = "darkred")
        }
        
        # Legend and grid
        legend("topright",
               legend = c("***: p < 0.0001", "**: p < 0.001", paste("*: p <", p_threshold)),
               bty = "n", cex = 1.0, col = "darkred")
        grid(nx = 0, ny = NULL, lty = 2, col = adjustcolor("gray", alpha.f = 0.3)
        )
      }),
      filename = "candidate_genes_dnds_barplot.pdf",
      output_dir = output_dir,
      width = 18,
      height = 10
    )
  } else {
    cat("⚠️  No candidate genes, skipping dN/dS bar plot\n")
  }
  
  # 2. Global dN/dS Estimation Plot
  if (!is.null(results_summary$global_dnds) && nrow(results_summary$global_dnds) >= 5) {
    safe_pdf(
      plot_code = quote({
        par(mar = c(14, 8, 6, 2) + 0.1)
        mutation_types <- c("Missense", "Nonsense", "Splice Site", "Truncating", "All Nonsynonymous")
        bar_cols <- adjustcolor(c("steelblue", "firebrick", "forestgreen", "orange", "purple"), alpha.f = 0.8)
        
        bar_plot <- barplot(
          height = results_summary$global_dnds$mle,
          names.arg = mutation_types,
          las = 2,
          col = bar_cols,
          border = "white",
          main = "Global dN/dS Estimation (Maximum Likelihood)",
          main.cex = 1.2,
          xlab = "",
          ylab = "dN/dS Ratio (MLE)",
          ylab.cex = 1.1,
          ylim = c(0, max(results_summary$global_dnds$cihigh) * 1.2),
          yaxt = "n"
        )
        
        # Add X-axis label
        mtext(
          text = "Mutation Types",
          side = 1,
          line = 10,
          cex = 1.1,
          font = 2
        )
        
        # Other elements
        y_ticks <- seq(0, ceiling(max(results_summary$global_dnds$cihigh) * 1.2), 
                      by = ceiling(max(results_summary$global_dnds$cihigh) * 1.2 / 6))
        axis(2, at = y_ticks, labels = as.character(y_ticks), cex.axis = 1.0, las = 1)
        arrows(x0 = bar_plot, y0 = results_summary$global_dnds$cilow, y1 = results_summary$global_dnds$cihigh,
               angle = 90, code = 3, length = 0.05, lwd = 1.5, col = "darkgray")
        abline(h = 1, lty = 2, col = "darkred", lwd = 1.5)
        legend("topright", legend = "Red dashed line: Neutral Evolution (dN/dS = 1)",
               bty = "n", cex = 1.0, col = "darkred", lty = 2)
        grid(nx = 0, ny = NULL, lty = 2, col = adjustcolor("gray", alpha.f = 0.3))
      }),
      filename = "global_dnds_estimation_plot.pdf",
      output_dir = output_dir,
      width = 12,
      height = 8
    )
  }
  
  # 3. Mutation Type Distribution Pie Chart
  if (!is.null(dndsout$annotmuts) && nrow(dndsout$annotmuts) > 0 && "impact" %in% colnames(dndsout$annotmuts)) {
    safe_pdf(
      plot_code = quote({
        mut_count <- table(dndsout$annotmuts$impact)
        small_types <- names(mut_count)[mut_count / sum(mut_count) < 0.05]
        if (length(small_types) > 0) {
          mut_count["Other (Rare)"] <- sum(mut_count[small_types])
          mut_count <- mut_count[!names(mut_count) %in% small_types]
        }
        
        pie_colors <- rainbow(length(mut_count))
        pie(
          mut_count,
          col = pie_colors,
          main = "Mutation Type Distribution (Based on 'impact' Column)",
          main.cex = 1.2,
          cex = 0.9
        )
        
        legend("bottomright",
               legend = paste(names(mut_count), " (", round(mut_count/sum(mut_count)*100, 1), "%)", sep = ""),
               fill = pie_colors,
               bty = "n",
               cex = 0.8,
               ncol = 2)
      }),
      filename = "mutation_type_distribution.pdf",
      output_dir = output_dir,
      width = 9,
      height = 9
    )
  } else {
    cat("⚠️  Skipping mutation pie chart (missing 'impact' column or data)\n")
  }
  
  cat(paste("\n🎉 Visualization Completed!\n"))
}


# ==============================
# Step 10: Main Analysis Pipeline
# ==============================
run_universal_dndscv_analysis <- function(
    base_dir = getwd(),
    species = "human",
    genome_build = "hg38",
    maf_file = NULL,
    cds_table = NULL,
    genome_fasta = NULL,
    p_threshold = 0.01
) {
  cat(paste("🚀 Starting Universal dNdScv Analysis for ", species, " (", genome_build, ")\n\n"))
  
  # Setup configuration
  config <- setup_analysis_config(species, genome_build)
  
  # Setup file paths
  file_paths <- setup_file_paths(base_dir, species, genome_build)
  
  # Override with provided file paths if given
  if (!is.null(maf_file)) file_paths$maf_file <- maf_file
  if (!is.null(cds_table)) file_paths$cds_table <- cds_table
  if (!is.null(genome_fasta)) file_paths$genome_fasta <- genome_fasta
  
  # For non-human species, preprocess CDS and build RefCDS
  refcds_path <- NULL
  if (!config$use_prebuilt_ref) {
    cat(paste("\n🧬 Building Reference Database for ", species, "\n"))
    clean_cds_path <- preprocess_cds_table(file_paths$cds_table)
    refcds_path <- build_species_refcds(
      clean_cds_path, 
      file_paths$genome_fasta, 
      file_paths$refcds_output,
      config$excludechrs,
      config$numcode
    )
  }
  
  # Preprocess MAF data
  cat(paste("\n📋 Preprocessing MAF Data\n"))
  mutations <- prepare_maf_data(
    file_paths$maf_file, 
    refcds_path, 
    config$use_prebuilt_ref
  )
  
  # Run dNdScv analysis
  cat(paste("\n🔬 Running dNdScv Analysis\n"))
  dndsout <- run_dndscv_analysis(
    mutations,
    refdb = config$refdb,
    refcds_path = refcds_path,
    use_prebuilt_ref = config$use_prebuilt_ref,
    numcode = config$numcode
  )
  
  # Filter candidate genes
  cat(paste("\n🔍 Filtering Candidate Genes\n"))
  results_summary <- filter_candidate_genes(
    dndsout, 
    file_paths$results_dir, 
    p_threshold
  )
  
  # Generate visualizations
  cat(paste("\n📊 Generating Visualizations\n"))
  visualize_results(dndsout, results_summary)
  
  cat(paste("\n🎊 Universal dNdScv Analysis Completed!\n"))
  cat(paste("📁 Results saved to: ", file_paths$results_dir, "\n"))
  
  return(list(
    mutations = mutations,
    dndsout = dndsout,
    results_summary = results_summary
  ))
}


# ==============================
# Usage Examples
# ==============================

# Example 1: Human analysis (hg38)
# results <- run_universal_dndscv_analysis(
#   base_dir = "/path/to/your/data",
#   species = "human",
#   genome_build = "hg38",
#   maf_file = "/path/to/your/human.maf",
#   p_threshold = 0.01
# )

# Example 2: Rat analysis (rn7)
# results <- run_universal_dndscv_analysis(
#   base_dir = "/path/to/your/data",
#   species = "rat",
#   genome_build = "rn7",
#   maf_file = "/path/to/your/rat.maf",
#   cds_table = "/path/to/your/rn7_cds.tsv",
#   genome_fasta = "/path/to/your/rn7_genome.fa",
#   p_threshold = 0.01
# )

# Example 3: Mouse analysis (mm39)
# results <- run_universal_dndscv_analysis(
#   base_dir = "/path/to/your/data",
#   species = "mouse",
#   genome_build = "mm39",
#   maf_file = "/path/to/your/mouse.maf",
#   cds_table = "/path/to/your/mm39_cds.tsv",
#   genome_fasta = "/path/to/your/mm39_genome.fa",
#   p_threshold = 0.01
# )
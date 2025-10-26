# Batch processing script for multiple samples + multi-sample joint analysis
# Function: Process all samples in a folder + generate joint CNV summary

# Install required packages
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("VariantAnnotation")) install.packages("VariantAnnotation")
if require("GenomicRanges")) install.packages("GenomicRanges")
if (!require("dplyr")) install.packages("dplyr")  # For joint statistics (using basic functions only)

# Load packages
library(ggplot2)
library(VariantAnnotation)
library(GenomicRanges)

# --------------------------
# 1. Configuration parameters
# --------------------------
input_dir <- "./delly_cnv_results"        # Folder containing all samples 
output_root <- "./result_plot"  # Root directory for all outputs
cn_normal <- 2                 # Normal copy number (diploid)
cn_threshold <- 0.2            # Threshold for CNV calling (±0.2)
genome_version <- "hg38"        # Reference genome version

# Create output directories
if (!dir.exists(output_root)) {
  dir.create(output_root, recursive = TRUE)
}
joint_plot_dir <- file.path(output_root, "joint_analysis")  # For multi-sample plots
if (!dir.exists(joint_plot_dir)) {
  dir.create(joint_plot_dir)
}

# --------------------------
# 2. Identify all samples in input folder
# --------------------------
# Get all VCF files in input directory
vcf_files <- list.files(input_dir, pattern = "\\.vcf$|\\.vcf\\.gz$", full.names = TRUE)
if (length(vcf_files) == 0) {
  stop("No VCF files found in input directory: ", input_dir)
}

# Match paired .cov files (assume cov filename matches VCF prefix, e.g., S1.vcf → S1.cov)
sample_ids <- sub(pattern = "\\.vcf$|\\.vcf\\.gz$", replacement = "", x = basename(vcf_files))
cov_files <- file.path(input_dir, paste0(sample_ids, ".cov.gz"))

# Check for missing .cov files
missing_cov <- which(!file.exists(cov_files))
if (length(missing_cov) > 0) {
  stop("Missing .cov files for samples: ", paste(sample_ids[missing_cov], collapse = ", "))
}

cat("Found", length(sample_ids), "samples to process:\n")
cat(paste(sample_ids, collapse = ", "), "\n\n")

# --------------------------
# 3. Process each sample (loop)
# --------------------------
# Initialize global statistics table for joint analysis
global_stats <- data.frame(
  sample = character(),
  chrom = character(),
  total_cnv = integer(),
  dup_count = integer(),
  del_count = integer(),
  stringsAsFactors = FALSE
)

# Process each sample
for (i in seq_along(sample_ids)) {
  sample_name <- sample_ids[i]
  cat("=====================================\n")
  cat("Processing sample:", sample_name, "\n")
  
  # Create sample-specific output directory
  sample_output <- file.path(output_root, sample_name)
  if (!dir.exists(sample_output)) {
    dir.create(sample_output)
  }
  sample_plot_dir <- file.path(sample_output, "chrom_plots")
  if (!dir.exists(sample_plot_dir)) {
    dir.create(sample_plot_dir)
  }
  
  # --------------------------
  # 3.1 Process coverage file
  # --------------------------
  cov_file <- cov_files[i]
  cov_data <- read.delim(cov_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cn_col <- grep("CN", colnames(cov_data), ignore.case = TRUE, value = TRUE)
  
  # Remove NA values
  cat("  Processing coverage data...\n")
  cov_na <- sum(is.na(cov_data$chr) | is.na(cov_data$start) | is.na(cov_data$end) | is.na(cov_data[[sample_name]]))
  if (cov_na > 0) {
    cat(paste0("    Removed ", cov_na, " rows with NA values\n"))
    cov_data <- cov_data[!is.na(cov_data$chr) & !is.na(cov_data$start) & !is.na(cov_data$end) & !is.na(cov_data[[sample_name]]), ]
  }
  
  # Prepare coverage data
  cov_clean <- data.frame(
    chrom = cov_data$chr,
    start = cov_data$start,
    end = cov_data$end,
    cn = cov_data[[cn_col]],  # Use sample-specific CN column
    stringsAsFactors = FALSE
  )
  cov_clean$mid_pos <- (cov_clean$start + cov_clean$end) / 1e6
  cov_clean$cn_status <- ifelse(
    cov_clean$cn > cn_normal + cn_threshold, "DUP",
    ifelse(
      cov_clean$cn < cn_normal - cn_threshold, "DEL", "NORMAL"
    )
  )
  
  all_chroms <- unique(cov_clean$chrom)
  cat("    Detected chromosomes:", paste(all_chroms, collapse = ", "), "\n")
  
  # --------------------------
  # 3.2 Process VCF file
  # --------------------------
  vcf_file <- vcf_files[i]
  vcf <- readVcf(vcf_file, genome = genome_version)
  vcf_info <- data.frame(
    chrom = as.character(seqnames(vcf)),
    start = start(vcf),
    end = end(vcf),
    id = rownames(vcf),
    stringsAsFactors = FALSE
  )
  
  # Remove NA values
  cat("  Processing VCF data...\n")
  vcf_na <- sum(is.na(vcf_info$chrom) | is.na(vcf_info$start) | is.na(vcf_info$end))
  if (vcf_na > 0) {
    cat(paste0("    Removed ", vcf_na, " rows with NA values\n"))
    vcf_info <- vcf_info[!is.na(vcf_info$chrom) & !is.na(vcf_info$start) & !is.na(vcf_info$end), ]
  }
  
  # --------------------------
  # 3.3 Generate per-chromosome plots
  # --------------------------
  sample_stats <- data.frame(
    chrom = character(),
    total_cnv = integer(),
    dup_count = integer(),
    del_count = integer(),
    stringsAsFactors = FALSE
  )
  
  for (chrom in all_chroms) {
    # Filter data for current chromosome
    cov_plot <- cov_clean[cov_clean$chrom == chrom, ]
    if (nrow(cov_plot) == 0) {
      cat(paste0("    Warning: No coverage data for chromosome ", chrom, "\n"))
      next
    }
    
    vcf_chrom <- vcf_info[vcf_info$chrom == chrom, , drop = FALSE]
    vcf_plot <- data.frame()
    
    if (nrow(vcf_chrom) > 0) {
      # Match VCF regions with coverage data
      cov_gr <- GRanges(
        seqnames = cov_plot$chrom,
        ranges = IRanges(start = cov_plot$start, end = cov_plot$end),
        cn = cov_plot$cn
      )
      vcf_gr <- GRanges(
        seqnames = vcf_chrom$chrom,
        ranges = IRanges(start = vcf_chrom$start, end = vcf_chrom$end)
      )
      
      overlaps <- findOverlaps(vcf_gr, cov_gr)
      vcf_chrom$mean_cn <- NA
      if (length(overlaps) > 0) {
        overlap_df <- data.frame(
          queryHits = queryHits(overlaps),
          cn = cov_gr$cn[subjectHits(overlaps)]
        )
        vcf_cn <- aggregate(cn ~ queryHits, data = overlap_df, FUN = mean)
        vcf_chrom$mean_cn[vcf_cn$queryHits] <- vcf_cn$cn
      }
      
      # Filter invalid CNV calls
      vcf_chrom <- vcf_chrom[!is.na(vcf_chrom$mean_cn), , drop = FALSE]
      if (nrow(vcf_chrom) > 0) {
        vcf_chrom$cnv_type <- ifelse(
          vcf_chrom$mean_cn > cn_normal + cn_threshold, "DUP",
          ifelse(
            vcf_chrom$mean_cn < cn_normal - cn_threshold, "DEL", "UNKNOWN"
          )
        )
        vcf_chrom$start_mb <- vcf_chrom$start / 1e6
        vcf_chrom$end_mb <- vcf_chrom$end / 1e6
        vcf_plot <- vcf_chrom
      }
    }
    
    # Generate plot
    p <- ggplot() +
      geom_line(
        data = cov_plot,
        aes(x = mid_pos, y = cn, color = cn_status),
        linewidth = ifelse(cov_plot$cn_status == "NORMAL", 0.3, 0.5),
        alpha = 0.8
      ) +
      {if (nrow(vcf_plot) > 0) geom_rect(
        data = vcf_plot[vcf_plot$cnv_type != "UNKNOWN", , drop = FALSE],
        aes(
          xmin = start_mb, xmax = end_mb,
          ymin = -Inf, ymax = Inf,
          fill = cnv_type
        ),
        alpha = 0.4,
        color = aes(color = cnv_type),
        size = 0.3
      )} +
      geom_hline(yintercept = cn_normal, color = "black", linetype = "dashed", linewidth = 0.5) +
      scale_color_manual(
        values = c("DUP" = "#E41A1C", "DEL" = "#377EB8", "NORMAL" = "#999999"),
        name = "Copy Number Status"
      ) +
      scale_fill_manual(
        values = c("DUP" = "#E41A1C", "DEL" = "#377EB8"),
        name = "VCF CNV Type"
      ) +
      scale_color_manual(
        values = c("DUP" = "#E41A1C", "DEL" = "#377EB8"),
        guide = "none"
      ) +
      labs(
        title = paste0(sample_name, " - Chromosome ", chrom),
        x = "Position (Mb)",
        y = "Predicted Copy Number"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(
      filename = file.path(sample_plot_dir, paste0(sample_name, "_", chrom, "_cnv.pdf")),
      plot = p,
      width = 12, height = 6
    )
    
    # Update sample statistics
    total <- if (nrow(vcf_plot) > 0) nrow(vcf_plot) else 0
    dup <- if (nrow(vcf_plot) > 0) sum(vcf_plot$cnv_type == "DUP") else 0
    del <- if (nrow(vcf_plot) > 0) sum(vcf_plot$cnv_type == "DEL") else 0
    sample_stats <- rbind(sample_stats, data.frame(
      chrom = chrom,
      total_cnv = total,
      dup_count = dup,
      del_count = del,
      stringsAsFactors = FALSE
    ))
  }
  
  # Generate genome-wide summary plot for individual sample
  sample_long <- data.frame(
    chrom = rep(sample_stats$chrom, 2),
    cnv_type = c(rep("DUP", nrow(sample_stats)), rep("DEL", nrow(sample_stats))),
    count = c(sample_stats$dup_count, sample_stats$del_count),
    stringsAsFactors = FALSE
  )
  
  # Plot genome-wide CNV distribution for individual sample
  sample_genome_plot <- ggplot(sample_long, aes(x = chrom, y = count)) +
    geom_bar(aes(fill = cnv_type), stat = "identity", position = "dodge", width = 0.7) +
    geom_text(
      data = sample_stats,
      aes(x = chrom, y = total_cnv + 1, label = total_cnv),
      color = "black", size = 3.5, show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c("DUP" = "#E41A1C", "DEL" = "#377EB8"),
      name = "Legend",
      labels = c("DUP: Red (copy number gain)", "DEL: Blue (copy number loss)")
    ) +
    labs(
      title = paste0(sample_name, " - Genome-wide CNV Summary"),
      subtitle = "Black numbers: Total CNV per chromosome",
      x = "Chromosome",
      y = "Number of CNV Regions"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # Save genome-wide summary plot for individual sample
  ggsave(
    filename = file.path(sample_output, paste0(sample_name, "_genome_wide_summary.pdf")),
    plot = sample_genome_plot,
    width = 16, height = 8, device = "pdf"
  )
  
  # Save sample-specific statistics
  write.csv(
    sample_stats,
    file.path(sample_output, paste0(sample_name, "_cnv_stats.csv")),
    row.names = FALSE
  )
  
  # Update global statistics for joint analysis
  global_stats <- rbind(global_stats, data.frame(
    sample = sample_name,
    chrom = sample_stats$chrom,
    total_cnv = sample_stats$total_cnv,
    dup_count = sample_stats$dup_count,
    del_count = sample_stats$del_count,
    stringsAsFactors = FALSE
  ))
  
  cat("  Sample", sample_name, "processing complete\n\n")
}

# --------------------------
# 4. Generate multi-sample joint analysis results
# --------------------------
cat("Generating multi-sample joint analysis...\n")

# --------------------------
# 4.1 Total CNV count per sample (bar plot)
# --------------------------
total_per_sample <- aggregate(
  total_cnv ~ sample,
  data = global_stats,
  FUN = sum
)

p_total <- ggplot(total_per_sample, aes(x = sample, y = total_cnv, fill = sample)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    aes(label = total_cnv),
    vjust = -0.3, size = 3.5
  ) +
  labs(
    title = "Total CNV Count per Sample",
    x = "Sample",
    y = "Total Number of CNVs"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  filename = file.path(joint_plot_dir, "total_cnv_per_sample.pdf"),
  plot = p_total,
  width = 10, height = 6
)

# --------------------------
# 4.2 DUP/DEL distribution across samples (stacked bar plot)
# --------------------------
dup_del_per_sample <- aggregate(
  cbind(dup_count, del_count) ~ sample,  
  data = global_stats,
  FUN = sum
)

# Reshape for plotting
dup_del_long <- data.frame(
  sample = rep(dup_del_per_sample$sample, 2),
  cnv_type = c(rep("DUP", nrow(dup_del_per_sample)), rep("DEL", nrow(dup_del_per_sample))),
  count = c(dup_del_per_sample$dup_count, dup_del_per_sample$del_count),
  stringsAsFactors = FALSE
)

p_dup_del <- ggplot(dup_del_long, aes(x = sample, y = count, fill = cnv_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c("DUP" = "#E41A1C", "DEL" = "#377EB8"),
    name = "CNV Type"
  ) +
  labs(
    title = "DUP/DEL Distribution Across Samples",
    x = "Sample",
    y = "Number of CNVs"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(joint_plot_dir, "dup_del_distribution.pdf"),
  plot = p_dup_del,
  width = 10, height = 6
)

# --------------------------
# 4.3 Chromosome-specific CNV burden (heatmap-like plot)
# --------------------------
# Reshape data for heatmap
chrom_burden <- global_stats[, c("sample", "chrom", "total_cnv")]

p_chrom <- ggplot(chrom_burden, aes(x = chrom, y = sample, fill = total_cnv)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(name = "Total CNVs") +
  labs(
    title = "Chromosome-Specific CNV Burden Across Samples",
    x = "Chromosome",
    y = "Sample"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(joint_plot_dir, "chromosome_cnv_burden.pdf"),
  plot = p_chrom,
  width = 12, height = 8
)

# --------------------------
# 5. Save global statistics
# --------------------------
write.csv(
  global_stats,
  file.path(output_root, "multi_sample_cnv_stats.csv"),
  row.names = FALSE
)

cat("\nAll samples processed successfully!\n")
cat("Results saved to:\n")
cat("  - Per-sample plots and stats: ", output_root, "\n")
cat("  - Multi-sample joint analysis: ", joint_plot_dir, "\n")
cat("  - Global statistics table: ", file.path(output_root, "multi_sample_cnv_stats.csv"), "\n")
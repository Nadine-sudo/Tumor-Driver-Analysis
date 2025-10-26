#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import re
from collections import Counter, defaultdict
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import glob
import gzip
from typing import Dict, List, Tuple, Optional

# --------------------------
# 1. Configuration & Utility Functions
# --------------------------
# Define species mapping - customize this according to your samples
SPECIES_MAPPING = {
    # Add your sample_id patterns and corresponding species here
    # Example:
    # "Rat": "rat",
    # "Human": "human",
    # "Mouse": "mouse"
    "Rat": "rat",
    "Cell": "cell"
}

def set_plot_style():
    """Set global plot style for consistency (English labels)"""
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'font.size': 11,
        'font.family': 'Arial',
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.dpi': 300,
        'figure.figsize': (10, 6)
    })

def is_vcf_header(line: str) -> bool:
    """Check if a line in VCF is a header line (##, #CHROM, or empty)"""
    line_stripped = line.strip()
    return (line_stripped.startswith("##") or 
            line_stripped.startswith("#CHROM") or 
            len(line_stripped) == 0)

def get_sample_vcfs(input_path: str) -> Dict[str, str]:
    """
    Get a dictionary of {sample_id: vcf_file_path} from input.
    - If input_path is a FOLDER: Auto-detect all .vcf/.vcf.gz files.
    - If input_path is a TEXT FILE: Read sample-VCF pairs (one per line: sample_id\tvcf_path).
    """
    sample_vcf_dict = {}
    
    if os.path.isdir(input_path):
        # Case 1: Input is a folder (auto-find VCFs)
        vcf_files = glob.glob(os.path.join(input_path, "*.vcf")) + glob.glob(os.path.join(input_path, "*.vcf.gz"))
        if not vcf_files:
            raise FileNotFoundError(f"No VCF files (.vcf/.vcf.gz) found in folder: {input_path}")
        
        # Extract sample ID from VCF filename (e.g., "Rat_Sample1_delly.vcf" → "Rat_Sample1")
        for vcf in vcf_files:
            filename = os.path.basename(vcf).replace(".vcf", "").replace(".vcf.gz", "")
            # Simplify sample ID (remove suffixes like "_delly" or "_sv")
            sample_id = filename.split("_delly")[0].split("_sv")[0].split("_SV")[0]
            sample_vcf_dict[sample_id] = vcf
        print(f"Auto-detected {len(sample_vcf_dict)} samples from folder: {input_path}")
    
    elif os.path.isfile(input_path) and input_path.endswith(".txt"):
        # Case 2: Input is a text file (sample-VCF pairs)
        with open(input_path, "r") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue  # Skip comments/empty lines
                parts = line.split("\t")
                if len(parts) != 2:
                    raise ValueError(f"Invalid format in line {line_num} of {input_path}: Use 'sample_id\tvcf_path'")
                sample_id, vcf_path = parts
                if not os.path.exists(vcf_path):
                    raise FileNotFoundError(f"VCF for {sample_id} not found: {vcf_path}")
                sample_vcf_dict[sample_id] = vcf_path
        print(f"Loaded {len(sample_vcf_dict)} samples from text file: {input_path}")
    
    else:
        raise ValueError(f"Input must be a folder or a .txt sample list: {input_path}")
    
    return sample_vcf_dict

# --------------------------
# 2. Batch VCF Parsing (Multiple Samples)
# --------------------------
def parse_delly_vcf(vcf_path: str, sample_id: str) -> pd.DataFrame:
    """
    Parse a single Delly VCF file and add sample ID to the output DataFrame.
    Returns: DataFrame with SV data + 'sample_id' + 'species' column.
    Species classification logic: uses user-defined SPECIES_MAPPING; fallback to "unknown"
    """
    # Read non-header lines (actual SV records)
    variant_lines = []
    open_func = open if not vcf_path.endswith(".gz") else gzip.open
    with open_func(vcf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if not is_vcf_header(line):
                variant_lines.append(line.strip())
    
    # Handle empty VCF (no SVs detected)
    if not variant_lines:
        print(f"⚠️ No SV records found for sample: {sample_id} (VCF: {os.path.basename(vcf_path)})")
        return pd.DataFrame(columns=["sample_id", "chr", "pos", "svtype", "length", "pe_support", "qual", "filter", "species"])
    
    # Parse each SV record
    sv_records = []
    for line in variant_lines:
        fields = line.split("\t")
        if len(fields) < 8:
            continue  # Skip malformed lines (fewer than 8 VCF columns)
        
        # Extract core VCF fields
        chrom = fields[0]
        pos = int(fields[1]) if fields[1].isdigit() else 0
        qual = float(fields[5]) if fields[5].replace(".", "").isdigit() else 0.0
        filter_status = fields[6] if fields[6] != "." else "PASS"
        info_str = fields[7]
        
        # Parse INFO field (extract Delly-specific SV attributes)
        info = defaultdict(str)
        for attr in info_str.split(";"):
            if "=" in attr:
                key, value = attr.split("=", 1)
                info[key.upper()] = value  # Standardize keys to uppercase (e.g., svtype → SVTYPE)
        
        # Extract SV-specific metrics
        sv_type = info.get("SVTYPE", "UNKNOWN")
        sv_len = info.get("SVLEN", "0")
        # Clean SV length (handle lists/negative values: e.g., "-1200,1300" → 1200)
        sv_len = abs(int(sv_len.split(",")[0])) if sv_len.replace("-", "").replace(",", "").isdigit() else 0
        pe_support = int(info.get("PE", "0")) if info.get("PE", "0").isdigit() else 0
        
        # --------------------------
        # Core: Species classification logic (user-defined)
        # --------------------------
        species = "unknown"  # Default species
        for pattern, species_name in SPECIES_MAPPING.items():
            if pattern in sample_id:
                species = species_name
                break
        
        # Append record with sample ID and species
        sv_records.append({
            "sample_id": sample_id,
            "chr": chrom,
            "pos": pos,
            "svtype": sv_type,
            "length": sv_len,
            "pe_support": pe_support,
            "qual": qual,
            "filter": filter_status,
            "species": species
        })
    
    # Convert to DataFrame and validate
    sv_df = pd.DataFrame(sv_records)
    print(f"✅ Parsed {len(sv_df)} SVs for sample: {sample_id} (species: {species}, VCF: {os.path.basename(vcf_path)})")
    return sv_df

def batch_parse_vcfs(sample_vcf_dict: dict) -> pd.DataFrame:
    """Batch parse all VCFs in sample_vcf_dict and return a combined DataFrame"""
    combined_sv_df = []
    for sample_id, vcf_path in sample_vcf_dict.items():
        sample_df = parse_delly_vcf(vcf_path, sample_id)
        combined_sv_df.append(sample_df)
    
    # Combine all samples into one DataFrame
    if not combined_sv_df:
        raise ValueError("No SV records found across all samples (check Delly output validity)")
    
    combined_df = pd.concat(combined_sv_df, ignore_index=True)
    if len(combined_df) == 0:
        raise ValueError("No SV records found across all samples (check Delly output validity)")
    
    # Print species distribution statistics (for chromosome SV plot separation verification)
    species_count = combined_df["species"].value_counts()
    print(f"\n📊 Species distribution statistics (for chromosome SV plot separation):")
    for sp, count in species_count.items():
        print(f"  - {sp}: {count} SV records (involving {combined_df[combined_df['species']==sp]['sample_id'].nunique()} samples)")
    return combined_df

# --------------------------
# 3. Multiple-Sample Visualizations (All plots saved as PDF)
# --------------------------
def plot_total_svs_per_sample(combined_df: pd.DataFrame, output_dir: str):
    """Plot total number of SVs per sample (ranked bar chart)"""
    sv_counts = combined_df["sample_id"].value_counts().sort_values(ascending=True)
    
    fig, ax = plt.subplots(figsize=(10, max(6, len(sv_counts) * 0.3)))  # Dynamic height based on sample count
    bars = ax.barh(sv_counts.index, sv_counts.values, color="#2E86AB")
    
    # Add value labels to bars
    for bar in bars:
        width = bar.get_width()
        ax.text(width + max(sv_counts.values)*0.01, bar.get_y() + bar.get_height()/2, 
                str(int(width)), ha="left", va="center", fontweight="bold")
    
    ax.set_xlabel("Total Number of Structural Variants (SVs)")
    ax.set_ylabel("Sample ID")
    ax.set_title("Total SV Count Per Sample (Ranked)")
    ax.grid(axis="x", alpha=0.3)
    
    # Save as PDF format
    output_path = os.path.join(output_dir, "total_svs_per_sample.pdf")
    plt.tight_layout()
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close()
    print(f"📈 Saved total SV count plot (PDF): {output_path}")

def plot_sv_type_by_sample(combined_df: pd.DataFrame, output_dir: str):
    """Plot SV type distribution for each sample (grouped bar chart)"""
    # Pivot data: rows = samples, columns = SV types, values = count
    sv_type_pivot = combined_df.groupby(["sample_id", "svtype"]).size().unstack(fill_value=0)
    
    # Plot grouped bars
    fig, ax = plt.subplots(figsize=(max(12, len(sv_type_pivot.index)*0.5), 7))
    bar_width = 0.15
    x = np.arange(len(sv_type_pivot.index))
    colors = ["#A23B72", "#F18F01", "#C73E1D", "#2E86AB", "#8E4EC6", "#3D9970", "#FF6B6B"]  # More distinct colors
    
    # Plot each SV type as a separate bar group
    for i, (sv_type, color) in enumerate(zip(sv_type_pivot.columns, colors[:len(sv_type_pivot.columns)])):
        ax.bar(
            x + (i - len(sv_type_pivot.columns)/2 + 0.5) * bar_width,
            sv_type_pivot[sv_type],
            width=bar_width,
            label=sv_type,
            color=color
        )
    
    # Customize plot
    ax.set_xlabel("Sample ID")
    ax.set_ylabel("Number of SVs")
    ax.set_title("SV Type Distribution Across Samples")
    ax.set_xticks(x)
    ax.set_xticklabels(sv_type_pivot.index, rotation=45, ha="right")
    ax.legend(title="SV Type", bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.grid(axis="y", alpha=0.3)
    
    # Save as PDF format
    output_path = os.path.join(output_dir, "sv_type_by_sample.pdf")
    plt.tight_layout()
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close()
    print(f"📈 Saved SV type distribution plot (PDF): {output_path}")

def plot_chrom_sv_by_species(combined_df: pd.DataFrame, output_dir: str, top_chroms: int = 10):
    """
    Only this plot is separated by species (chrom_sv_by_sample_top10), saved as PDF
    Other plots do not distinguish species, maintaining original logic
    """
    # Get all unique species from the data
    all_species = combined_df["species"].unique().tolist()
    
    # Filter species with data
    valid_species = [sp for sp in all_species if len(combined_df[combined_df["species"]==sp]) > 0]
    
    if len(valid_species) == 0:
        print("⚠️ No valid species data found for chromosome SV plot")
        return
    
    for species in valid_species:
        species_df = combined_df[combined_df["species"] == species]
        # Get TOP N chromosomes with most SVs for this species
        top_chrom_list = species_df["chr"].value_counts().head(top_chroms).index.tolist()
        filtered_df = species_df[species_df["chr"].isin(top_chrom_list)]
        
        # Plot by sample (1 sample per row)
        samples = sorted(species_df["sample_id"].unique())
        n_samples = len(samples)
        
        if n_samples == 0:
            print(f"⚠️ No samples found for species {species}")
            continue
            
        fig, axes = plt.subplots(nrows=n_samples, ncols=1, figsize=(12, max(6, 3*n_samples)), sharex=True)
        axes = [axes] if n_samples == 1 else axes  # Handle single sample case
        
        for ax, sample in zip(axes, samples):
            # Count SVs per chromosome for current sample
            chrom_counts = filtered_df[filtered_df["sample_id"] == sample]["chr"].value_counts().reindex(top_chrom_list, fill_value=0)
            bars = ax.bar(chrom_counts.index, chrom_counts.values, color="#2E86AB", alpha=0.7)
            
            # Label only non-zero values (to avoid clutter)
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2, height + 0.1,
                            str(int(height)), ha="center", va="bottom", fontsize=9)
            
            ax.set_ylabel(f"SV Count\n(Sample: {sample})")
            ax.grid(axis="y", alpha=0.3)
            # Adjust y-axis range (to avoid label overflow)
            ax.set_ylim(0, max(chrom_counts.values) * 1.1 if max(chrom_counts.values) > 0 else 1)
        
        # Unified x-axis labels
        axes[-1].set_xlabel("Chromosome")
        axes[-1].set_xticks(range(len(top_chrom_list)))
        axes[-1].set_xticklabels(top_chrom_list)
        
        # Figure title (including species and sample count)
        fig.suptitle(f"SV Distribution Across Top {top_chroms} Chromosomes (Species: {species}, {n_samples} Samples)", 
                     y=1.02, fontsize=16)
        
        # Save as PDF, filename retains original format "chrom_sv_by_sample_top10" with species identifier
        output_path = os.path.join(output_dir, f"chrom_sv_by_sample_top{top_chroms}_{species}.pdf")
        plt.tight_layout()
        plt.savefig(output_path, format="pdf", bbox_inches="tight")
        plt.close()
        print(f"📈 Saved chromosome SV plot (PDF) for {species}: {output_path}")

def plot_sv_length_distribution(combined_df: pd.DataFrame, output_dir: str):
    """Plot distribution of SV lengths by SV type"""
    # Remove zero-length SVs for log scale
    non_zero_df = combined_df[combined_df["length"] > 0]
    
    if len(non_zero_df) == 0:
        print("⚠️ No non-zero length SVs found for length distribution plot")
        return
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Get unique SV types
    sv_types = non_zero_df["svtype"].unique()
    colors = ["#A23B72", "#F18F01", "#C73E1D", "#2E86AB", "#8E4EC6", "#3D9970", "#FF6B6B"]
    
    # Plot histogram for each SV type
    for i, sv_type in enumerate(sv_types):
        data = non_zero_df[non_zero_df["svtype"] == sv_type]["length"]
        if len(data) > 0:
            ax.hist(data, bins=np.logspace(np.log10(min(data)), np.log10(max(data)), 50), 
                   alpha=0.7, label=sv_type, color=colors[i % len(colors)], edgecolor='black', linewidth=0.2)
    
    ax.set_xscale('log')
    ax.set_xlabel("SV Length (bp)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of SV Lengths by SV Type")
    ax.legend(title="SV Type")
    ax.grid(alpha=0.3)
    
    # Save as PDF
    output_path = os.path.join(output_dir, "sv_length_distribution.pdf")
    plt.tight_layout()
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close()
    print(f"📈 Saved SV length distribution plot (PDF): {output_path}")

# --------------------------
# 4. Main Workflow (Multiple Samples)
# --------------------------
def main(input_path: str, output_dir: str = "delly_multi_sample_results", top_chroms: int = 10):
    """
    Main workflow:
    1. Load sample-VCF mapping
    2. Batch parse VCFs (annotate species, only for chromosome SV plot)
    3. Generate visualizations (only chromosome SV plot separated by species; all plots saved as PDF)
    4. Save results (data files only, no species separation)
    """
    # Initialize plot style and output directory
    set_plot_style()
    os.makedirs(output_dir, exist_ok=True)
    print(f"===== Starting Delly Multiple-Sample Analysis =====")
    print(f"Input: {input_path}")
    print(f"Output Directory: {os.path.abspath(output_dir)}")
    print(f"Top Chromosomes to Plot: {top_chroms}")
    print(f"Settings: Only chromosome SV plot split by species; All plots saved as PDF\n")
    
    # Step 1: Get sample-VCF mapping
    sample_vcf_dict = get_sample_vcfs(input_path)
    
    # Step 2: Batch parse VCFs (annotate species)
    combined_sv_df = batch_parse_vcfs(sample_vcf_dict)
    
    # Step 3: Generate visualizations (only chromosome SV plot separated by species, all plots saved as PDF)
    print(f"\n===== Generating Visualizations (All PDF) =====")
    plot_total_svs_per_sample(combined_df=combined_sv_df, output_dir=output_dir)  # No species separation
    plot_sv_type_by_sample(combined_df=combined_sv_df, output_dir=output_dir)      # No species separation
    plot_chrom_sv_by_species(combined_df=combined_sv_df, output_dir=output_dir, top_chroms=top_chroms)  # Separated by species
    plot_sv_length_distribution(combined_df=combined_sv_df, output_dir=output_dir)  # New plot: SV length distribution
    
    # Step 4: Save result data (no species separation, maintain original structure)
    print(f"\n===== Saving Results (No Species Split) =====")
    # Save combined CSV (includes species column for traceability)
    combined_csv = os.path.join(output_dir, "combined_sv_analysis.csv")
    combined_sv_df.to_csv(combined_csv, index=False, sep='\t')
    print(f"💾 Saved combined cohort data: {combined_csv}")
    
    # Save sample-specific CSVs (not separated by species, all stored in sample_specific)
    sample_dir = os.path.join(output_dir, "sample_specific")
    os.makedirs(sample_dir, exist_ok=True)
    for sample_id in combined_sv_df["sample_id"].unique():
        sample_df = combined_sv_df[combined_sv_df["sample_id"] == sample_id]
        sample_csv = os.path.join(sample_dir, f"{sample_id}_sv_results.csv")
        sample_df.to_csv(sample_csv, index=False, sep='\t')
    print(f"💾 Saved sample-specific CSVs to: {sample_dir}")
    
    print(f"\n===== Analysis Completed Successfully =====")
    return combined_sv_df

# --------------------------
# 5. Execution Entry (RUN THIS!)
# --------------------------
if __name__ == "__main__":
    # --------------------------
    # USER CONFIGURATION (EDIT THIS SECTION!)
    # --------------------------
    # Option 1: VCF folder path (recommended)
    INPUT_PATH = "/delly/sv_vcf"  # Replace with your VCF folder path
    
    # Option 2: Sample-VCF list file (format: sample_id\tvcf_path)
    # INPUT_PATH = "/delly/sample_list.txt"
    
    # Output directory (automatically created)
    OUTPUT_DIR = "delly_multi_sample_results_pdf"
    
    # Number of TOP chromosomes to show in chromosome SV plot (fixed at 10, matching "chrom_sv_by_sample_top10")
    TOP_CHROMOSOMES = 10
    
    # --------------------------
    # Execute analysis
    # --------------------------
    try:
        combined_results = main(
            input_path=INPUT_PATH,
            output_dir=OUTPUT_DIR,
            top_chroms=TOP_CHROMOSOMES
        )
        
        # Print final statistical summary
        print("\n===== Final Analysis Summary =====")
        print(f"Total samples analyzed: {combined_results['sample_id'].nunique()}")
        print(f"Total SVs detected: {len(combined_results)}")
        print(f"Chromosome SV plot split by species: {sorted(combined_results['species'].unique())}")
        print(f"All plots saved as PDF in: {os.path.abspath(OUTPUT_DIR)}")
        
    except Exception as e:
        print(f"\n❌ Analysis failed: {str(e)}")
        import traceback
        traceback.print_exc()
        import sys
        sys.exit(1)

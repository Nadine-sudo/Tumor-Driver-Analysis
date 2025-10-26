#!/bin/bash

# Delly somatic structural variant calling pipeline
# This script processes tumor-normal paired BAM files using Delly

# Set variables
REFERENCE="reference.fasta"
EXCLUDE_REGIONS="human.hg38.excl.tsv"
SAMPLES_FILE="samples.tsv"
INPUT_DIR="."
OUTPUT_DIR="./sv_vcf"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if required files exist
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference genome file '$REFERENCE' not found!"
    exit 1
fi

if [ ! -f "$EXCLUDE_REGIONS" ]; then
    echo "Error: Exclude regions file '$EXCLUDE_REGIONS' not found!"
    exit 1
fi

# Check if tumor BAM files exist
if [ ! -d "./tumor_bams" ] || [ -z "$(ls -A ./tumor_bams/*.bam 2>/dev/null)" ]; then
    echo "Error: No tumor BAM files found in ./tumor_bams/"
    exit 1
fi

# Check if normal BAM files exist
if [ ! -d "./normal_bams" ] || [ -z "$(ls -A ./normal_bams/*.bam 2>/dev/null)" ]; then
    echo "Error: No normal BAM files found in ./normal_bams/"
    exit 1
fi

echo "Starting Delly somatic SV analysis..."
echo "Reference genome: $REFERENCE"
echo "Exclude regions: $EXCLUDE_REGIONS"
echo "Samples file: $SAMPLES_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "----------------------------------------"

# Process each tumor BAM file
for tumor_bam in ./tumor_bams/*.bam; do
    # Extract sample name from BAM file name
    sample_name=$(basename "$tumor_bam" .bam)
    echo "Processing sample: $sample_name"
    
    # Define corresponding normal BAM
    normal_bam="./normal_bams/$(basename "$tumor_bam")"
    
    # Check if normal BAM exists
    if [ ! -f "$normal_bam" ]; then
        echo "Warning: Normal BAM '$normal_bam' not found for tumor '$tumor_bam'. Skipping..."
        continue
    fi
    
    # Define output files
    bcf_file="$OUTPUT_DIR/${sample_name}.bcf"
    filtered_bcf="$OUTPUT_DIR/${sample_name}.pre.bcf"
    vcf_file="$OUTPUT_DIR/${sample_name}.vcf"
    
    echo "  Running Delly call..."
    # Run Delly call
    delly call \
        -g "$REFERENCE" \
        -x "$EXCLUDE_REGIONS" \
        -o "$bcf_file" \
        "$tumor_bam" "$normal_bam"
    
    # Check if Delly call was successful
    if [ $? -ne 0 ] || [ ! -f "$bcf_file" ]; then
        echo "Error: Delly call failed for sample $sample_name"
        continue
    fi
    
    echo "  Running Delly filter..."
    # Run Delly filter
    delly filter \
        -f somatic \
        -o "$filtered_bcf" \
        -s "$SAMPLES_FILE" \
        "$bcf_file"
    
    # Check if Delly filter was successful
    if [ $? -ne 0 ] || [ ! -f "$filtered_bcf" ]; then
        echo "Error: Delly filter failed for sample $sample_name"
        continue
    fi
    
    echo "  Converting to VCF..."
    # Convert BCF to VCF
    bcftools view "$filtered_bcf" -O v -o "$vcf_file"
    
    # Check if conversion was successful
    if [ $? -ne 0 ] || [ ! -f "$vcf_file" ]; then
        echo "Error: BCF to VCF conversion failed for sample $sample_name"
        continue
    fi
    
    echo "  Completed sample: $sample_name"
    echo "  ----------------------------------------"
done

echo "Delly analysis completed!"
echo "Running SV plot script..."

# Check if any VCF files were generated
if [ ! -d "$OUTPUT_DIR" ] || [ -z "$(ls -A $OUTPUT_DIR/*.vcf 2>/dev/null)" ]; then
    echo "Warning: No VCF files found in $OUTPUT_DIR. Skipping plot generation."
else
    # Run the Python plotting script
    if command -v python3 &> /dev/null; then
        python3 sv_plot.py
    else
        echo "Error: Python not found. Please install Python to generate plots."
        exit 1
    fi
fi

echo "Pipeline finished!"
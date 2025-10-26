#!/bin/bash

# Pipeline for calculating mappability (optional) and running CNV analysis
# Usage: ./run.sh [-m] [-r reference.fa] [-b bam_pattern] [-o output_dir]
#   -m: Run mappability calculation (optional)
#   -r: Reference genome 
#   -b: BAM file pattern 
#   -o: Output directory (default: cnv)
#   -h: Show help

# Default values
RUN_MAPPABILITY=false
REFERENCE="hg38.fa"
BAM_PATTERN="bam_file/*.bam"
OUTPUT_DIR="cnv"
MAPPABILITY_FILE="map.fa.gz"

# Parse command line arguments
while getopts "mr:b:o:h" opt; do
    case $opt in
        m)
            RUN_MAPPABILITY=true
            ;;
        r)
            REFERENCE="$OPTARG"
            ;;
        b)
            BAM_PATTERN="$OPTARG"
            ;;
        o)
            OUTPUT_DIR="$OPTARG"
            ;;
        h)
            echo "Usage: $0 [-m] [-r reference.fa] [-b bam_pattern] [-o output_dir]"
            echo "  -m: Run mappability calculation (optional)"
            echo "  -r: Reference genome (default: $REFERENCE)"
            echo "  -b: BAM file pattern (default: $BAM_PATTERN)"
            echo "  -o: Output directory (default: $OUTPUT_DIR)"
            echo "  -h: Show help"
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done

# Function to check if required tools are available
check_tool() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: $1 is not installed or not in PATH"
        exit 1
    fi
}

# Function to check if file exists
check_file() {
    if [ ! -f "$1" ]; then
        echo "Error: $1 not found!"
        exit 1
    fi
}

# Step 1: Calculate mappability using yeast genome (optional)
if [ "$RUN_MAPPABILITY" = true ]; then
    echo "=== Step 1: Calculating mappability with yeast genome (OPTIONAL) ==="
    
    # Check for required files
    check_file "sacCer3.fa"
    
    # Check for required tools
    check_tool "dicey"
    check_tool "bwa"
    check_tool "samtools"
    
    echo "Chopping genome for Dicey..."
    dicey chop sacCer3.fa
    
    echo "Indexing genome with BWA..."
    bwa index sacCer3.fa
    
    echo "Aligning reads..."
    bwa mem sacCer3.fa read1.fq.gz read2.fq.gz | samtools sort -@ 4 -o srt.bam -
    
    echo "Indexing BAM file..."
    samtools index srt.bam 
    
    echo "Calculating mappability..."
    dicey mappability2 srt.bam 
    
    echo "Compressing and indexing mappability file..."
    if [ -f "map.fa.gz" ]; then
        gunzip map.fa.gz && bgzip map.fa && samtools faidx map.fa.gz
        MAPPABILITY_FILE="map.fa.gz"
    else
        echo "Warning: map.fa.gz not generated. Using default mappability file if available."
    fi
    
    # Clean up intermediate files
    rm -f srt.bam srt.bam.bai sacCer3.fa.sa sacCer3.fa.pac sacCer3.fa.bwt sacCer3.fa.ann sacCer3.fa.amb read1.fq.gz read2.fq.gz
else
    echo "=== Skipping Step 1: Mappability calculation ==="
    echo "Using existing mappability file: $MAPPABILITY_FILE"
fi

# Check if mappability file exists
if [ ! -f "$MAPPABILITY_FILE" ]; then
    echo "Warning: Mappability file $MAPPABILITY_FILE not found!"
    echo "Proceeding without mappability correction..."
    MAPPABILITY_PARAM=""
else
    MAPPABILITY_PARAM="-m $MAPPABILITY_FILE"
    echo "Using mappability file: $MAPPABILITY_FILE"
fi

# Step 2: CNV analysis on samples
echo "=== Step 2: Running CNV analysis on samples ==="

# Check if reference genome exists
check_file "$REFERENCE"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Find BAM files
echo "Searching for BAM files matching pattern: $BAM_PATTERN"
shopt -s nullglob
bam_files=($BAM_PATTERN)

if [ ${#bam_files[@]} -eq 0 ]; then
    echo "Warning: No BAM files found matching pattern: $BAM_PATTERN"
    echo "Skipping CNV analysis..."
else
    echo "Found ${#bam_files[@]} BAM files"
    
    # Check for required tools
    check_tool "delly"
    
    # Process each BAM file
    for file in "${bam_files[@]}"; do
        filename=$(basename "$file")
        prefix="${filename%.*}"  # Remove extension
        
        echo "Processing sample: $prefix"
        
        # Run Delly CNV calling
        echo "  Running Delly CNV calling..."
        if [ -n "$MAPPABILITY_PARAM" ]; then
            delly cnv -a -g "$REFERENCE" -m $MAPPABILITY_PARAM -c "$OUTPUT_DIR/${prefix}.cov.gz" -o out.bcf "$file"
        else
            echo "Error: mappability file not found!"
        fi
        
        # Check if Delly completed successfully
        if [ ! -f "out.bcf" ]; then
            echo "  Warning: Delly failed for sample $prefix. Skipping..."
            continue
        fi
        
        # Filter and convert results
        echo "  Filtering and converting results..."
        
        bcftools filter -i 'QUAL >= 20' -O b -o filtered.cnv.bcf out.bcf
        bcftools view filtered.cnv.bcf -O v -o filtered.cnv.vcf
           
        # Move final results to output directory
        mv filtered.cnv.vcf "$OUTPUT_DIR/${prefix}.vcf"
        
        # Clean up temporary files
        rm -f out.bcf filtered.cnv.bcf
        
        echo "  Completed sample: $prefix"
    done
fi

# Step 3: Generate visualizations
echo "=== Step 3: Generating visualizations ==="

# Check if CNV files were generated
if [ ! -d "$OUTPUT_DIR" ] || [ -z "$(ls -A "$OUTPUT_DIR"/*.vcf 2>/dev/null)" ]; then
    echo "Warning: No CNV VCF files found in $OUTPUT_DIR/ directory. Skipping visualization..."
else
    # Check if R script exists
    if [ ! -f "cnv_plot.R" ]; then
        echo "Warning: cnv_plot.R not found. Skipping visualization..."
    else
        echo "Running R script for visualization..."
        if command -v Rscript &> /dev/null; then
            # Update input directory in R script if needed
            sed -i "s|input_dir <- \".*\"|input_dir <- \"./$OUTPUT_DIR\"|" cnv_plot.R 2>/dev/null || true
            Rscript cnv_plot.R
        else
            echo "Warning: Rscript not found. Please install R to generate visualizations."
        fi
    fi
fi

echo "Pipeline completed!"

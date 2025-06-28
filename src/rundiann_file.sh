#!/bin/bash

# Initialize variables
inputs=()
library=""
outputdir=""
output="result_peptide.tsv"
extra_params=()
diann_path="/usr/diann/1.8.1/diann-1.8.1"  # Default path

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dir)
            inputs+=("--dir $2")
            shift 2
            ;;
        --f)
            inputs+=("--f $2")
            shift 2
            ;;
        --lib)
            library="$2"
            shift 2
            ;;
        --out)
            output="$2"
            shift 2
            ;;
        --output-dir)  # Parameter to explicitly set output directory
            outputdir="$2"
            shift 2
            ;;
        --diann-path)  # NEW: Parameter for DIANN path
            diann_path="$2"
            shift 2
            ;;
        --threads|--verbose|--qvalue|--matrix-qvalue|--mass-acc|--mass-acc-ms1|--double-search|--no-prot-inf|--rt-profiling|--pg-level|--report-lib-info|--matrices|--reanalyse)
            # Store any additional parameters with their values
            param_name="$1"
            param_value="$2"
            extra_params+=("$param_name $param_value")
            shift 2
            ;;
        *)
            # For unrecognized parameters, just add them directly
            extra_params+=("$1")
            shift
            ;;
    esac
done

# Check if output directory is set
if [ -z "$outputdir" ]; then
    echo "Error: Output directory not specified. Please use --output-dir parameter."
    exit 1
fi

# Check if output directory exists
if [ ! -d "$outputdir" ]; then
    echo "Creating output directory: $outputdir"
    mkdir -p "$outputdir"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create output directory: $outputdir"
        exit 1
    fi
fi

# Change to output directory
cd "$outputdir" || { echo "Error: Cannot change to directory $outputdir"; exit 1; }

# Build command using the specified DIANN path
diann_cmd="$diann_path --lib $library --out $output"

# Add inputs
for input in "${inputs[@]}"; do
    diann_cmd="$diann_cmd $input"
done

# Add all extra parameters
for param in "${extra_params[@]}"; do
    diann_cmd="$diann_cmd $param"
done

# Execute command
echo "Executing: $diann_cmd"
eval "$diann_cmd"
#!/bin/bash
# Usage:
#   ./make_extend.sh --mdname <mdname_name> --extend <extend_time_ns>
# Example:
   ./make_extend.sh --mdname md_100 --extend 200

# Embed completion section
if [[ $1 == --complete ]]; then
    echo "-mdname -extend -h"
    exit 0
fi

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --mdname)
            mdname="$2"
            shift 2
            ;;
        --extend)
            extend_ns="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 --mdname <mdname_name> --extend <extend_time_ns>"
            echo "Example: $0 --mdname md_100 --extend 200"
            exit 0
            ;;
        *)
            echo "Unknown parameter: $1"
            echo "Usage: $0 --mdname <mdname_name> --extend <extend_time_ns>"
            exit 1
            ;;
    esac
done

# Check input
if [ -z "$mdname" ] || [ -z "$extend_ns" ]; then
    echo "Error: missing required arguments."
    echo "Usage: $0 --mdname <mdname_name> --extend <extend_time_ns>"
    exit 1
fi

# Convert ns to ps for gmx convert-tpr
extend_ps=$(echo "$extend_ns * 1000" | bc)

# Derive new TPR name
new_tpr="${mdname%_*}_extend.tpr"

# Step 1: Extend the TPR file
gmx convert-tpr -s ${mdname}.tpr -extend ${extend_ps} -o ${new_tpr}
if [ $? -ne 0 ]; then
    echo "Error: gmx convert-tpr failed"
    exit 1
fi

# Step 2: Continue the simulation (appending to existing files)
gmx mdrun -deffnm ${mdname} -s ${new_tpr} -cpi ${mdname}.cpt -v


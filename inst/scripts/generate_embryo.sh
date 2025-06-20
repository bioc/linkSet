
#!/bin/bash

# This script downloads the mouse Embryo body enhancer BED file from EnhancerAtlas
# Source: http://www.enhanceratlas.org/data/download/enhancer/mm/Embryo_body.bed

# Define the URL and output file name
URL="http://www.enhanceratlas.org/data/download/enhancer/mm/Embryo_body.bed"
OUTPUT="Embryo_body.bed"

# Check if the file already exists
if [ -f "$OUTPUT" ]; then
    echo "File $OUTPUT already exists. Skipping download."
else
    echo "Downloading Embryo body enhancer BED file..."
    curl -fSL "$URL" -o "$OUTPUT"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to download $URL"
        exit 1
    fi
    echo "Download complete! File saved as $OUTPUT"
fi

#!/bin/bash

# Get the input file name from the first argument
input_file="$1"

# Extract the base name (without extension) using parameter expansion
base_name="${input_file%.*}"

# Specify the new extension (for example, .newext)

# Form the output file name by concatenating the base name with the new extension
output_file="${base_name}".nc

# Now use the output file name in your command
# Example command: replace `command` with your actual program
#gdal_translate -of NetCDF "$input_file" "$output_file"

gdalwarp $input_file $output_file -t_srs "+proj=longlat +ellps=WGS84"


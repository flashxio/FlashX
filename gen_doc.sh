#!/bin/bash

echo 'roxygen2::roxygenize("Rpkg/")' | R --no-save
mkdir -p docs/FlashGraphR
for file in `ls flash-graph/Rpkg/man/`
do
	input_file="Rpkg/man/$file"
	output_file="docs/FlashGraphR/${file}.txt"
	echo "tools::Rd2txt(\"$input_file\", out=\"$output_file\", options=list(underline_titles=FALSE))" | R --no-save
done

#!/bin/bash

echo 'roxygen2::roxygenize("flash-graph/Rpkg/")' | R --no-save
mkdir -p doc/FlashGraphR
for file in `ls flash-graph/Rpkg/man/`
do
	R CMD Rd2txt "flash-graph/Rpkg/man/$file" -o "doc/FlashGraphR/$file"
done

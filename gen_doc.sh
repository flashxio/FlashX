#!/bin/bash

echo 'roxygen2::roxygenize("Rpkg/")' | R --no-save
echo 'roxygen2::roxygenize("FlashGraphR/")' | R --no-save

#!/bin/sh

R CMD build FlashR-learn
version=`R --version | grep "R version" | awk '{print $3}' | awk -F. '{print $1 "." $2}'`

R CMD INSTALL --no-configure FlashRLearn_0.1-0.tar.gz
rm FlashRLearn_0.1-0.tar.gz

#!/bin/bash

Rscript --vanilla -e "source('_targets.R'); tar_make_clustermq(workers = 16)" > output.txt


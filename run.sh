#!/bin/bash

Rscript --vanilla -e "source('_targets.R'); tar_make()" > output.txt


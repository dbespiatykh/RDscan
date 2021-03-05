#!/usr/bin/env RScript --vanilla

library(tidyverse)
library(openxlsx)
library(optparse)
 
option_list <- list(
  make_option(c("-T", "--threshold"), type = "double", default = NULL, 
              help = "Input treshold value"),
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Input SV table"),
  make_option(c("-x", "--xlsx"), type = "character", default = NULL, 
              help = "Input .xlsx filename"),
  make_option(c("-t", "--tsv"), type = "character", default = NULL, 
              help = "Input .tsv filename"),
  make_option(c("-b", "--bin"), type = "character", default = NULL, 
              help = "Input .bin.tsv filename")
)
 
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$threshold)){
  print_help(opt_parser)
  stop("Treshold value should be supplied", call.=FALSE)
}

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Input table should be supplied", call.=FALSE)
}

if (is.null(opt$xlsx)){
  print_help(opt_parser)
  stop("Output .xlsx filename should be supplied", call.=FALSE)
}

if (is.null(opt$tsv)){
  print_help(opt_parser)
  stop("Output .tsv filename should be supplied", call.=FALSE)
}

if (is.null(opt$bin)){
  print_help(opt_parser)
  stop("Output .bin.tsv filename should be supplied", call.=FALSE)
}

threshold <- opt$threshold
input <- opt$input
xlsx_output <- opt$xlsx
tab_output <- opt$tsv
bin_output <- opt$bin

posStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
rule <- sprintf("<=%s", threshold)

long <- read.table(input, header = TRUE, check.names = FALSE, sep = "\t")
wide <- spread(long, key = RD, value = depth)
binary <- wide %>% mutate_if(is.numeric, ~1 * (. <= threshold))
wb <- openxlsx::createWorkbook("RD_results")
openxlsx::addWorksheet(wb, "RD", gridLines = TRUE)
writeData(wb, sheet = 1, wide)
conditionalFormatting(wb = wb,
                      sheet = 'RD',
                      cols = 1:ncol(wide),
                      rows = 1:nrow(wide)+1,
                      rule = rule,
                      type = "expression",
                      style = posStyle
)
suppressMessages(openxlsx::saveWorkbook(wb, xlsx_output, overwrite = TRUE))
write.table(wide, tab_output, sep = "\t", row.names = FALSE)
write.table(binary, bin_output, sep = "\t", row.names = FALSE)
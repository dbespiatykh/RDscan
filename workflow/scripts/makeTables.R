#!/usr/bin/env RScript --vanilla

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

threshold <- snakemake@params[["threshold"]]
input <- snakemake@input[[1]]
xlsx_output <- snakemake@output[[1]]
tab_output <- snakemake@output[[2]]
bin_output <- snakemake@output[[3]]

posStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
rule <- sprintf("<=%s", threshold)

long <-
  read.table(input,
             header = TRUE,
             check.names = FALSE,
             sep = "\t")
wide <- spread(long, key = RD, value = depth)
binary <- wide %>% mutate_if(is.numeric, ~ 1 * (. <= threshold))
wb <- openxlsx::createWorkbook("RD_results")
openxlsx::addWorksheet(wb, "RD", gridLines = TRUE)
writeData(wb, sheet = 1, wide)
conditionalFormatting(
  wb = wb,
  sheet = 'RD',
  cols = 1:ncol(wide),
  rows = 1:nrow(wide) + 1,
  rule = rule,
  type = "expression",
  style = posStyle
)
suppressMessages(openxlsx::saveWorkbook(wb, xlsx_output, overwrite = TRUE))
write.table(wide, tab_output, sep = "\t", row.names = FALSE)
write.table(binary, bin_output, sep = "\t", row.names = FALSE)

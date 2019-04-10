#!/usr/bin/R
# USAGE: Rscript create_kraken_mqc_file.R
# compute mqc file for multiQC reports from kraken genus report
file <- "mini_kraken/genus_map_report_filtered.tsv"
t <- read.delim(file,sep="\t",header=T,comment.char = "#",colClasses="character",row.names=1, fill=TRUE)

rows <- length(t[,1])

# compute log10 and set negative (-Inf) values to zero
for (i in 1:rows) {
  row_log10 <- log10(as.numeric(t[i,]))
  row_log10[row_log10 < 0] <- 0
  t[i,] <- row_log10

}

# write mqc report (for multiQC)
mqc_report = "mini_kraken_summary/kraken_genus_map_report_mqc.txt"
title <- "# title: 'Kraken report (genus level)'"
desc <- "# description: 'kraken-mpa-report --db $minikrakenDB mini_kraken/*sequences.kraken >mini_kraken/mpa_report.tsv'"
section <- "# section: 'Custom Data File'"
format <- "# format: 'tsv'"
plot_type <- "# plot_type: 'heatmap'"
pconfig <- "# pconfig:"
id <- "#    id: 'samples'"
ylab <- "#    ylab: 'clades'"

#writeLines(title, mqc_report)
#writeLines(desc, mqc_report)
#writeLines(section, mqc_report)
#writeLines(format, mqc_report)
#writeLines(plot_type, mqc_report)
#writeLines(pconfig, mqc_report)
#writeLines(id, mqc_report)
#writeLines(ylab, mqc_report)
column_names <- c("clade", colnames(t))
column_names <- gsub('\\.', '-', column_names)
s_column_names <- paste(column_names,collapse="\t")

writeLines(c(title, desc, section, format, plot_type, pconfig, id, ylab, s_column_names), mqc_report)
#colnames(t) <- c("clade",colnames(t))
write.table(t, file=mqc_report, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE, append=TRUE, qmethod = c("escape"))

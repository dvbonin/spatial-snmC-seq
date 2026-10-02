library(optparse)
library(data.table)

# Merged BED of the RepeatMasker classes that collapse many near-identical copies onto
# one reference position: satellite arrays, simple repeats, low-complexity stretches and
# rDNA. SINE/LINE/LTR are left in -- they cover roughly half the genome and are mostly
# mappable, so masking them would cost far more CpGs than the pileups are worth.
#
# rmsk.txt.gz is the UCSC table: 6 genoName, 7 genoStart (0-based), 8 genoEnd, 12 repClass.

option_list <- list(
  make_option("--rmsk", type = "character"),
  make_option("--classes", type = "character",
              default = "Satellite,Simple_repeat,Low_complexity,rRNA"),
  make_option("--output", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

cl <- strsplit(opt$classes, ",")[[1]]
d <- fread(cmd = sprintf("zcat %s", opt$rmsk), select = c(6, 7, 8, 12),
           col.names = c("chrom", "start", "end", "repClass"))
d <- d[repClass %in% cl]
setorder(d, chrom, start)

# merge overlapping and abutting intervals per chromosome
d[, grp := cumsum(start > shift(cummax(end), fill = -1L)), by = chrom]
m <- d[, .(start = min(start), end = max(end)), by = .(chrom, grp)][, grp := NULL]

fwrite(m, opt$output, sep = "\t", col.names = FALSE)
message(sprintf("%s intervals in classes {%s} -> %s merged, %.1f Mb",
                format(nrow(d), big.mark = ","), opt$classes,
                format(nrow(m), big.mark = ","), sum(as.numeric(m$end - m$start)) / 1e6))

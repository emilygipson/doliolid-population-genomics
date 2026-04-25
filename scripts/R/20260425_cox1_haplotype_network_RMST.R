# 20260425_cox1_haplotype_network.R
# rmst, mitobim 119

library(ape)
library(pegas)

fasta_file  <- "C:\\Users\\emily\\OneDrive - University of Georgia\\----------manuscripts----------\\genetics\\mitochondrial\\working directories\\4.25.26 workdir\\mitobim_unfilt_cox1_aligned.fasta"
meta_file   <- "C:\\Users\\emily\\OneDrive - University of Georgia\\----------manuscripts----------\\genetics\\mitochondrial\\working directories\\4.25.26 workdir\\COORDINATE-CORRECTD-133_samplelist-metadata_wFieldData_plus_EG16_20260423.csv"
id_map_file <- "C:\\Users\\emily\\OneDrive - University of Georgia\\----------manuscripts----------\\genetics\\mitochondrial\\working directories\\4.25.26 workdir\\CLAUDE_eg_sample_id_map_20260424.tsv"
out_dir     <- "C:\\Users\\emily\\OneDrive - University of Georgia\\----------manuscripts----------\\genetics\\mitochondrial\\working directories\\4.25.26 workdir"

aln <- read.FASTA(fasta_file, type = "DNA")

# translate EG labels to real sample IDs
id_map <- read.delim(id_map_file, stringsAsFactors = FALSE)
hit <- match(names(aln), id_map$alignment_label)
new_names <- ifelse(is.na(hit), names(aln), id_map$real_sample_id[hit])
eg_unmapped <- new_names[grepl("^EG", new_names) & new_names == names(aln)]
if (length(eg_unmapped) > 0) stop("EG labels not in map: ", paste(eg_unmapped, collapse = ", "))
names(aln) <- new_names

# drop divergent sample if present
if ("DM18_27" %in% names(aln)) aln <- aln[names(aln) != "DM18_27"]

# strict ACGT column filter
aln_mat <- as.character(as.matrix(aln))
keep_col <- apply(aln_mat, 2, function(x) all(tolower(x) %in% c("a", "c", "g", "t")))
aln_dna <- as.DNAbin(aln_mat[, keep_col, drop = FALSE])

# metadata + year column
meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
names(meta)[names(meta) == "sample ID"] <- "sample_id"
meta$year <- vapply(strsplit(meta$collect_date, "/"), function(x) x[3], character(1))
meta_match <- meta[match(rownames(aln_mat[, keep_col, drop = FALSE]), meta$sample_id), ]
if (any(is.na(meta_match$sample_id))) stop("metadata mismatch")

# rmst
h <- haplotype(aln_dna)
d <- dist.dna(h, model = "N", pairwise.deletion = TRUE)
net <- rmst(d, quiet = TRUE)
hap_freq <- vapply(attr(h, "index"), length, integer(1))

# pie table by year
make_pie <- function(h, group_vec) {
  rows <- do.call(rbind, lapply(seq_along(attr(h, "index")), function(i) {
    data.frame(hap = rownames(h)[i],
               group = group_vec[attr(h, "index")[[i]]],
               stringsAsFactors = FALSE)
  }))
  tab <- table(rows$hap, rows$group)
  tab[rownames(h), , drop = FALSE]
}
pie_year <- make_pie(h, meta_match$year)

col_year <- c("2016" = "#E41A1C",
              "2017" = "#377EB8",
              "2022" = "#4DAF4A",
              "2024" = "#984EA3")
pie_year <- pie_year[, names(col_year), drop = FALSE]

png(file.path(out_dir, "20260425_cox1_rmst_n119_year.png"),
    width = 8, height = 8, units = "in", res = 300)
par(mar = c(1, 1, 3, 1))
plot(net,
     size = hap_freq, pie = pie_year, bg = col_year,
     legend = FALSE, fast = FALSE, show.mutation = 2,
     scale.ratio = 1, labels = FALSE, cex = 0.7)
legend("topright", legend = names(col_year), fill = col_year, bty = "n", cex = 0.9)
title(main = "cox1 RMST, n=119 (strict ACGT), year")
dev.off()

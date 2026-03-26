# 20260319_mitogenome_map.R
# Circular mitogenome map for D. gegenbauri (DD_21_05, 15,329 bp)
# For lab meeting presentation

# ---- Parse GFF ----
gff <- read.delim("Galaxy6-[MITOS2 on dataset 4_ GFF].gff",
                  header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

# Extract gene name from attributes
get_name <- function(attr) {
    m <- regmatches(attr, regexpr("Name=[^;]+", attr))
    if(length(m) == 0) return(NA)
    sub("Name=", "", m)
}
gff$name <- sapply(gff$attributes, get_name)

# Get genes (type == "gene") and tRNAs/rRNAs from ncRNA_gene
pcg <- gff[gff$source == "mitos" & gff$type == "gene", ]
trna_rrna <- gff[gff$source == "mitfi" & gff$type == "ncRNA_gene", ]

# Classify trna vs rrna
trna_rrna$genetype <- ifelse(grepl("rrn", trna_rrna$name), "rRNA", "tRNA")

# Combine
genes <- data.frame(
    name  = pcg$name,
    start = pcg$start,
    end   = pcg$end,
    type  = "PCG",
    stringsAsFactors = FALSE
)

noncoding <- data.frame(
    name  = trna_rrna$name,
    start = trna_rrna$start,
    end   = trna_rrna$end,
    type  = trna_rrna$genetype,
    stringsAsFactors = FALSE
)

all_genes <- rbind(genes, noncoding)

# ---- Genome parameters ----
genome_len <- 15329

# Convert position to angle (radians, starting at top, clockwise)
pos_to_angle <- function(pos) {
    (pi/2) - (pos / genome_len) * 2 * pi
}

# Draw an arc (filled polygon between two radii)
draw_arc <- function(start, end, r_inner, r_outer, col, border = "white", lwd = 0.5) {
    a1 <- pos_to_angle(start)
    a2 <- pos_to_angle(end)
    # Generate sequence of angles (clockwise means decreasing angle)
    angles <- seq(a1, a2, length.out = max(50, abs(end - start) / 5))
    x_outer <- r_outer * cos(angles)
    y_outer <- r_outer * sin(angles)
    x_inner <- r_inner * cos(rev(angles))
    y_inner <- r_inner * sin(rev(angles))
    polygon(c(x_outer, x_inner), c(y_outer, y_inner), col = col, border = border, lwd = lwd)
}

# ---- Colors ----
col_pcg  <- "#4A90D9"   # blue
col_rrna <- "#E8A838"   # gold/orange
col_trna <- "#5AAA46"   # green
col_bg   <- "gray92"    # background ring

# ---- Plot ----
png("mitogenome_map.png", width = 8, height = 8, units = "in", res = 300)
par(mar = c(0, 0, 0, 0))
plot(NA, xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), asp = 1,
     axes = FALSE, xlab = "", ylab = "")

# Background ring
r_inner_gene <- 0.75
r_outer_gene <- 1.05
r_inner_trna <- 1.06
r_outer_trna <- 1.15

# Draw background ring for genes
angles_full <- seq(0, 2 * pi, length.out = 500)
polygon(r_outer_gene * cos(angles_full), r_outer_gene * sin(angles_full),
        col = col_bg, border = "gray70", lwd = 0.5)
polygon(r_inner_gene * cos(angles_full), r_inner_gene * sin(angles_full),
        col = "white", border = "gray70", lwd = 0.5)

# Draw PCGs and rRNAs on main ring
for(i in 1:nrow(all_genes)) {
    g <- all_genes[i, ]
    if(g$type == "PCG") {
        draw_arc(g$start, g$end, r_inner_gene, r_outer_gene, col_pcg, border = "white", lwd = 0.8)
    } else if(g$type == "rRNA") {
        draw_arc(g$start, g$end, r_inner_gene, r_outer_gene, col_rrna, border = "white", lwd = 0.8)
    } else if(g$type == "tRNA") {
        draw_arc(g$start, g$end, r_inner_trna, r_outer_trna, col_trna, border = "white", lwd = 0.5)
    }
}

# ---- Gene labels ----
# Label PCGs and rRNAs at midpoint of arc
label_genes <- all_genes[all_genes$type %in% c("PCG", "rRNA"), ]
for(i in 1:nrow(label_genes)) {
    g <- label_genes[i, ]
    mid <- (g$start + g$end) / 2
    a <- pos_to_angle(mid)
    r_label <- (r_inner_gene + r_outer_gene) / 2
    
    # Determine rotation angle for text
    deg <- a * 180 / pi
    # Adjust so text reads correctly
    if(deg < -90 | deg > 90) {
        deg <- deg + 180
    }
    
    # Gene name formatting
    gname <- g$name
    
    # Font size based on gene length
    gene_len <- g$end - g$start
    cex_val <- ifelse(gene_len > 800, 0.7, ifelse(gene_len > 400, 0.6, 0.5))
    
    text(r_label * cos(a), r_label * sin(a), gname,
         srt = deg, cex = cex_val, col = "white", font = 2)
}

# Label tRNAs on outer edge
trna_genes <- all_genes[all_genes$type == "tRNA", ]
for(i in 1:nrow(trna_genes)) {
    g <- trna_genes[i, ]
    mid <- (g$start + g$end) / 2
    a <- pos_to_angle(mid)
    r_label <- r_outer_trna + 0.06
    
    deg <- a * 180 / pi
    if(deg < -90 | deg > 90) {
        deg <- deg + 180
    }
    
    # Clean up tRNA name: trnG_1 -> G, trnM_0 -> M
    short <- sub("trn", "", g$name)
    short <- sub("_[0-9]+", "", short)
    
    text(r_label * cos(a), r_label * sin(a), short,
         srt = deg, cex = 0.45, col = "gray30")
}

# ---- Tick marks (every 1 kb) ----
for(kb in seq(0, 15000, by = 1000)) {
    a <- pos_to_angle(kb)
    r1 <- r_outer_trna + 0.01
    r2 <- r_outer_trna + 0.05
    segments(r1 * cos(a), r1 * sin(a), r2 * cos(a), r2 * sin(a),
             col = "gray50", lwd = 0.5)
    # Label every 2 kb
    if(kb %% 2000 == 0 & kb > 0) {
        r3 <- r_outer_trna + 0.12
        text(r3 * cos(a), r3 * sin(a), paste0(kb/1000, "k"),
             cex = 0.45, col = "gray40")
    }
}

# ---- Center text ----
text(0, 0.08, expression(italic("D. gegenbauri")), cex = 1.5)
text(0, -0.08, "15,329 bp", cex = 1.1, col = "gray30")
text(0, -0.22, "13 PCGs, 2 rRNAs, 23 tRNAs", cex = 0.8, col = "gray50")

# ---- Legend ----
legend(-1.35, -0.9,
       legend = c("Protein-coding", "rRNA", "tRNA"),
       fill = c(col_pcg, col_rrna, col_trna),
       border = "white",
       cex = 0.85,
       bty = "n")

dev.off()
cat("Saved: mitogenome_map.png\n")

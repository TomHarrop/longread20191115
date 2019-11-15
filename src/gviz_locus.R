#!/usr/bin/env Rscript

library(Gviz)
library(GenomicRanges)

options(ucscChromosomeNames=FALSE)

# load the txdb
txdb_file <- "data/txdb.sqlite"
txdb <- AnnotationDbi::loadDb(txdb_file)

# set up transcripts, exons and coding sequences
tx <- transcriptsBy(txdb, "gene")
cds <- cdsBy(txdb, "gene")

# find Csd in the txdb
csd_transcripts <- tx$Csd
csd_cds <- cds$Csd
csd_cds$cds_name <- "complementary sex determiner"
csd_cds$id <- csd_cds$cds_name

# define ranges to plot
csd_chr <- as.character(seqnames(csd_transcripts))[[1]]
csd_start <- min(start(csd_cds))
csd_end <- max(end(csd_cds))

# set up track schemes
set1 <- viridis::viridis_pal()(4)

my_scheme <- list(
    shape = "smallArrow",
    background.title = "transparent",
    fontface.title = 1,
    col.title = set1[1],
    col.frame = set1[1],
    thinBoxFeature = c("lincRNA"),
    fontcolor.group = set1[1],
    cex.title = 1,
    cex.group = 0.5,
    showTitle=FALSE
)

grt_scheme <- as.list(c(
    my_scheme,
    col.frame = set1[2],
    fontcolor.group = set1[2],
    fontcolor.title = set1[2],
    size = 1,
    col = set1[2],
    col.line = set1[2],
    fill = set1[2],
    fontface.group = 3,
    transcriptAnnotation = "symbol"
    cex.title = 0
))

annot_scheme <- as.list(c(
    my_scheme,
    col.frame = set1[3],
    fontcolor.group = set1[3],
    fontcolor.title = set1[3],
    size = 0.2,
    col = set1[3],
    col.line = set1[3],
    fill = set1[3],
    groupAnnotation = "group",
    shape = "smallArrow"
))

# the genome
gat <- GenomeAxisTrack(showTitle = FALSE,
                       col = set1[4],
                       fontcolor = set1[4],
                       cex.title = 1,
                       col.title = set1[4],
                       name = csd_chr,
                       cex = 0.5)

# make tracks
csd_grt <- GeneRegionTrack(txdb,
                           chromosome = csd_chr,
                           start = csd_start,
                           end = csd_end,
                           name = "NCBI transcripts")
displayPars(csd_grt) <- grt_scheme

# annotation track for primers and HVR
annot <- AnnotationTrack(start = c(11771976,
                                   11771976 - 100,
                                   11772216 + 80),
                         end = c(11772216,
                                 11771976 - 80,
                                 11772216 + 100),
                         chromosome = csd_chr,
                         strand = c("+", "+", "-"),
                         group = c("HVR ", "Csd \namplicon ", "Csd \namplicon "))
displayPars(annot) <- annot_scheme

# join tracks 
grid.newpage()
extend = 0
plotTracks(list(annot, csd_grt, gat),
           sizes = c(0.2, 1, 0.1),
           cex.main = 1,
           fontface.main = 4,
           extend.right = extend,
           extend.left = extend)

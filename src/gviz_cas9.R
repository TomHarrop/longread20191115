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
csd_start <- min(start(csd_cds)) - 1e4
csd_end <- max(end(csd_cds)) + 1e3

# set up track schemes
set1 <- viridis::viridis_pal()(7)

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
    size = 0.5,
    col = set1[2],
    col.line = set1[2],
    fill = set1[2],
    fontface.group = 3,
    transcriptAnnotation = "symbol",
    cex.title = 0))

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

aln_scheme <- as.list(c(
    my_scheme,
    size = 0.5,
    fontcolor.title = set1[4],
    type = c("coverage", "pileup")))

# the genome
gat <- GenomeAxisTrack(showTitle = FALSE,
                       col = set1[4],
                       fontcolor = set1[5],
                       cex.title = 1,
                       col.title = set1[5],
                       name = csd_chr,
                       cex = 0.5)

# sequence track for reference
fasta_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
fa <- Biostrings::readDNAStringSet(fasta_file)
names(fa) <- gsub("^([^[:space:]]+).*", "\\1", names(fa)) 
st <- SequenceTrack(fa, chromosome = csd_chr)

# make tracks
csd_grt <- GeneRegionTrack(txdb,
                           chromosome = csd_chr,
                           start = csd_start,
                           end = csd_end,
                           name = "NCBI transcripts")
displayPars(csd_grt) <- grt_scheme

# alignment track (too big, have to subset)
at1 <- AlignmentsTrack(range = "data/mapped_sorted.bam",
                       referenceSequence = st,
                       isPaired = TRUE,
                       name = "All individuals",
                       chromosome = csd_chr,
                       start = csd_start,
                       end = csd_end)
displayPars(at1) <- aln_scheme

# annotation track for primers and HVR
cas9_starts <- c(11768444,
                 11780722,
                 11780351)
annot <- AnnotationTrack(start = cas9_starts,
                         end = cas9_starts + 23,
                         chromosome = csd_chr,
                         strand = c("+", "-", "-"),
                         group = c("upstream_1 ", "downstream_1 ", "downstream_2 "))
displayPars(annot) <- annot_scheme


# join tracks and highlight hypervariable region
ht1 <- HighlightTrack(trackList = list(at1, annot, csd_grt, gat),
                      start = c(11771976, cas9_starts),
                      end = c(11772216, cas9_starts + 23),
                      chromosome = csd_chr,
                      col = c(set1[7], rep(set1[6], 3)),
                      fill = c(set1[7], rep(set1[6], 3)))

# join tracks 
pdf("output/csd-cas9.pdf",
    width = 5.5,
    height = 2.95,
    pointsize = 8)

plotTracks(list(ht1),
           cex.main = 1,
           fontface.main = 4,
           extend.right = extend,
           extend.left = extend)
dev.off()
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
set1 <- viridis::viridis_pal()(6)

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
    size = 0.2,
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
    size = 0.1,
    fontcolor.title = set1[4],
    type = c("pileup"),
    showTitle = TRUE,
    cex.title = 0.5))

# the genome
gat <- GenomeAxisTrack(showTitle = FALSE,
                       col = set1[5],
                       fontcolor = set1[5],
                       cex.title = 1,
                       col.title = set1[5],
                       name = csd_chr,
                       cex = 0.5,
                       size = 0.2)

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
at1 <- AlignmentsTrack(range = "output/DR01_51_subset.bam",
                       name = "DR01 Nanopore",
                       referenceSequence = st,
                       isPaired = FALSE,
                       chromosome = csd_chr,
                       start = csd_start,
                       end = csd_end)

at2 <- AlignmentsTrack(range = "/cifs/ro_deardenlabarchive/tomharrop/projects/apis-tp-201910/output/030_process-aln/DR01_marked.bam",
                       name = "DR01 Illumina",
                       referenceSequence = st,
                       isPaired = TRUE,
                       chromosome = csd_chr,
                       start = csd_start,
                       end = csd_end)

displayPars(at1) <- aln_scheme
displayPars(at2) <- aln_scheme

# join tracks and highlight hypervariable region
ht1 <- HighlightTrack(trackList = list(at1, at2, csd_grt),
                      start = c(11771976),
                      end = c(11772216),
                      chromosome = csd_chr,
                      col = c(set1[6]),
                      fill = c(set1[6]))

# join tracks 
pdf("output/DR01.pdf",
    width = 5.5,
    height = 2.95,
    pointsize = 8)

plotTracks(list(ht1),
           cex.main = 1,
           fontface.main = 4)
dev.off()

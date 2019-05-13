# Author: Feng Li
##################### R script to draw trackes of gene model digrams ##################

library(GenomicFeatures)
library(ggbio)
# Draw gene diagram for AvrSr35 and AvrSr50 locus in Ug99 assembly
# load Ug99 annotation gff3 file
# keep only one transcript for drawing the diagram of gene models; leave the alternate transcripts out
# cat ../Puccinia_graminis_f._sp._tritici_Ug99.gff3 | grep -v 'ID=PGTUg99_.*-T2' | grep -v 'ID=PGTUg99_.*-T3' > Puccinia_graminis_f._sp._tritici_Ug99_oneTranscript.gff3
ug99_txdb = makeTxDbFromGFF("~/Downloads/Puccinia_graminis_f._sp._tritici_Ug99_oneTranscript.gff3",
                              dataSource="Ug99 gff3 from funannotate",
                              organism="Puccinia graminis")
# extract features overlapping tig2147 and tig2125 in Ug99 where they have AvrSr35 and AvrSr50
# Method ref: Bioinformatics Data Skills Book Pge 311 and https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf


head(seqlevels(ug99_txdb))

# get only the genomic features on a specific contig
seqlevels(ug99_txdb) <- "tig00002147"
seqlevels(ug99_txdb)
# now work with tig2147
# group exons by transcripts
tig2147_exons <- exonsBy(ug99_txdb, "tx")
length(tig2147_exons)
#all(unlist(seqnames(tig2147_exons)) == "tig00002147")
# give the coordinates where AvrSr35/AvrSr50 is located on this contig, and draw a plot
avr_region <- GRanges("tig00002147", IRanges(946000, 1017000))
avr_region_tx <- transcriptsByOverlaps(ug99_txdb, avr_region)
autoplot(avr_region_tx, label.color="black", color="grey43", fill="grey43")
jpeg("~/Downloads/R_to_draw_chro_gene_structure/gene_model_tig2147_geneOnly.jpg",
     width = 960*64, height = 960, units = "px", pointsize = 12,  quality = 300)
autoplot(avr_region_tx, label.color="black", color="grey43", fill="grey43")
dev.off()


# get only the genomic features on a specific contig - itg2125
seqlevels(ug99_txdb) <- "tig00002125"
seqlevels(ug99_txdb)
# now work with tig2147
# group exons by transcripts
tig2125_exons <- exonsBy(ug99_txdb, "tx")
length(tig2125_exons)
#all(unlist(seqnames(tig2147_exons)) == "tig00002147")
# give the coordinates where AvrSr35/AvrSr50 is located on this contig, and draw a plot
avr_region <- GRanges("tig00002125", IRanges(894000, 938000))
avr_region_tx <- transcriptsByOverlaps(ug99_txdb, avr_region)
autoplot(avr_region_tx, label.color="black", color="grey43", fill="grey43")
jpeg("~/Downloads/R_to_draw_chro_gene_structure/gene_model_tig2125_geneOnly.jpg",
     width = 960*64, height = 960, units = "px", pointsize = 12,  quality = 300)
autoplot(avr_region_tx, label.color="black", color="grey43", fill="grey43")
dev.off()


# Draw gene diagram for AvrSr35 and AvrSr50 locus in Pgt210 assembly

# loadannotation gff3 file; annotation pgt210_afterBP
# keep only one transcript for drawing the diagram of gene models; leave the alternate transcripts out
# cat ../Puccinia_graminis_tritici_21-0.gff3 | grep -v 'ID=PGT21_.*-T2' | grep -v 'ID=PGT21_.*-T3' > Puccinia_graminis_tritici_21-0_oneTranscript.gff3

pgt210_txdb = makeTxDbFromGFF("~/Downloads/Puccinia_graminis_tritici_21-0_oneTranscript.gff3",
                            dataSource="21-0 gff3 from funannotate",
                            organism="Puccinia graminis")
head(seqlevels(pgt210_txdb))

# get only the genomic features on a specific contig
seqlevels(pgt210_txdb) <- "tig00001259_02"
seqlevels(pgt210_txdb)
# now work with tig1259_02
# group exons by transcripts
tig1259_02_exons <- exonsBy(pgt210_txdb, "tx")
length(tig1259_02_exons)

# give the coordinates where AvrSr35/AvrSr50 is located on this contig, and draw a plot
avr_region <- GRanges("tig00001259_02", IRanges(214000, 264000))
avr_region_tx <- transcriptsByOverlaps(pgt210_txdb, avr_region)
autoplot(avr_region_tx, label.color="black", color="grey43", fill="grey43")
jpeg("~/Downloads/R_to_draw_chro_gene_structure/gene_model_tig1259_02_geneOnly.jpg",
     width = 960*64, height = 960, units = "px", pointsize = 12,  quality = 300)
autoplot(avr_region_tx, label.color="black", color="grey43", fill="grey43")
dev.off()


# get only the genomic features on a specific contig - tig150
seqlevels(pgt210_txdb) <- "tig00000150"
seqlevels(pgt210_txdb)
# now work with tig150
# group exons by transcripts
tig150_exons <- exonsBy(pgt210_txdb, "tx")
length(tig150_exons)

# give the coordinates where AvrSr35/AvrSr50 is located on this contig, and draw a plot
avr_region <- GRanges("tig00000150", IRanges(56000, 70000))
avr_region_tx <- transcriptsByOverlaps(pgt210_txdb, avr_region)
autoplot(avr_region_tx, label.color="black", color="grey43", fill="grey43")
jpeg("~/Downloads/R_to_draw_chro_gene_structure/gene_model_tig150_geneOnly.jpg",
     width = 960*64, height = 960, units = "px", pointsize = 12,  quality = 300)
autoplot(avr_region_tx, label.color="black", color="grey43", fill="grey43")
dev.off()


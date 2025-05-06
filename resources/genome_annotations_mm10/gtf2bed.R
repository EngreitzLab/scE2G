setwd("/maps/projects/ralab_nnfc-AUDIT/people/lpm537/software/scE2G_pipeline/250420_mm10/scE2G/resources/genome_annotations_mm10")

library(rtracklayer)
library(dplyr)
library(GenomeInfoDb)
library(GenomicRanges)


gtf <- import("/maps/projects/ralab_nnfc-AUDIT/people/lpm537/project/Sperm/processed/Masahiro/scmulti/genes.gtf.gz", format = "gtf")

# calculate transcript lenths
gtf.exons <- gtf[gtf$type == "exon"]
transcript_ids <- gtf.exons$transcript_id
transcript_lengths <- gtf.exons %>%
  split(., .$transcript_id) %>%
  lapply(reduce) %>%            
  sapply(function(x) sum(width(x))) 

gtf.transcript <- gtf[gtf$type == "transcript"]
gtf.transcript$length = transcript_lengths[gtf.transcript$transcript_id]

gtf.tss <- gtf.transcript
start(gtf.tss)[as.character(strand(gtf.tss)) == "-"] = end(gtf.tss)[as.character(strand(gtf.tss)) == "-"]
end(gtf.tss)[as.character(strand(gtf.tss)) == "+"] = start(gtf.tss)[as.character(strand(gtf.tss)) == "+"]

gtf.tss$gene_tss = paste(gtf.tss$gene_id,
                         seqnames(gtf.tss),
                         start(gtf.tss),
                         sep = "_")
gene_tss.count = table(gtf.tss$gene_tss)
gtf.tss$gene_tss.count = gene_tss.count[gtf.tss$gene_tss]


gtf.tss.sort = gtf.tss[order(gtf.tss$length,decreasing = T)]
gtf.tss.sort = gtf.tss.sort[order(gtf.tss.sort$gene_tss.count,decreasing = T)]
gtf.tss.sort = gtf.tss.sort[order(gtf.tss.sort$gene_id,decreasing = T)]
gtf.tss.sort

gtf.tss.uniq = gtf.tss.sort[!duplicated(gtf.tss.sort$gene_id)]
gtf.tss.uniq

saveRDS(gtf.tss.uniq,
        "gtf.tss.uniq.rds")

df.gencode.vM25.protein_coding = read.delim("../../ENCODE_rE2G/ABC/reference/mm10_vM25/gencode.vM25.protein_coding.genes.bed")
gene_name.dup  <- gtf.tss.uniq$gene_name[duplicated(gtf.tss.uniq$gene_name)]
gene_id.dup <- gtf.tss.uniq[gtf.tss.uniq$gene_name %in% gene.dup]$gene_id
gene_id.rm <- gene_id.dup[!gene_id.dup %in% df.gencode.vM25.protein_coding$Ensembl_ID]

gtf.tss.uniq.2 = gtf.tss.uniq[!gtf.tss.uniq$gene_id %in% gene_id.rm]

table(duplicated(gtf.tss.uniq.2$gene_name))
gene_name.dup2  <- gtf.tss.uniq.2$gene_name[duplicated(gtf.tss.uniq.2$gene_name)]
gtf.tss.uniq.2[gtf.tss.uniq.2$gene_name %in% gene_name.dup2]
df.gencode.vM25.protein_coding[df.gencode.vM25.protein_coding$name %in% gene_name.dup2,]

gtf.tss.uniq.3 = gtf.tss.uniq.2[!gtf.tss.uniq.2$gene_id %in% c("ENSMUSG00000118554",
                                                               "ENSMUSG00000118409")]
table(duplicated(gtf.tss.uniq.3$gene_name))



gtf.tss.uniq.filter = sort(gtf.tss.uniq.3)
gtf.tss.uniq.filter = gtf.tss.uniq.filter[gtf.tss.uniq.filter$gene_type == "protein_coding"]
gtf.tss.uniq.filter = 
  gtf.tss.uniq.filter[seqnames(gtf.tss.uniq.filter) %in% paste("chr",c(1:19,"X","Y"),sep="")]
seqlevels(gtf.tss.uniq.filter) = paste("chr",c(1:19,"X","Y"),sep="")

gtf.tss500.uniq.filter = gtf.tss.uniq.filter
start(gtf.tss500.uniq.filter) = pmax(start(gtf.tss500.uniq.filter) - 250, 1)
end(gtf.tss500.uniq.filter) = end(gtf.tss500.uniq.filter) + 250

chrom_info <- getChromInfoFromUCSC("mm10")
chrom_sizes <- setNames(chrom_info$size, chrom_info$chrom)
seqlengths(gtf.tss500.uniq.filter) <- chrom_sizes[names(seqlengths(gtf.tss500.uniq.filter))]
gtf.tss500.uniq.filter <- trim(gtf.tss500.uniq.filter)


df.out = data.frame("chr" = seqnames(gtf.tss500.uniq.filter),
                    "start" = start(gtf.tss500.uniq.filter),
                    "end" = end(gtf.tss500.uniq.filter),
                    "name" = gtf.tss500.uniq.filter$gene_name,
                    "score" = 0,
                    "strand" = strand(gtf.tss500.uniq.filter),
                    "Ensembl_ID" = gtf.tss500.uniq.filter$gene_id,
                    "gene_type" = gtf.tss500.uniq.filter$gene_type)

cat("#",file = "CollapsedGeneBounds.mm10.gencode.M23.TSS500bp.bed")
write.table(df.out,
            "CollapsedGeneBounds.mm10.gencode.M23.TSS500bp.bed",
            sep = "\t",
            quote = F,
            row.names = F,
            append = T)

names(gtf.transcript) = gtf.transcript$transcript_id
gtf.transcript.output = gtf.transcript[gtf.tss500.uniq.filter$transcript_id]
df.out = data.frame("chr" = seqnames(gtf.transcript.output),
                    "start" = start(gtf.transcript.output),
                    "end" = end(gtf.transcript.output),
                    "name" = gtf.transcript.output$gene_name,
                    "score" = 0,
                    "strand" = strand(gtf.transcript.output),
                    "Ensembl_ID" = gtf.transcript.output$gene_id,
                    "gene_type" = gtf.transcript.output$gene_type)

cat("#",file = "CollapsedGeneBounds.mm10.gencode.M23.bed")
write.table(df.out,
            "CollapsedGeneBounds.mm10.gencode.M23.bed",
            sep = "\t",
            quote = F,
            row.names = F,
            append = T)

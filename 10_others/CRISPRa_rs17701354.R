MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
MPRA_DIR <- "/Volumes/evo_immuno_pop/MPRA"
FIGURE_DIR <- sprintf("%s/figures/MPRA_count_exp6_analysis3/", MPRA_DIR)

library(data.table)
library(ComplexHeatmap)
    library(ggplot2)
    CHR <- 'chr3'

    mySTART <- 45.0e6
    myEND <- 46.7e6
    mySTART <- 45.7e6
    myEND <- 46.7e6
    Gene_annotation_DT <- fread(sprintf("%s/data/Gencode_annotated.txt", MPRA_DIR))
    Gene_annotation_DT[, gene_id := gsub("\\..*", "", gene_id)]
    DEgenes <- fread(sprintf("%s/data/03_CRISPRa/B7075-10_Li_2025-12-15_v1/rnaseq_v1/post_analysis/rnadiff_v2/rnadiff.csv", MPRA_DIR), skip = 1)
    DEgenes[is.na(log2FoldChangeNotShrinked),log2FoldChangeNotShrinked:=0]
    DEgenes[is.na(lfcSENotShrinked),lfcSENotShrinked:=0]
    DEgenes[is.na(baseMean),baseMean:=0]
    DEgenes[is.na(stat),stat:=0]
    DEgenes[is.na(pvalue),pvalue:=1]
    DEgenes[is.na(padj),padj:=1]

    DEgenes <- merge(DEgenes, Gene_annotation_DT[,.(gene_id,chrom=seqid,gene_name,strand,start,end)], by='gene_id',suffix=c('','_new'))


    locus_genes <- DEgenes[chrom==CHR & start >= mySTART & end <= myEND,.(gene_id,gene_name=gene_name_new,strand,start,end,log2FoldChangeNotShrinked,lfcSENotShrinked,pvalue,padj)]

Table <- fread(sprintf("%s/data/03_CRISPRa/B7075-10_Li_2025-12-15_v1/rnaseq_v1/post_analysis/rnadiff_v2/counts/counts_vst_norm.csv", MPRA_DIR))


Table <- fread(sprintf("%s/data/03_CRISPRa/B7075-10_Li_2025-12-15_v1/rnaseq_v1/post_analysis/rnadiff_v2/counts/counts_normed.csv", MPRA_DIR))
Table <- melt(Table)
Table[, raw_count := value]
Table[, value := log2(value + 1)]
Table[, variable := gsub("rs17713054_sgRNA", "rs17713054_sgRNA_rep", variable)]
Table[, rep := gsub("CRISPRa_(Neg|Pos)_(GAL4|rs17713054)_sgRNA_rep([1-8])_S([0-9]+)", "rep\\3", variable)]
Table[, status := gsub("CRISPRa_(Neg|Pos)_(GAL4|rs17713054)_sgRNA_rep([1-8])_S([0-9]+)", "\\1_\\2", variable)]
Table[, value_scaled := scale(value), by = V1]
TableExpr <- dcast(Table, V1 ~ status + rep, value.var = "value_scaled")

Table <- dcast(Table, V1 + rep ~ status, value.var = "value")
Table[, FC := Pos_rs17713054 - Neg_GAL4]

Table[V1 %in% DEgenes[padj < 0.05, gene_id], ]

Table_locus <- merge(expand.grid(V1=locus_genes$gene_id,rep=c('rep1','rep2','rep3','rep4','rep5','rep6')),Table[V1%in%locus_genes$gene_id,],all.x=TRUE)
Table_locus <- as.data.table(Table_locus)
Table_locus[is.na(Pos_rs17713054),Pos_rs17713054:=0]
Table_locus[is.na(Neg_GAL4),Neg_GAL4:=0]
Table_locus[is.na(FC),FC:=0]
Table_locus <- merge(Table_locus, locus_genes[,.(gene_id, gene_name)], by.x='V1', by.y='gene_id')

# color_celltype  <- c("HepG2" = "#AE89C4", "A549-ACE2" = "#81C16E", "A549" = "#62cf2f", "K562" = "#BF9454")
p <- ggplot(locus_genes, aes(x=factor(gene_name,locus_genes[order(start),gene_name]), y=log2FoldChangeNotShrinked)) 
p <- p + geom_bar(stat='identity', fill=color_celltype['A549-ACE2']) + ylab('log2 Fold Change')
p <- p + geom_errorbar(aes(ymin=log2FoldChangeNotShrinked-1.96*lfcSENotShrinked, ymax=log2FoldChangeNotShrinked+1.96*lfcSENotShrinked), width=.3) 
p <- p + theme_plot(rotate.x=90, lpos='none', fontsize=11) + xlab('locus genes') 
p <- p + geom_point(data=Table_locus,aes(x=factor(gene_name,locus_genes[order(start),gene_name]), y=FC), alpha=.5,size=.4)

dir.create(sprintf('%s/07_cripsra/', FIGURE_DIR))
pdf(sprintf('%s/07_cripsra/barplot_locus_genes.pdf', FIGURE_DIR), height = 3, width = 4)
print(p)
dev.off()

DE_heatmap_data <- merge(DEgenes[, .(gene_id, gene_name)], dcast(Table[V1 %in% DEgenes[padj < 0.05, gene_id], ], V1 ~ rep, value.var = "FC"), by.x = "gene_id", by.y = "V1")
DE_heatmap_mat <- as.matrix(DE_heatmap_data[, .(rep1, rep2, rep3, rep4, rep5, rep6)])
pdf(sprintf("%s/data/03_CRISPRa/B7075-10_Li_2025-12-15_v1/rnaseq_v1/post_analysis/rnadiff_v2/DE_heatmap.pdf", MPRA_DIR), height = 10, width = 5)
pheatmap(DE_heatmap_mat)
dev.off()

DE_heatmap_data <- merge(DEgenes[, .(gene_id, gene_name)], TableExpr[V1 %in% DEgenes[padj < 0.05, gene_id], ], by.x = "gene_id", by.y = "V1")
DE_heatmap_mat <- as.matrix(DE_heatmap_data[, -c("gene_id", "gene_name")])
row.names(DE_heatmap_mat) <- DE_heatmap_data$gene_name

pdf(sprintf("%s/data/03_CRISPRa/B7075-10_Li_2025-12-15_v1/rnaseq_v1/post_analysis/rnadiff_v2/RNA_heatmap.pdf", MPRA_DIR), height = 10, width = 5)
pheatmap(DE_heatmap_mat, cluster_cols = FALSE)
dev.off()

random_genes <- sample(DEgenes[, gene_id], 50)
DE_heatmap_data <- merge(DEgenes[, .(gene_id, gene_name)], TableExpr[V1 %in% random_genes, ], by.x = "gene_id", by.y = "V1")
DE_heatmap_mat <- as.matrix(DE_heatmap_data[, -c("gene_id", "gene_name")])
row.names(DE_heatmap_mat) <- DE_heatmap_data$gene_name

pdf(sprintf("/pasteur/helix/projects/IGSR/RNA_heatmap_random.pdf", MPRA_DIR), height = 10, width = 5)
pheatmap(DE_heatmap_mat, cluster_cols = FALSE)
dev.off()

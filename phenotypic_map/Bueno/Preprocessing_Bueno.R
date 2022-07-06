##################################################################
############# Pre-processing Bueno cohort ########################
##################################################################

require(data.table)
require(openxlsx)
require(DESeq2)
require(rtracklayer)

## Load annotation
# Gene annotation GENCODE v37
ref.37 = rtracklayer::import("/path/to/gencode/annotation/file/gencode.v33lift37.annotation.gtf.gz", format="gtf")
ref.37 = as.data.frame(ref.37)
ref.37.gene = ref.37[which(ref.37$type == "gene"),]
ref.37.gene$seqnames = gsub("^chr", "", ref.37.gene$seqnames)
ref.37.gene = as.data.table(ref.37.gene)
setnames(ref.37.gene,c("gene_name","gene_id"),c("geneSymbol","id"))

# Samples annotation
samp_info = read.xlsx("/path/to/supplementary/table1/41588_2016_BFng3520_MOESM229_ESM.xlsx", sheet = 2 ,startRow = 2)
samp_info$Tumor.ID[which(!startsWith(samp_info$Tumor.ID,"M"))] = paste0("M",samp_info$Tumor.ID[which(!startsWith(samp_info$Tumor.ID,"M"))])
samp_info$`Histology-reduced.(WHO.categories.based.on.diagnosis.reported.in.the.surgical.pathology.report)`[which(samp_info$Tumor.ID == "M78PT")] = "MME"
samp_info$`%.Sarcomatoid`[which(samp_info$Tumor.ID == "M78PT")] = 0

## Pre-processing of RNA-seq data
# Load raw read counts
expr_coun = read.csv("/path/to/gene_count_matrix_1pass.csv",row.names = 1)
colnames(expr_coun)[which(startsWith(colnames(expr_coun),"X"))] = sapply(strsplit(colnames(expr_coun)[which(startsWith(colnames(expr_coun),"X"))],"X"),"[[",2)
colnames(expr_coun)[which(!startsWith(colnames(expr_coun),"M"))] = paste0("M",colnames(expr_coun)[which(!startsWith(colnames(expr_coun),"M"))])
expr_coun = expr_coun[,which(colnames(expr_coun) %in% samp_info$Tumor.ID[which(samp_info$`preop.Treatment.(1=.treated;.0=.untreated)` == 0 & samp_info$`RNA-seq` == "Y")])]

## Normalisation
deseqexpr = DESeqDataSetFromMatrix(expr_coun, colData = data.frame(colnames(expr_coun)), design = ~ 1, tidy = F)
vstexpr = varianceStabilizingTransformation(deseqexpr, blind = T)

expr_annot = read.table("/path/to/602PT_pass1_gene_abund.tab",header = T,sep = "\t")[,1:6]
annot_ordered = expr_annot[sapply(rownames(vstexpr), function(x) which(expr_annot$Gene.ID==x)[1]),]
vstexpr_nosex = vstexpr[!annot_ordered$Reference %in% c("chrM","chrX","chrY"),]
vstexpr_nosex = assay(vstexpr_nosex)

# Minimal level filter based on FPKM
expr_fpkm = read.csv("/path/to/gene_FPKM_matrix.csv",row.names = 1)
colnames(expr_fpkm) = sapply(strsplit(colnames(expr_fpkm),"FPKM."),"[[",2)
colnames(expr_fpkm)[which(!startsWith(colnames(expr_fpkm),"M"))] = paste0("M",colnames(expr_fpkm)[which(!startsWith(colnames(expr_fpkm),"M"))])
expr_fpkm = expr_fpkm[,which(colnames(expr_fpkm) %in% samp_info$Tumor.ID[which(samp_info$`preop.Treatment.(1=.treated;.0=.untreated)` == 0 & samp_info$`RNA-seq` == "Y")])]

FPKM.diff = apply(expr_fpkm, 1, max, na.rm=TRUE) - apply(expr_fpkm, 1, min, na.rm=TRUE)
names(FPKM.diff) = rownames(expr_fpkm)
FPKM.diff = FPKM.diff[which(FPKM.diff >= 1)]

vstexpr_nosex_red = vstexpr_nosex[which(rownames(vstexpr_nosex) %in% names(FPKM.diff)),]

# Select the 5,000 most variable genes
vv = apply(vstexpr_nosex_red,1,var)
cv = cumsum(sort(vv,decreasing = T))/sum(vv)
vstexpr_nosex_red = vstexpr_nosex_red[order(match(rownames(vstexpr_nosex_red),names(cv))),]
D_expr_red = vstexpr_nosex_red[1:5000,]

## Pre-processing of DNA gene alterations data
som_mut = read.table(file = "/path/to/som_mut_annot.txt", header = T, sep="\t") # Bueno Supplementary Table 5 annotated with annovar

# Remove non-chimionaif samples
som_mut = som_mut[which(som_mut$Tumor.ID %in% samp_info$Tumor.ID[which(samp_info$`preop.Treatment.(1=.treated;.0=.untreated)` == 0 & samp_info$Exome == "Y")]),]

# Keep only the coding mutations
som_mut = som_mut[which(som_mut$Func.wgEncodeGencodeBasicV33lift37 %in% c("exonic","exonic;splicing","splicing","ncRNA_exonic;splicing")),]

# Remove sex-chromosomes
som_mut = som_mut[which(!som_mut$Chr %in% c("chrX","chrY","chrM")),]

# Cases where two genes are annotated in Gene.ensGene
colnames(som_mut)[6:10] = c("Func.ensGene","Gene.ensGene","GeneDetail.ensGene","ExonicFunc.ensGene","AAChange.ensGene")

#Need to use ENSEMBL IDs to make identify the recurrently altered genes
som_mut$ENS_ID = sapply(1:nrow(som_mut), function(i) if(grepl("ENST",som_mut$AAChange.ensGene[i])){strsplit(som_mut$AAChange.ensGene[i],":")[[1]][2]}else{strsplit(som_mut$GeneDetail.ensGene[i],":")[[1]][1]})
som_mut$ENS_ID[which(is.na(som_mut$ENS_ID))] = som_mut$Gene.ensGene[which(is.na(som_mut$ENS_ID))]

som_mut$ENS_ID[which(!startsWith(som_mut$ENS_ID,"ENST"))] = sapply(som_mut$ENS_ID[which(!startsWith(som_mut$ENS_ID,"ENST"))], function(x) unique(ref.37.gene$id[which(ref.37.gene$geneSymbol == x)]))
som_mut$Gene.id1 = sapply(1:nrow(som_mut), function(i) if(startsWith(som_mut$ENS_ID[i],"ENST")){unique(ref.37$gene_id[which(ref.37$transcript_id == som_mut$ENS_ID[i])])}else{som_mut$ENS_ID[i]})
som_mut$Gene1 = sapply(som_mut$Gene.id1, function(x) ref.37.gene$geneSymbol[which(ref.37.gene$id == x)])
som_mut$Gene2 = sapply(1:nrow(som_mut), function(i) if(grepl(";",som_mut$Gene.ensGene[i]) & som_mut$Func.ensGene[i] %in% c("exonic","exonic;splicing","splicing")){ifelse(length(unique(unlist(strsplit(som_mut$Gene.ensGene[i],";")))) > 1, unlist(strsplit(som_mut$Gene.ensGene[i],";"))[which(!unlist(strsplit(som_mut$Gene.ensGene[i],";")) == som_mut$Gene1[i])],NA)}else{NA})
som_mut$Gene.id2 = sapply(som_mut$Gene2, function(x) if(!is.na(x)){ref.37.gene$id[which(ref.37.gene$geneSymbol == x)][1]}else{NA})
som_mut$Gene.id2[which(som_mut$Gene2 == "ECE2")] = "ENSG00000145194.18_3"

# Keep only the expressed genes
FPKM.max = apply(expr_fpkm,1,max,na.rm=TRUE)
names(FPKM.max) = rownames(expr_fpkm)
FPKM.max = FPKM.max[which(FPKM.max >= 0.01)]

som_mut$Gene1[which(!sub('\\_.*', '', som_mut$Gene.id1) %in% names(FPKM.max))] = NA
som_mut$Gene.id1[which(!sub('\\_.*', '', som_mut$Gene.id1) %in% names(FPKM.max))] = NA

som_mut$Gene2[which(!is.na(som_mut$Gene.id2) & !sub('\\_.*', '', som_mut$Gene.id2)  %in% names(FPKM.max))] = NA
som_mut$Gene.id2[which(!is.na(som_mut$Gene.id2) & !sub('\\_.*', '', som_mut$Gene.id2)  %in% names(FPKM.max))] = NA
som_mut = som_mut[-which(is.na(som_mut$Gene1) & is.na(som_mut$Gene2)),]

gene_switch = which(is.na(som_mut$Gene1) & !is.na(som_mut$Gene2))
som_mut$Gene.id1[gene_switch] = som_mut$Gene.id2[gene_switch]
som_mut$Gene1[gene_switch] = som_mut$Gene2[gene_switch]
som_mut$Gene.id2[gene_switch] = NA
som_mut$Gene2[gene_switch] = NA

# Keep only protein-coding and lncRNA genes
som_mut = som_mut[which(som_mut$Gene.id1 %in% ref.37.gene$id[which(ref.37.gene$gene_type %in% c("protein_coding","lncRNA"))] | !is.na(som_mut$Gene.id2)),]

# Conversion from annovar to maftools
som_mut$variant_id = paste0(som_mut$Tumor.ID,"_",som_mut$Chr,":",som_mut$Start,"-",som_mut$End)
colnames(som_mut)[6:10] = c("Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene")
write.table(som_mut, file = "/intermadiary/path/som_mut.txt", row.names=FALSE,col.names = T, sep="\t")

var.annovar.maf = annovarToMaf(annovar = "/intermadiary/path/som_mut.txt", Center = 'Bueno', refBuild = 'hg19',
                               tsbCol = 'Tumor.ID', table = 'ensGene', basename = "/intermadiary/path/var_annovar_maf", ens2hugo = T)

som_mut$Variant_Classification = sapply(1:nrow(som_mut), function(i) var.annovar.maf$Variant_Classification[which(var.annovar.maf$variant_id == som_mut$variant_id[i])])
som_mut = som_mut[-which(som_mut$Variant_Classification == "Missense_Mutation" & (is.na(som_mut$REVEL) | as.numeric(som_mut$REVEL) < 0.5)),] 

# Building DNA Alt data set
samples_mut = sort(unique(som_mut$Tumor.ID))
genes = unique(som_mut$Gene.id1)
D_mut = matrix(0, length(genes), length(samples_mut), dimnames = list(genes, samples_mut))
for(i in 1:length(genes)){
  D_mut[i,] = sapply(colnames(D_mut), function(x) if(x %in% som_mut$Tumor.ID[which(som_mut$Gene.id1 == genes[i])]){1}else{0})
}
D_alt = D_mut[which(rowSums(D_mut) > 2),]

## Homogeneisation
path.MOFA = "/path/to/preprocessed/matrices/"
All.samples = sort(unique(c(colnames(D_expr_red), colnames(D_alt))))

D_alt_MOFA = matrix(NA,nrow(D_alt),length(All.samples[which(!All.samples%in%colnames(D_alt))]),dimnames = list(rownames(D_alt),All.samples[which(!All.samples%in%colnames(D_alt))]))
D_alt_MOFA = as.data.frame(D_alt_MOFA)
D_alt_MOFA = cbind(D_alt,D_alt_MOFA)
D_alt_MOFA = D_alt_MOFA[,order(match(colnames(D_alt_MOFA),All.samples))]
D_alt_MOFA = as.matrix(D_alt_MOFA)
save(D_alt_MOFA, paste0(path.MOFA,"D_alt_MOFA.RData"))

D_expr_MOFA = matrix(NA,nrow(D_expr_red),length(All.samples[which(!All.samples%in%colnames(D_expr_red))]),dimnames = list(rownames(D_expr_red ),All.samples[which(!All.samples%in%colnames(D_expr_red))]))
D_expr_MOFA = as.data.frame(D_expr_MOFA)
D_expr_MOFA = cbind(D_expr_red ,D_expr_MOFA)
D_expr_MOFA = D_expr_MOFA[,order(match(colnames(D_expr_MOFA),All.samples))]
D_expr_MOFA = as.matrix(D_expr_MOFA)
save(D_expr_MOFA, paste0(path.MOFA,"D_expr_MOFA.RData"))

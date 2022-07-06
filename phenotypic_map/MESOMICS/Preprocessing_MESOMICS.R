##################################################################
############# Pre-processing MESOMICS cohort #####################
##################################################################

require(data.table)
require(openxlsx)
require(DESeq2)
require(rtracklayer)

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
sexcpgs = rownames(Locations)[Locations$chr%in%c("chrX","chrY","chrM")]

devtools::install_github("walaj/roverlaps")
library(roverlaps)

## Load annotation
# Gene annotation GENCODE v33
ref.33 = rtracklayer::import("/path/to/gencode/annotation/file/gencode.v33.annotation.gtf", format="gtf")
ref.33 = as.data.frame(ref.33)
ref.33.gene = ref.33[which(ref.33$type == "gene"),]
ref.33.gene$seqnames = gsub("^chr", "", ref.33.gene$seqnames)
ref.33.gene = as.data.table(ref.33.gene)
setnames(ref.33.gene,c("gene_name","gene_id"),c("geneSymbol","id"))

# Samples annotation
Data.Clin = read.xlsx("/path/to/supplementary/table1/SupplementaryTable1-SamplesOverview.xlsx", sheet = 1 ,startRow = 2)

# CpGs annotation
Annotation = read.delim("/path/to/EPIC850K/array/manifest/b5/Infinium_MethylationEPIC_v1b5_manifest.txt")

## Pre-processing of RNA-seq data
# Load raw read counts
expr_all = read.csv("/path/to/gene_count_matrix_1pass.csv",row.names = 1)
expr_all = expr_all[,which(colnames(expr_all) %in% Data.Clin$Cohort)]
expr_all = expr_all[,-which(colnames(expr_all) == "MESO_054_T")]

## Normalisation
deseqexpr = DESeqDataSetFromMatrix(expr_all, colData = data.frame(colnames(expr_all)), design = ~ 1, tidy = F)
vstexpr = varianceStabilizingTransformation(deseqexpr, blind = T)

expr_annot = read.table("/path/to/MESO_001_T_pass1_gene_abund.tab",header = T,sep = "\t")[,1:6]
annot_ordered = expr_annot[sapply(rownames(vstexpr), function(x) which(expr_annot$Gene.ID==x)[1]),]
vstexpr_nosex = vstexpr[which(!annot_ordered$Reference %in% c("chrM","chrX","chrY")),]
vstexpr_nosex = assay(vstexpr_nosex)

# Minimal level filter based on FPKM
Gene_exp = read.csv("/path/to/gene_FPKM_matrix.csv",row.names = 1)
Gene_exp = Gene_exp[,which(startsWith(colnames(Gene_exp),"FPKM"))]
colnames(Gene_exp) = sapply(strsplit(colnames(Gene_exp),"FPKM."),"[[",2)
Gene_exp = Gene_exp[,which(colnames(Gene_exp) %in% colnames(expr_all))]

FPKM.diff = apply(Gene_exp,1,max,na.rm=TRUE) - apply(Gene_exp,1,min,na.rm=TRUE)
names(FPKM.diff) = rownames(Gene_exp)
FPKM.diff = FPKM.diff[which(FPKM.diff >= 1)]

vstexpr_nosex_red = vstexpr_nosex[which(rownames(vstexpr_nosex) %in% names(FPKM.diff)),]

# Select the 5,000 most variable genes
vv = apply(vstexpr_nosex_red,1,var)
cv = cumsum(sort(vv,decreasing = T))/sum(vv)
vstexpr_nosex_red = vstexpr_nosex_red[order(match(rownames(vstexpr_nosex_red),names(cv))),]
D_expr_redB = vstexpr_nosex_red[1:5000,]

## Pre-processing of DNA methylation array data
# Load M-values matrix
D_met = read.table("/path/to/NormFilteredMTable_noInf_139.csv",header = TRUE, sep = ",")
rownames(D_met) = D_met$X
D_met = D_met[,-1]
D_met = D_met[,which(colnames(D_met) %in% Data.Clin$Cohort)]
any(rownames(D_met) %in% sexcpgs)

# Load beta-values matrix
D_beta = read.table("/path/to/NormalisedFilteredBetaTable_139.csv",header = TRUE, sep = ",")
rownames(D_beta) = D_beta$X
D_beta = D_beta[,-1]
D_beta = D_beta[,which(colnames(D_beta) %in% Data.Clin$Cohort)]
any(rownames(D_beta) %in% sexcpgs)

# Minimal level filter based on beta-values
Beta.diff = apply(D_beta,1,max,na.rm=TRUE) - apply(D_beta,1,min,na.rm=TRUE)
names(Beta.diff) = rownames(D_beta)
Beta.diff = Beta.diff[which(Beta.diff >= 0.1)]
D_met_red = D_met[which(rownames(D_met) %in% names(Beta.diff)),]

# Select the 5,000 most variable Promoter CpGs
D_met.pro = D_met_red[which(rownames(D_met_red) %in% as.character(Annotation$IlmnID[which(Annotation$class == "Promoter")])),]
vv = apply(D_met.pro,1,var)
cv = cumsum(sort(vv,decreasing = T))/sum(vv)
D_met.pro = D_met.pro[order(match(rownames(D_met.pro),names(cv))),]
D_met.proB = D_met.pro[1:5000,]

# Select the 5,000 most variable Gene Body CpGs
D_met.bod = D_met_red[which(rownames(D_met_red) %in% as.character(Annotation$IlmnID[which(Annotation$class == "Body")])),]
vv = apply(D_met.bod,1,var)
cv = cumsum(sort(vv,decreasing = T))/sum(vv)
D_met.bod = D_met.bod[order(match(rownames(D_met.bod),names(cv))),]
D_met.bodB = D_met.bod[1:5000,]

# Select the 5,000 most variable Enhancer CpGs
D_met.enh = D_met_red[which(rownames(D_met_red) %in% as.character(Annotation$IlmnID[which(Annotation$class == "Enhancer")])),]
vv = apply(D_met.enh,1,var)
cv = cumsum(sort(vv,decreasing = T))/sum(vv)
D_met.enh = D_met.enh[order(match(rownames(D_met.enh),names(cv))),]
D_met.enhB = D_met.enh[1:5000,]

## Pre-processing of CN data
dataset.seg = read.table("/path/to/supplementary/table9/mesomics.discovery_cohort.purple.cnv.somatic.filtered_27042022.txt", header = T, sep = "\t", stringsAsFactors = F) #Table S9
colnames(dataset.seg)[which(colnames(dataset.seg) == "sample")] = "ID_IARC"

dataset.pur = data.frame(ID_IARC = Data.Clin$Cohort[which(!is.na(Data.Clin$WGS))],
                         Purity = Data.Clin$Purity.WGS.PURPLE[which(!is.na(Data.Clin$WGS))], stringsAsFactors = F)

dataset.seg$chromosome = gsub("^chr", "", dataset.seg$chromosome)
dataset.seg = dataset.seg[which(!dataset.seg$chromosome %in% c("X","Y")),]
dataset.seg$segmentID = paste0("chr",dataset.seg$chromosome,":",dataset.seg$start,"-",dataset.seg$end,"_",sapply(strsplit(dataset.seg$ID_IARC,"_"),"[[",2))
colnames(dataset.seg)[which(colnames(dataset.seg) == "id")] = "ID"

# Gene annotation
dataset.seg = as.data.table(dataset.seg)
setnames(dataset.seg,"chromosome","seqnames")

rm(dataset.seg_tmp)
for(i in 1:nrow(dataset.seg)){ 
  rm(d1,ro)
  if(!exists("dataset.seg_tmp")){
    d1 = dataset.seg[i,]
    ro = roverlaps::roverlaps(d1, ref.33.gene)
    ro[, geneSymbol := ref.33.gene$geneSymbol[subject.id]]
    d1[, geneSymbol := paste(unique(ro$geneSymbol), collapse=",")]
    ro[, id := ref.33.gene$id[subject.id]]
    d1[, id := paste(unique(ro$id), collapse=",")]
    dataset.seg_tmp = d1
  }else{
    d1 = dataset.seg[i,]
    ro = roverlaps::roverlaps(d1, ref.33.gene)
    ro[, geneSymbol := ref.33.gene$geneSymbol[subject.id]]
    d1[, geneSymbol := paste(unique(ro$geneSymbol), collapse=",")]
    ro[, id := ref.33.gene$id[subject.id]]
    d1[, id := paste(unique(ro$id), collapse=",")]
    temp_data = d1
    
    dataset.seg_tmp = rbind(dataset.seg_tmp, temp_data)
    rm(temp_data)
  }
}
dataset.seg = dataset.seg_tmp
dataset.seg$purity = sapply(dataset.seg$ID_IARC, function(x) unique(dataset.pur$Purity[which(dataset.pur$ID_IARC == x)]))
dataset.seg = as.data.frame(dataset.seg)

# Build Total dataset
samples_CNV = sort(unique(dataset.seg$ID_IARC))
D_cnv = c()
chrs = 1:22
rm(chr,g,i,x,ev.chr,o,ev)
for(chr in chrs){
  dataset.seg_tmp = dataset.seg[dataset.seg$seqnames == chr,]
  dataset.seg_tmp = as.data.frame(dataset.seg_tmp %>% group_by(start) %>% arrange(ID_IARC))
  
  genes = sort(unique(unlist(strsplit(dataset.seg_tmp$id,","))))
  
  events = sapply(genes, function(g) paste(which(dataset.seg_tmp$id == g | startsWith(dataset.seg_tmp$id, paste0(g,",")) | 
                                                   grepl(paste0(",",g,","),dataset.seg_tmp$id) | endsWith(dataset.seg_tmp$id,paste0(",",g))), collapse = ","))
  events.com = sapply(names(events), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events[g],","))), function(i) if(ref.33.gene$start[which(ref.33.gene$id == g)]>=dataset.seg_tmp$start[i] & ref.33.gene$end[which(ref.33.gene$id == g)]<=dataset.seg_tmp$end[i]){i})), collapse = ",") ) # Segments where the gene is complete
  events.cut = sapply(names(events), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events[g],","))), function(i) if(!ref.33.gene$start[which(ref.33.gene$id == g)]>=dataset.seg_tmp$start[i] | !ref.33.gene$end[which(ref.33.gene$id == g)]<=dataset.seg_tmp$end[i]){i})), collapse = ",") ) # Segments where the gene is cut
  
  # For segments where the gene is cut
  events.cut = sapply(names(events.cut), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events.cut[g],","))), function(i) if(ref.33.gene$end[which(ref.33.gene$id == g)]>dataset.seg_tmp$end[i] & (!(i+1) %in% as.numeric(unlist(strsplit(events.cut[g],","))) | ((i+1) %in% as.numeric(unlist(strsplit(events.cut[g],","))) & dataset.seg_tmp$ID_IARC[i+1] != dataset.seg_tmp$ID_IARC[i]))){paste0(i,",NA")}else{i})), collapse = ",")) #Segments without changes at the right of the segment which cut the gene
  events.cut = sapply(names(events.cut), function(g) paste0(unlist(sapply(unlist(strsplit(events.cut[g],",")), function(x) if(ref.33.gene$start[which(ref.33.gene$id == g)]<ifelse(x=="NA",0,dataset.seg_tmp$start[as.numeric(x)]) & (!ifelse(x=="NA",0,as.numeric(x)-1) %in% unlist(strsplit(events.cut[g],",")) | (ifelse(x=="NA",0,as.numeric(x)-1) %in% unlist(strsplit(events.cut[g],",")) & ifelse(x=="NA",0,dataset.seg_tmp$ID_IARC[as.numeric(x)-1]) != ifelse(x=="NA",0,dataset.seg_tmp$ID_IARC[as.numeric(x)])))){paste0("NA,",x)}else{x})), collapse = ",")) #Segments without changes at the left of the segment which cut the gene
  
  events.all = sapply(names(events.com), function(g) paste0(paste0(events.com[g],"//"), events.cut[g])) # Merge the segments for which the gene is complete and those for which the gene is cut
  
  events_dup = events.all[which( events.all %in% events.all[which(duplicated(events.all))] )] # Genes for which the sequence of segment index is not unique (need to be grouped)
  
  evs = c(names(events.all)[which(!events.all %in% events_dup)], unlist(sapply(unique(events_dup[which(duplicated(events_dup))]), function(e) paste0("chr",chr,":",min(ref.33.gene$start[which(ref.33.gene$id %in% names(events_dup)[which(events_dup == e)])]),"-",max(ref.33.gene$end[which(ref.33.gene$id %in% names(events_dup)[which(events_dup == e)])])))))
  names(evs)[1:length(which(!events.all%in%events_dup))] = events.all[which(!events.all%in%events_dup)]
  
  # Make sure that genes are grouped because of genome proximity
  evs.chr = evs[which(startsWith(evs,"chr"))] 
  for(ev.chr in evs.chr){
    other = which(dataset.seg_tmp$start >= as.numeric(sapply(strsplit(sapply(strsplit(ev.chr,":"),"[[",2),"-"),"[[",1)) & dataset.seg_tmp$end <= as.numeric(sapply(strsplit(sapply(strsplit(ev.chr,":"),"[[",2),"-"),"[[",2)))
    other = other[which(other %in% unlist(strsplit(unlist(strsplit(names(evs),"//")),",")) & !other %in% unlist(strsplit(unlist(strsplit(names(evs.chr)[which(evs.chr == ev.chr)],"//")),",")))]
    
    if(length(other)!=0){
      start = ref.33.gene$start[which(ref.33.gene$id %in% names(events.all)[which(events.all==names(evs)[which(evs == ev.chr)])])]
      names(start) = ref.33.gene$id[which(ref.33.gene$id %in% names(events.all)[which(events.all==names(evs)[which(evs == ev.chr)])])]
      
      for(o in other){
        events.all[which(names(events.all) %in% names(start)[which(start < dataset.seg_tmp$start[o])])] = paste0("T,",events.all[which(names(events.all) %in% names(start)[which(start < dataset.seg_tmp$start[o])])])
        events.all[which(names(events.all) %in% names(start)[which(start > dataset.seg_tmp$start[o])])] = paste0("F,",events.all[which(names(events.all) %in% names(start)[which(start > dataset.seg_tmp$start[o])])])
      }
    }
  }
  events_dup = events.all[which( events.all %in% events.all[which(duplicated(events.all))] )]
  
  evs = c(names(events.all)[which(!events.all %in% events_dup)], unlist(sapply(unique(events_dup[which(duplicated(events_dup))]), function(e) paste0(names(events_dup)[which(events_dup == e)], collapse = ","))))
  names(evs)[1:length(which(!events.all %in% events_dup))] = events.all[which(!events.all %in% events_dup)]
  
  D_cnv_tmp = matrix(NA, length(evs), length(samples_CNV), dimnames = list(evs, samples_CNV))
  for(ev in evs){ 
    g = sapply(strsplit(ev,","),"[[",1)
    id = dataset.seg_tmp$ID_IARC[as.numeric(unlist(strsplit(strsplit(names(evs)[which(evs == ev)],"//")[[1]][2],",")))]
    dataset.seg_gene = dataset.seg_tmp[which(dataset.seg_tmp$id == g | startsWith(dataset.seg_tmp$id,paste0(g,",")) | 
                                               grepl(paste0(",",g,","),dataset.seg_tmp$id) | endsWith(dataset.seg_tmp$id,paste0(",",g)) ), c("ID_IARC","purity","copyNumber.corrected")]
    
    stopifnot(length(unique(dataset.seg_gene$ID_IARC[which(dataset.seg_gene$copyNumber.corrected != 2)])) > 2)
    if(length(unique(dataset.seg_gene$ID_IARC))>1){
      D_cnv_tmp[ev,sort(unique(dataset.seg_gene$ID_IARC))] = sapply(sort(unique(dataset.seg_gene$ID_IARC)), function(x) if(x %in% id & length(which(dataset.seg_gene$ID_IARC==x))==1){mean(c(dataset.seg_gene$purity[which(dataset.seg_gene$ID_IARC == x)]*dataset.seg_gene$copyNumber.corrected[which(dataset.seg_gene$ID_IARC == x)]+(1-dataset.seg_gene$purity[which(dataset.seg_gene$ID_IARC == x)])*2,2))}
                                                                    else{mean(dataset.seg_gene$purity[which(dataset.seg_gene$ID_IARC == x)]*dataset.seg_gene$copyNumber.corrected[which(dataset.seg_gene$ID_IARC == x)]+(1-dataset.seg_gene$purity[which(dataset.seg_gene$ID_IARC == x)])*2)})
    }else{
      D_cnv_tmp[ev,] = NA
    }
  }
  D_cnv = rbind(D_cnv, D_cnv_tmp)
}

# Select the 5,000 most variable segments
vv = apply(D_cnv,1,var,na.rm=T)
cv = cumsum(sort(vv,decreasing = T))/sum(vv,na.rm=T)
D_cnv_red = D_cnv[order(match(rownames(D_cnv),names(cv))),]
D_cnv_red = D_cnv_red[1:5000,]

# Build Minor dataset
D_loh = c()
rm(chr,g,i,x,ev.chr,o,ev)
for(chr in chrs){
  dataset.seg_tmp = dataset.seg[dataset.seg$seqnames == chr,]
  dataset.seg_tmp = as.data.frame(dataset.seg_tmp %>% group_by(start) %>% arrange(ID_IARC))
  
  genes = sort(unique(unlist(strsplit(dataset.seg_tmp$id,","))))
  
  events = sapply(genes, function(g) paste(which(dataset.seg_tmp$id == g | startsWith(dataset.seg_tmp$id, paste0(g,",")) | 
                                                   grepl(paste0(",",g,","),dataset.seg_tmp$id) | endsWith(dataset.seg_tmp$id,paste0(",",g))), collapse = ","))
  events.com = sapply(names(events), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events[g],","))), function(i) if(ref.33.gene$start[which(ref.33.gene$id == g)]>=dataset.seg_tmp$start[i] & ref.33.gene$end[which(ref.33.gene$id == g)]<=dataset.seg_tmp$end[i]){i})), collapse = ",") ) # Segments where the gene is complete
  events.cut = sapply(names(events), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events[g],","))), function(i) if(!ref.33.gene$start[which(ref.33.gene$id == g)]>=dataset.seg_tmp$start[i] | !ref.33.gene$end[which(ref.33.gene$id == g)]<=dataset.seg_tmp$end[i]){i})), collapse = ",") ) # Segments where the gene is cut
  
  # For segments where the gene is cut
  events.cut = sapply(names(events.cut), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events.cut[g],","))), function(i) if(ref.33.gene$end[which(ref.33.gene$id == g)]>dataset.seg_tmp$end[i] & (!(i+1) %in% as.numeric(unlist(strsplit(events.cut[g],","))) | ((i+1) %in% as.numeric(unlist(strsplit(events.cut[g],","))) & dataset.seg_tmp$ID_IARC[i+1] != dataset.seg_tmp$ID_IARC[i]))){paste0(i,",NA")}else{i})), collapse = ",")) # Segments without changes at the right of the segment which cut the gene
  events.cut = sapply(names(events.cut), function(g) paste0(unlist(sapply(unlist(strsplit(events.cut[g],",")), function(x) if(ref.33.gene$start[which(ref.33.gene$id == g)]<ifelse(x=="NA",0,dataset.seg_tmp$start[as.numeric(x)]) & (!ifelse(x=="NA",0,as.numeric(x)-1) %in% unlist(strsplit(events.cut[g],",")) | (ifelse(x=="NA",0,as.numeric(x)-1) %in% unlist(strsplit(events.cut[g],",")) & ifelse(x=="NA",0,dataset.seg_tmp$ID_IARC[as.numeric(x)-1]) != ifelse(x=="NA",0,dataset.seg_tmp$ID_IARC[as.numeric(x)])))){paste0("NA,",x)}else{x})), collapse = ",")) # Segments without changes at the left of the segment which cut the gene
  
  events.all = sapply(names(events.com), function(g) paste0(paste0(events.com[g],"//"), events.cut[g])) # Merge the segments for which the gene is complete and those for which the gene is cut
  
  events_dup = events.all[which( events.all %in% events.all[which(duplicated(events.all))] )] # Genes for which the sequence of segment index is not unique (need to be grouped)
  
  evs = c(names(events.all)[which(!events.all %in% events_dup)],unlist(sapply(unique(events_dup[which(duplicated(events_dup))]), function(e) paste0("chr",chr,":",min(ref.33.gene$start[which(ref.33.gene$id %in% names(events_dup)[which(events_dup == e)])]),"-",max(ref.33.gene$end[which(ref.33.gene$id %in% names(events_dup)[which(events_dup == e)])])))))
  names(evs)[1:length(which(!events.all%in%events_dup))] = events.all[which(!events.all%in%events_dup)]
  
  # Make sure that genes are grouped because of genome proximity
  evs.chr = evs[which(startsWith(evs,"chr"))] 
  for(ev.chr in evs.chr){
    other = which(dataset.seg_tmp$start >= as.numeric(sapply(strsplit(sapply(strsplit(ev.chr,":"),"[[",2),"-"),"[[",1)) & dataset.seg_tmp$end <= as.numeric(sapply(strsplit(sapply(strsplit(ev.chr,":"),"[[",2),"-"),"[[",2)))
    other = other[which(other %in% unlist(strsplit(unlist(strsplit(names(evs),"//")),",")) & !other %in% unlist(strsplit(unlist(strsplit(names(evs.chr)[which(evs.chr == ev.chr)],"//")),",")))]
    
    if(length(other)!=0){
      start = ref.33.gene$start[which(ref.33.gene$id %in% names(events.all)[which(events.all==names(evs)[which(evs == ev.chr)])])]
      names(start) = ref.33.gene$id[which(ref.33.gene$id %in% names(events.all)[which(events.all==names(evs)[which(evs == ev.chr)])])]
      
      for(o in other){
        events.all[which(names(events.all) %in% names(start)[which(start < dataset.seg_tmp$start[o])])] = paste0("T,",events.all[which(names(events.all) %in% names(start)[which(start < dataset.seg_tmp$start[o])])])
        events.all[which(names(events.all) %in% names(start)[which(start > dataset.seg_tmp$start[o])])] = paste0("F,",events.all[which(names(events.all) %in% names(start)[which(start > dataset.seg_tmp$start[o])])])
      }
    }
  }
  events_dup = events.all[which( events.all %in% events.all[which(duplicated(events.all))] )]
  
  evs = c(names(events.all)[which(!events.all %in% events_dup)], unlist(sapply(unique(events_dup[which(duplicated(events_dup))]), function(e) paste0(names(events_dup)[which(events_dup == e)], collapse = ","))))
  names(evs)[1:length(which(!events.all %in% events_dup))] = events.all[which(!events.all %in% events_dup)]
  
  D_loh_tmp = matrix(NA, length(evs), length(samples_CNV), dimnames = list(evs, samples_CNV))
  for(ev in evs){ 
    g = sapply(strsplit(ev,","),"[[",1)
    id = dataset.seg_tmp$ID_IARC[as.numeric(unlist(strsplit(strsplit(names(evs)[which(evs == ev)],"//")[[1]][2],",")))]
    dataset.seg_gene = dataset.seg_tmp[which(dataset.seg_tmp$id == g | startsWith(dataset.seg_tmp$id,paste0(g,",")) | 
                                               grepl(paste0(",",g,","),dataset.seg_tmp$id) | endsWith(dataset.seg_tmp$id,paste0(",",g)) ), c("ID_IARC","purity","minorAlleleCopyNumber.corrected")]
    
    stopifnot(length(unique(dataset.seg_gene$ID_IARC[which(dataset.seg_gene$minorAlleleCopyNumber.corrected != 1)])) > 2)
    if(length(unique(dataset.seg_gene$ID_IARC))>1){
      D_loh_tmp[ev,sort(unique(dataset.seg_gene$ID_IARC))] = sapply(sort(unique(dataset.seg_gene$ID_IARC)), function(x) if(x %in% id & length(which(dataset.seg_gene$ID_IARC==x))==1){mean(c(dataset.seg_gene$purity[which(dataset.seg_gene$ID_IARC == x)]*dataset.seg_gene$minorAlleleCopyNumber.corrected[which(dataset.seg_gene$ID_IARC == x)]+(1-dataset.seg_gene$purity[which(dataset.seg_gene$ID_IARC == x)])*1,1))}
                                                                    else{mean(dataset.seg_gene$purity[which(dataset.seg_gene$ID_IARC == x)]*dataset.seg_gene$minorAlleleCopyNumber.corrected[which(dataset.seg_gene$ID_IARC == x)]+(1-dataset.seg_gene$purity[which(dataset.seg_gene$ID_IARC == x)])*1)}) # (slide 24)
    }else{
      D_loh_tmp[ev,] = NA
    }
  }
  D_loh = rbind(D_loh, D_loh_tmp)
}

# Select the 5,000 most variable segments
vv = apply(D_loh,1,var,na.rm=T)
cv = cumsum(sort(vv,decreasing = T))/sum(vv,na.rm=T)
D_loh_red = D_loh[order(match(rownames(D_loh),names(cv))),]
D_loh_red = D_loh_red[1:5000,]

## Pre-processing of DNA gene alterations data
# From SNVs
var.annovar.maf = read.table("/path/to/damaging/SNVs/var_annovar_maf_corr.txt",header = T,sep = "\t")
dataset = var.annovar.maf[which(var.annovar.maf$Variant_Classification != "Unknown"),]
id = sort(unique(paste0("MESO_",sapply(strsplit(dataset$Tumor_Sample_Barcode,"_"),"[[",2))))

samples_mut = sort(unique(dataset$Tumor_Sample_Barcode))
genes = unique(dataset$Gene.id1)
D_mut = matrix(0, length(genes), length(samples_mut), dimnames = list(genes, samples_mut))
for(i in 1:length(genes)){
  D_mut[i,] = sapply(colnames(D_mut), function(x) if(x %in% dataset$Tumor_Sample_Barcode[which(dataset$Gene.id1 == genes[i])]){1}else{0})
}
D_mut_all = cbind(D_mut, matrix(0, nrow(D_mut), 2, dimnames = list(rownames(D_mut), c("MESO_001_T","MESO_090_T"))))

# From SVs
dataset.sv = read.table("/path/to/damaging/SVs/dataset.sv_mofa.correct.txt",header = T,sep = "\t")
samples_sv = sort(unique(dataset.sv$Sample_name))
genes = unique(dataset.sv$ID)

D_sv = matrix(0, length(genes), length(samples_sv), dimnames = list(genes, samples_sv))
for(i in 1:length(genes)){
  D_sv[i,] = sapply(colnames(D_sv), function(x) if(genes[i] %in% dataset.sv$ID[which(dataset.sv$Sample_name == x)]){1}else{0})
}

# Merge SNVs and SVs
genes = unique(c(rownames(D_sv), rownames(D_mut_all)))
D_alt = matrix(0, length(genes), length(samples_sv), dimnames = list(genes, samples_sv))
D_mut_all = D_mut_all[,order(match(colnames(D_mut_all), samples_sv))]

for(i in 1:length(genes)){
  D_alt[i,] = sapply(colnames(D_alt), function(x) if(ifelse(genes[i] %in% rownames(D_mut_all), D_mut_all[which(rownames(D_mut_all) == genes[i]), which(colnames(D_mut_all) == x)] == 1,F) | ifelse(genes[i] %in% rownames(D_sv), D_sv[which(rownames(D_sv) == genes[i]), which(colnames(D_sv) == x)] == 1,F)){1}else{0})
}
D_alt = D_alt[which(rowSums(D_alt) > 2),]

## Homogeneisation
path.MOFA = "/path/to/preprocessed/matrices/"
All.samples = sort(unique(c(colnames(D_expr_redB), colnames(D_met.proB), colnames(D_alt))))

D_alt_MOFA = matrix(NA,nrow(D_alt),length(All.samples[which(!All.samples%in%colnames(D_alt))]),dimnames = list(rownames(D_alt),All.samples[which(!All.samples%in%colnames(D_alt))]))
D_alt_MOFA = as.data.frame(D_alt_MOFA)
D_alt_MOFA = cbind(D_alt,D_alt_MOFA)
D_alt_MOFA = D_alt_MOFA[,order(match(colnames(D_alt_MOFA),All.samples))]
D_alt_MOFA = as.matrix(D_alt_MOFA)
save(D_alt_MOFA, paste0(path.MOFA,"D_alt_MOFA.RData"))

D_loh_MOFA = matrix(NA,nrow(D_loh_red),length(All.samples[which(!All.samples%in%colnames(D_loh_red))]),dimnames = list(rownames(D_loh_red),All.samples[which(!All.samples%in%colnames(D_loh_red))]))
D_loh_MOFA = as.data.frame(D_loh_MOFA)
D_loh_MOFA = cbind(D_loh_red,D_loh_MOFA)
D_loh_MOFA = D_loh_MOFA[,order(match(colnames(D_loh_MOFA),All.samples))]
D_loh_MOFA = as.matrix(D_loh_MOFA)
save(D_loh_MOFA, paste0(path.MOFA,"D_loh_MOFA.RData"))

D_cnv_MOFA = matrix(NA,nrow(D_cnv_red),length(All.samples[which(!All.samples%in%colnames(D_cnv_red))]),dimnames = list(rownames(D_cnv_red),All.samples[which(!All.samples%in%colnames(D_cnv_red))]))
D_cnv_MOFA = as.data.frame(D_cnv_MOFA)
D_cnv_MOFA = cbind(D_cnv_red,D_cnv_MOFA)
D_cnv_MOFA = D_cnv_MOFA[,order(match(colnames(D_cnv_MOFA),All.samples))]
D_cnv_MOFA = as.matrix(D_cnv_MOFA)
save(D_cnv_MOFA, paste0(path.MOFA,"D_cnv_MOFA.RData"))

D_exprB_MOFA = matrix(NA,nrow(D_expr_redB),length(All.samples[which(!All.samples%in%colnames(D_expr_redB))]),dimnames = list(rownames(D_expr_redB),All.samples[which(!All.samples%in%colnames(D_expr_redB))]))
D_exprB_MOFA = as.data.frame(D_exprB_MOFA)
D_exprB_MOFA = cbind(D_expr_redB,D_exprB_MOFA)
D_exprB_MOFA = D_exprB_MOFA[,order(match(colnames(D_exprB_MOFA),All.samples))]
D_exprB_MOFA = as.matrix(D_exprB_MOFA)
save(D_exprB_MOFA, paste0(path.MOFA,"D_expr_MOFA.RData"))

D_met.proB_MOFA = matrix(NA,nrow(D_met.proB),length(All.samples[which(!All.samples%in%colnames(D_met.proB))]),dimnames = list(rownames(D_met.proB),All.samples[which(!All.samples%in%colnames(D_met.proB))]))
D_met.proB_MOFA = as.data.frame(D_met.proB_MOFA)
D_met.proB_MOFA = cbind(D_met.proB,D_met.proB_MOFA)
D_met.proB_MOFA = D_met.proB_MOFA[,order(match(colnames(D_met.proB_MOFA),All.samples))]
D_met.proB_MOFA = as.matrix(D_met.proB_MOFA)
save(D_met.proB_MOFA, paste0(path.MOFA,"D_met.proB_MOFA.RData"))

D_met.bodB_MOFA = matrix(NA,nrow(D_met.bodB),length(All.samples[which(!All.samples%in%colnames(D_met.bodB))]),dimnames = list(rownames(D_met.bodB),All.samples[which(!All.samples%in%colnames(D_met.bodB))]))
D_met.bodB_MOFA = as.data.frame(D_met.bodB_MOFA)
D_met.bodB_MOFA = cbind(D_met.bodB,D_met.bodB_MOFA)
D_met.bodB_MOFA = D_met.bodB_MOFA[,order(match(colnames(D_met.bodB_MOFA),All.samples))]
D_met.bodB_MOFA = as.matrix(D_met.bodB_MOFA)
save(D_met.bodB_MOFA, paste0(path.MOFA,"D_met.bodB_MOFA.RData"))

D_met.enhB_MOFA = matrix(NA,nrow(D_met.enhB),length(All.samples[which(!All.samples%in%colnames(D_met.enhB))]),dimnames = list(rownames(D_met.enhB),All.samples[which(!All.samples%in%colnames(D_met.enhB))]))
D_met.enhB_MOFA = as.data.frame(D_met.enhB_MOFA)
D_met.enhB_MOFA = cbind(D_met.enhB,D_met.enhB_MOFA)
D_met.enhB_MOFA = D_met.enhB_MOFA[,order(match(colnames(D_met.enhB_MOFA),All.samples))]
D_met.enhB_MOFA = as.matrix(D_met.enhB_MOFA)
save(D_met.enhB_MOFA, paste0(path.MOFA,"D_met.enhB_MOFA.RData"))


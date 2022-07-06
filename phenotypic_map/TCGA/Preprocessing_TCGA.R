##################################################################
############# Pre-processing TCGA cohort #########################
##################################################################

require(data.table)
require(openxlsx)
require(DESeq2)
require(rtracklayer)

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
sexcpgs = rownames(Locations)[Locations$chr%in%c("chrX","chrY","chrM")]

## Load annotation
# Gene annotation GENCODE v19 and v37
ref.19 = rtracklayer::import("/path/to/gencode/annotation/file/gencode.v19.annotation.gtf.gz", format="gtf")
ref.19 = as.data.frame(ref.19)
ref.19.gene = ref.19[which(ref.19$type == "gene"),]
ref.19.gene$seqnames = gsub("^chr", "", ref.19.gene$seqnames)
ref.19.gene = as.data.table(ref.19.gene)
setnames(ref.19.gene,c("gene_name","gene_id"),c("geneSymbol","id"))

ref.37 = rtracklayer::import("/path/to/gencode/annotation/file/gencode.v33lift37.annotation.gtf.gz", format="gtf")
ref.37 = as.data.frame(ref.37)
ref.37.gene = ref.37[which(ref.37$type == "gene"),]
ref.37.gene$seqnames = gsub("^chr", "", ref.37.gene$seqnames)
ref.37.gene = as.data.table(ref.37.gene)
setnames(ref.37.gene,c("gene_name","gene_id"),c("geneSymbol","id"))

# Samples annotation
load(file = "path/to/groupsBuenoTCGAred_LM.RData") # merged clinical data for Bueno and TCGA cohorts 
Data.TCGA = read.xlsx("/path/to/supplementary/table1/205173_2_supp_5073239_pg14sl.xlsx", sheet = 2 ,startRow = 1)

# CpGs annotation
Annotation = read.delim("/path/to/EPIC850K/array/manifest/b5/Infinium_MethylationEPIC_v1b5_manifest.txt")

## Pre-processing of RNA-seq data
# Load raw read counts
expr_coun = read.csv("/path/to/gene_count_matrix_1pass.csv",row.names = 1)
colnames(expr_coun) = sapply(colnames(expr_coun), function(x) paste0(strsplit(x,"[.]")[[1]], collapse = "-"))
colnames(expr_coun) = gsub('.{1}$', '', colnames(expr_coun))
expr_coun = expr_coun[,which(colnames(expr_coun) %in% groupsBuenoTCGAred$Sample)]

## Normalisation
deseqexpr = DESeqDataSetFromMatrix(expr_coun, colData = data.frame(colnames(expr_coun)), design = ~ 1, tidy = F)
vstexpr = varianceStabilizingTransformation(deseqexpr, blind = T)

expr_annot = read.table("/path/to/TCGA-3H-AB3K-01A_pass1_gene_abund.tab",header = T,sep = "\t")[,1:6]
annot_ordered = expr_annot[sapply(rownames(vstexpr), function(x) which(expr_annot$Gene.ID==x)[1]),]
vstexpr_nosex = vstexpr[!annot_ordered$Reference %in% c("chrM","chrX","chrY"),]
vstexpr_nosex = assay(vstexpr_nosex)

# Minimal level filter based on FPKM
expr_fpkm = read.csv("/path/to/gene_FPKM_matrix.csv",row.names = 1)
colnames(expr_fpkm) = sapply(strsplit(colnames(expr_fpkm),"FPKM."),"[[",2)
colnames(expr_fpkm) = sapply(colnames(expr_fpkm), function(x) paste0(strsplit(x,"[.]")[[1]], collapse = "-"))
colnames(expr_fpkm) = gsub('.{1}$', '', colnames(expr_fpkm))
expr_fpkm = expr_fpkm[,which(colnames(expr_fpkm) %in% groupsBuenoTCGAred$Sample)]

FPKM.diff = apply(expr_fpkm, 1, max, na.rm=TRUE) - apply(expr_fpkm, 1, min, na.rm=TRUE)
names(FPKM.diff) = rownames(expr_fpkm)
FPKM.diff = FPKM.diff[which(FPKM.diff >= 1)]

vstexpr_nosex_red = vstexpr_nosex[which(rownames(vstexpr_nosex) %in% names(FPKM.diff)),]

# Select the 5,000 most variable genes
vv = apply(vstexpr_nosex_red,1,var)
cv = cumsum(sort(vv,decreasing = T))/sum(vv)
vstexpr_nosex_red = vstexpr_nosex_red[order(match(rownames(vstexpr_nosex_red),names(cv))),]
D_expr_red = vstexpr_nosex_red[1:5000,]

## Pre-processing of DNA methylation array data
# Load M-values matrix
D_met_all = read.table("/path/to/NormFilteredMTable_noInf_TCGA_87.csv",header = TRUE, sep = ",")
rownames(D_met_all) = D_met_all$X
D_met_all = D_met_all[,-1]
D_met_all = D_met_all[,which(startsWith(colnames(D_met_all),"TCGA"))]
colnames(D_met_all) = sapply(colnames(D_met_all), function(c) paste0(strsplit(c,"[.]")[[1]], collapse = "-"))
D_met_all = D_met_all[,which(colnames(D_met_all) %in% groupsBuenoTCGAred$Sample)]

# Load beta-values matrix
D_beta_all = read.table("/path/to/NormalisedFilteredBetaTable_TCGA_87.csv",header = TRUE, sep = ",")
rownames(D_beta_all) = D_beta_all$X
D_beta_all = D_beta_all[,-1]
D_beta_all = D_beta_all[,which(startsWith(colnames(D_beta_all),"TCGA"))]
colnames(D_beta_all) = sapply(colnames(D_beta_all), function(c) paste0(strsplit(c,"[.]")[[1]], collapse = "-"))
D_beta_all = D_beta_all[,which(colnames(D_beta_all) %in% groupsBuenoTCGAred$Sample)]

# Minimal level filter based on beta-values
Beta.diff = apply(D_beta_all,1,max,na.rm=TRUE) - apply(D_beta_all,1,min,na.rm=TRUE)
names(Beta.diff) = rownames(D_beta_all)
Beta.diff = Beta.diff[which(Beta.diff >= 0.1)]
D_met_red = D_met_all[which(rownames(D_met_all) %in% names(Beta.diff)),]

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
load(file = "/path/to/dataset.tcga_all.RData") # GDC TCGA MESO portal
dataset.tcga = dataset.tcga[which(!dataset.tcga$Chromosome %in% c("chrX","chrY")),]

dataset.tcga = as.data.table(dataset.tcga)
setnames(dataset.tcga,c("Chromosome","Start","End"),c("seqnames","start","end"))
dataset.tcga$seqnames = gsub("^chr", "", dataset.tcga$seqnames)

rm(dataset.tcga_tmp)
for(i in 1:nrow(dataset.tcga)){ 
  rm(d1,ro)
  if(!exists("dataset.tcga_tmp")){
    d1 = dataset.tcga[i,]
    ro = roverlaps::roverlaps(d1, ref.19.gene)
    ro[, geneSymbol := ref.19.gene$geneSymbol[subject.id]]
    d1[, geneSymbol := paste(unique(ro$geneSymbol), collapse=",")]
    ro[, id := ref.19.gene$id[subject.id]]
    d1[, id := paste(unique(ro$id), collapse=",")]
    dataset.tcga_tmp = d1
  }else{
    d1 = dataset.tcga[i,]
    ro = roverlaps::roverlaps(d1, ref.19.gene)
    ro[, geneSymbol := ref.19.gene$geneSymbol[subject.id]]
    d1[, geneSymbol := paste(unique(ro$geneSymbol), collapse=",")]
    ro[, id := ref.19.gene$id[subject.id]]
    d1[, id := paste(unique(ro$id), collapse=",")]
    temp_data = d1
    
    dataset.tcga_tmp = rbind(dataset.tcga_tmp, temp_data)
    rm(temp_data)
  }
}
dataset.tcga = dataset.tcga_tmp
dataset.tcga$purity = sapply(dataset.tcga$ID, function(x) Data.TCGA$purity[which(Data.TCGA$TCGA_barcode == x)])

dataset.tcga = as.data.frame(dataset.tcga)
samples_CNV = sort(unique(dataset.tcga$ID))
chrs = 1:22

# Build Total dataset
D_cnv = c()
rm(chr,g,i,x,ev.chr,o,ev)
for(chr in chrs){
  dataset.tcga_tmp = dataset.tcga[dataset.tcga$seqnames == chr,]
  dataset.tcga_tmp = as.data.frame(dataset.tcga_tmp %>% group_by(start) %>% arrange(ID))
  
  genes = sort(unique(unlist(strsplit(dataset.tcga_tmp$id,","))))
  
  events = sapply(genes, function(g) paste(which(dataset.tcga_tmp$id == g | startsWith(dataset.tcga_tmp$id, paste0(g,",")) | 
                                                   grepl(paste0(",",g,","),dataset.tcga_tmp$id) | endsWith(dataset.tcga_tmp$id,paste0(",",g))), collapse = ","))
  events.com = sapply(names(events), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events[g],","))), function(i) if(ref.19.gene$start[which(ref.19.gene$id == g)]>=dataset.tcga_tmp$start[i] & ref.19.gene$end[which(ref.19.gene$id == g)]<=dataset.tcga_tmp$end[i]){i})), collapse = ",") ) # Segments where the gene is complete
  events.cut = sapply(names(events), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events[g],","))), function(i) if(!ref.19.gene$start[which(ref.19.gene$id == g)]>=dataset.tcga_tmp$start[i] | !ref.19.gene$end[which(ref.19.gene$id == g)]<=dataset.tcga_tmp$end[i]){i})), collapse = ",") ) # Segments where the gene is cut
  
  # For segments where the gene is cut
  events.cut = sapply(names(events.cut), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events.cut[g],","))), function(i) if(ref.19.gene$end[which(ref.19.gene$id == g)]>dataset.tcga_tmp$end[i] & (!(i+1) %in% as.numeric(unlist(strsplit(events.cut[g],","))) | ((i+1) %in% as.numeric(unlist(strsplit(events.cut[g],","))) & dataset.tcga_tmp$ID[i+1] != dataset.tcga_tmp$ID[i]))){paste0(i,",NA")}else{i})), collapse = ",")) # Segments without changes at the right of the segment which cut the gene
  events.cut = sapply(names(events.cut), function(g) paste0(unlist(sapply(unlist(strsplit(events.cut[g],",")), function(x) if(ref.19.gene$start[which(ref.19.gene$id == g)]<ifelse(x=="NA",0,dataset.tcga_tmp$start[as.numeric(x)]) & (!ifelse(x=="NA",0,as.numeric(x)-1) %in% unlist(strsplit(events.cut[g],",")) | (ifelse(x=="NA",0,as.numeric(x)-1) %in% unlist(strsplit(events.cut[g],",")) & ifelse(x=="NA",0,dataset.tcga_tmp$ID[as.numeric(x)-1]) != ifelse(x=="NA",0,dataset.tcga_tmp$ID[as.numeric(x)])))){paste0("NA,",x)}else{x})), collapse = ",")) # Segments without changes at the left of the segment which cut the gene
  
  events.all = sapply(names(events.com), function(g) paste0(paste0(events.com[g],"//"), events.cut[g])) # Merge the segments for which the gene is complete and those for which the gene is cut
  
  events_dup = events.all[which( events.all %in% events.all[which(duplicated(events.all))] )] # Genes for which the sequence of segment indexs is not unique (need to be grouped)
  
  evs = c(names(events.all)[which(!events.all %in% events_dup)],unlist(sapply(unique(events_dup[which(duplicated(events_dup))]), function(e) paste0("chr",chr,":",min(ref.19.gene$start[which(ref.19.gene$id %in% names(events_dup)[which(events_dup == e)])]),"-",max(ref.19.gene$end[which(ref.19.gene$id %in% names(events_dup)[which(events_dup == e)])])))))
  names(evs)[1:length(which(!events.all%in%events_dup))] = events.all[which(!events.all%in%events_dup)]
  
  # Make sure that genes are grouped because of genome proximity
  evs.chr = evs[which(startsWith(evs,"chr"))] 
  for(ev.chr in evs.chr){
    other = which(dataset.tcga_tmp$start >= as.numeric(sapply(strsplit(sapply(strsplit(ev.chr,":"),"[[",2),"-"),"[[",1)) & dataset.tcga_tmp$end <= as.numeric(sapply(strsplit(sapply(strsplit(ev.chr,":"),"[[",2),"-"),"[[",2)))
    other = other[which(other %in% unlist(strsplit(unlist(strsplit(names(evs),"//")),",")) & !other %in% unlist(strsplit(unlist(strsplit(names(evs.chr)[which(evs.chr == ev.chr)],"//")),",")))] #(slide 30)
    
    if(length(other)!=0){
      start = ref.19.gene$start[which(ref.19.gene$id %in% names(events.all)[which(events.all==names(evs)[which(evs == ev.chr)])])]
      names(start) = ref.19.gene$id[which(ref.19.gene$id %in% names(events.all)[which(events.all==names(evs)[which(evs == ev.chr)])])]
      
      for(o in other){
        events.all[which(names(events.all) %in% names(start)[which(start < dataset.tcga_tmp$start[o])])] = paste0("T,",events.all[which(names(events.all) %in% names(start)[which(start < dataset.tcga_tmp$start[o])])])
        events.all[which(names(events.all) %in% names(start)[which(start > dataset.tcga_tmp$start[o])])] = paste0("F,",events.all[which(names(events.all) %in% names(start)[which(start > dataset.tcga_tmp$start[o])])])
      }
    }
  }
  events_dup = events.all[which( events.all %in% events.all[which(duplicated(events.all))] )]
  
  evs = c(names(events.all)[which(!events.all %in% events_dup)], unlist(sapply(unique(events_dup[which(duplicated(events_dup))]), function(e) paste0(names(events_dup)[which(events_dup == e)], collapse = ","))))
  names(evs)[1:length(which(!events.all %in% events_dup))] = events.all[which(!events.all %in% events_dup)]
  
  D_cnv_tmp = matrix(NA, length(evs), length(samples_CNV), dimnames = list(evs, samples_CNV))
  for(ev in evs){ 
    g = sapply(strsplit(ev,","),"[[",1)
    id = dataset.tcga_tmp$ID[as.numeric(unlist(strsplit(strsplit(names(evs)[which(evs == ev)],"//")[[1]][2],",")))]
    dataset.tcga_gene = dataset.tcga_tmp[which(dataset.tcga_tmp$id == g | startsWith(dataset.tcga_tmp$id,paste0(g,",")) | 
                                                 grepl(paste0(",",g,","),dataset.tcga_tmp$id) | endsWith(dataset.tcga_tmp$id,paste0(",",g)) ), c("ID","purity","Copy_Number")]
    
    dataset.tcga_gene$purity[which(dataset.tcga_gene$purity > 1)] = 1
    
    stopifnot(length(unique(dataset.tcga_gene$ID[which(dataset.tcga_gene$Copy_Number != 2)])) > 2)
    if(length(unique(dataset.tcga_gene$ID))>1){ 
      D_cnv_tmp[ev,sort(unique(dataset.tcga_gene$ID))] = sapply(sort(unique(dataset.tcga_gene$ID)), function(x) if(x %in% id & length(which(dataset.tcga_gene$ID==x))==1){mean(c(dataset.tcga_gene$purity[which(dataset.tcga_gene$ID == x)]*dataset.tcga_gene$Copy_Number[which(dataset.tcga_gene$ID == x)]+(1-dataset.tcga_gene$purity[which(dataset.tcga_gene$ID == x)])*2,2))}
                                                                else{mean(dataset.tcga_gene$purity[which(dataset.tcga_gene$ID == x)]*dataset.tcga_gene$Copy_Number[which(dataset.tcga_gene$ID == x)]+(1-dataset.tcga_gene$purity[which(dataset.tcga_gene$ID == x)])*2)}) 
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
  dataset.tcga_tmp = dataset.tcga[dataset.tcga$seqnames == chr,]
  dataset.tcga_tmp = as.data.frame(dataset.tcga_tmp %>% group_by(start) %>% arrange(ID))
  
  genes = sort(unique(unlist(strsplit(dataset.tcga_tmp$id,","))))
  
  events = sapply(genes, function(g) paste(which(dataset.tcga_tmp$id == g | startsWith(dataset.tcga_tmp$id, paste0(g,",")) | 
                                                   grepl(paste0(",",g,","),dataset.tcga_tmp$id) | endsWith(dataset.tcga_tmp$id,paste0(",",g))), collapse = ","))
  events.com = sapply(names(events), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events[g],","))), function(i) if(ref.19.gene$start[which(ref.19.gene$id == g)]>=dataset.tcga_tmp$start[i] & ref.19.gene$end[which(ref.19.gene$id == g)]<=dataset.tcga_tmp$end[i]){i})), collapse = ",") ) # Segments where the gene is complete
  events.cut = sapply(names(events), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events[g],","))), function(i) if(!ref.19.gene$start[which(ref.19.gene$id == g)]>=dataset.tcga_tmp$start[i] | !ref.19.gene$end[which(ref.19.gene$id == g)]<=dataset.tcga_tmp$end[i]){i})), collapse = ",") ) # Segments where the gene is cut
  
  # For segments where the gene is cut
  events.cut = sapply(names(events.cut), function(g) paste0(unlist(sapply(as.numeric(unlist(strsplit(events.cut[g],","))), function(i) if(ref.19.gene$end[which(ref.19.gene$id == g)]>dataset.tcga_tmp$end[i] & (!(i+1) %in% as.numeric(unlist(strsplit(events.cut[g],","))) | ((i+1) %in% as.numeric(unlist(strsplit(events.cut[g],","))) & dataset.tcga_tmp$ID[i+1] != dataset.tcga_tmp$ID[i]))){paste0(i,",NA")}else{i})), collapse = ",")) # Segments without changes at the right of the segment which cut the gene
  events.cut = sapply(names(events.cut), function(g) paste0(unlist(sapply(unlist(strsplit(events.cut[g],",")), function(x) if(ref.19.gene$start[which(ref.19.gene$id == g)]<ifelse(x=="NA",0,dataset.tcga_tmp$start[as.numeric(x)]) & (!ifelse(x=="NA",0,as.numeric(x)-1) %in% unlist(strsplit(events.cut[g],",")) | (ifelse(x=="NA",0,as.numeric(x)-1) %in% unlist(strsplit(events.cut[g],",")) & ifelse(x=="NA",0,dataset.tcga_tmp$ID[as.numeric(x)-1]) != ifelse(x=="NA",0,dataset.tcga_tmp$ID[as.numeric(x)])))){paste0("NA,",x)}else{x})), collapse = ",")) # Segments without changes at the left of the segment which cut the gene 
  
  events.all = sapply(names(events.com), function(g) paste0(paste0(events.com[g],"//"), events.cut[g])) # Merge the segments for which the gene is complete and those for which the gene is cut
  
  events_dup = events.all[which( events.all %in% events.all[which(duplicated(events.all))] )] # Genes for which the sequence of segment index is not unique (need to be grouped)
  
  evs = c(names(events.all)[which(!events.all %in% events_dup)],unlist(sapply(unique(events_dup[which(duplicated(events_dup))]), function(e) paste0("chr",chr,":",min(ref.19.gene$start[which(ref.19.gene$id %in% names(events_dup)[which(events_dup == e)])]),"-",max(ref.19.gene$end[which(ref.19.gene$id %in% names(events_dup)[which(events_dup == e)])])))))
  names(evs)[1:length(which(!events.all%in%events_dup))] = events.all[which(!events.all%in%events_dup)]
  
  # Make sure that genes are grouped because of genome proximity
  evs.chr = evs[which(startsWith(evs,"chr"))] 
  for(ev.chr in evs.chr){
    other = which(dataset.tcga_tmp$start >= as.numeric(sapply(strsplit(sapply(strsplit(ev.chr,":"),"[[",2),"-"),"[[",1)) & dataset.tcga_tmp$end <= as.numeric(sapply(strsplit(sapply(strsplit(ev.chr,":"),"[[",2),"-"),"[[",2)))
    other = other[which(other %in% unlist(strsplit(unlist(strsplit(names(evs),"//")),",")) & !other %in% unlist(strsplit(unlist(strsplit(names(evs.chr)[which(evs.chr == ev.chr)],"//")),",")))]
    
    if(length(other)!=0){
      start = ref.19.gene$start[which(ref.19.gene$id %in% names(events.all)[which(events.all==names(evs)[which(evs == ev.chr)])])]
      names(start) = ref.19.gene$id[which(ref.19.gene$id %in% names(events.all)[which(events.all==names(evs)[which(evs == ev.chr)])])]
      
      for(o in other){
        events.all[which(names(events.all) %in% names(start)[which(start < dataset.tcga_tmp$start[o])])] = paste0("T,",events.all[which(names(events.all) %in% names(start)[which(start < dataset.tcga_tmp$start[o])])]) 
        events.all[which(names(events.all) %in% names(start)[which(start > dataset.tcga_tmp$start[o])])] = paste0("F,",events.all[which(names(events.all) %in% names(start)[which(start > dataset.tcga_tmp$start[o])])])
      }
    }
  }
  events_dup = events.all[which( events.all %in% events.all[which(duplicated(events.all))] )]
  
  evs = c(names(events.all)[which(!events.all %in% events_dup)], unlist(sapply(unique(events_dup[which(duplicated(events_dup))]), function(e) paste0(names(events_dup)[which(events_dup == e)], collapse = ","))))
  names(evs)[1:length(which(!events.all %in% events_dup))] = events.all[which(!events.all %in% events_dup)]
  
  D_loh_tmp = matrix(NA, length(evs), length(samples_CNV), dimnames = list(evs, samples_CNV))
  for(ev in evs){ 
    g = sapply(strsplit(ev,","),"[[",1)
    id = dataset.tcga_tmp$ID[as.numeric(unlist(strsplit(strsplit(names(evs)[which(evs == ev)],"//")[[1]][2],",")))]
    dataset.tcga_gene = dataset.tcga_tmp[which(dataset.tcga_tmp$id == g | startsWith(dataset.tcga_tmp$id,paste0(g,",")) | 
                                                 grepl(paste0(",",g,","),dataset.tcga_tmp$id) | endsWith(dataset.tcga_tmp$id,paste0(",",g)) ), c("ID","purity","Minor_Copy_Number")]
    
    dataset.tcga_gene$purity[which(dataset.tcga_gene$purity > 1)] = 1
    
    stopifnot(length(unique(dataset.tcga_gene$ID[which(dataset.tcga_gene$Minor_Copy_Number != 1)])) > 2)
    if(length(unique(dataset.tcga_gene$ID))>1){ 
      D_loh_tmp[ev,sort(unique(dataset.tcga_gene$ID))] = sapply(sort(unique(dataset.tcga_gene$ID)), function(x) if(x %in% id & length(which(dataset.tcga_gene$ID==x))==1){mean(c(dataset.tcga_gene$purity[which(dataset.tcga_gene$ID == x)]*dataset.tcga_gene$Minor_Copy_Number[which(dataset.tcga_gene$ID == x)]+(1-dataset.tcga_gene$purity[which(dataset.tcga_gene$ID == x)])*1,1))}
                                                                else{mean(dataset.tcga_gene$purity[which(dataset.tcga_gene$ID == x)]*dataset.tcga_gene$Minor_Copy_Number[which(dataset.tcga_gene$ID == x)]+(1-dataset.tcga_gene$purity[which(dataset.tcga_gene$ID == x)])*1)}) 
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
som_mut = read.table("/path/to/som_mut_annot.txt", header = T, sep = "\t", stringsAsFactors = F) # TCGA mutations file (mc3.v0.2.8.PUBLIC.maf.gz) annotated with annovar

# Keep only the coding mutations
som_mut = som_mut[which(som_mut$Func.wgEncodeGencodeBasicV33lift37 %in% c("exonic","splicing","ncRNA_exonic;splicing") & som_mut$ExonicFunc.wgEncodeGencodeBasicV33lift37 != "synonymous SNV"),]

# Remove sex-chromosomes
som_mut = som_mut[which(!som_mut$Chr %in% c("X","Y")),]

# Cases where two genes are annotated in Gene.ensGene
colnames(som_mut)[6:10] = c("Func.ensGene","Gene.ensGene","GeneDetail.ensGene","ExonicFunc.ensGene","AAChange.ensGene")

# Need to use ENSEMBL IDs to make identify the recurrently altered genes
som_mut$ENS_ID = sapply(1:nrow(som_mut), function(i) if(grepl("ENST",som_mut$AAChange.ensGene[i])){strsplit(som_mut$AAChange.ensGene[i],":")[[1]][2]}else{strsplit(som_mut$GeneDetail.ensGene[i],":")[[1]][1]})
som_mut$ENS_ID[which(is.na(som_mut$ENS_ID))] = som_mut$Gene.ensGene[which(is.na(som_mut$ENS_ID))]
som_mut$ENS_ID[which(!startsWith(som_mut$ENS_ID,"ENST"))] = sapply(som_mut$ENS_ID[which(!startsWith(som_mut$ENS_ID,"ENST"))], function(x) unique(ref.37.gene$id[which(ref.37.gene$geneSymbol == x)]))

som_mut$Gene.id1 = sapply(1:nrow(som_mut), function(i) if(startsWith(som_mut$ENS_ID[i],"ENST")){unique(ref.37$gene_id[which(ref.37$transcript_id == som_mut$ENS_ID[i])])}else{som_mut$ENS_ID[i]})
som_mut$Gene1 = sapply(som_mut$Gene.id1, function(x) ref.37.gene$geneSymbol[which(ref.37.gene$id == x)])
som_mut$Gene2 = sapply(1:nrow(som_mut), function(i) if(grepl(";",som_mut$Gene.ensGene[i]) & som_mut$Func.ensGene[i] %in% c("exonic","exonic;splicing","splicing")){ifelse(length(unique(unlist(strsplit(som_mut$Gene.ensGene[i],";")))) > 1, unlist(strsplit(som_mut$Gene.ensGene[i],";"))[which(!unlist(strsplit(som_mut$Gene.ensGene[i],";")) == som_mut$Gene1[i])],NA)}else{NA})
som_mut$Gene.id2 = sapply(som_mut$Gene2, function(x) if(!is.na(x)){ref.37.gene$id[which(ref.37.gene$geneSymbol == x)][1]}else{NA})

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

samples_mut = sort(unique(som_mut$Tumor.ID))
genes = unique(som_mut$Gene.id1)
D_mut = matrix(0, length(genes), length(samples_mut), dimnames = list(genes, samples_mut))
for(i in 1:length(genes)){
  D_mut[i,] = sapply(colnames(D_mut), function(x) if(x %in% som_mut$Tumor.ID[which(som_mut$Gene.id1 == genes[i])]){1}else{0})
}
D_alt = D_mut[which(rowSums(D_mut) > 2),]

## Homogeneisation
path.MOFA = "/path/to/preprocessed/matrices/"
All.samples = sort(unique(c(colnames(D_expr_red), colnames(D_met.proB), colnames(D_cnv_red), colnames(D_alt))))

D_alt_MOFA = matrix(NA,nrow(D_alt),length(All.samples[which(!All.samples%in%colnames(D_alt))]),dimnames = list(rownames(D_alt),All.samples[which(!All.samples%in%colnames(D_alt))]))
D_alt_MOFA = as.data.frame(D_alt_MOFA)
D_alt_MOFA = cbind(D_alt,D_alt_MOFA)
D_alt_MOFA = D_alt_MOFA[,order(match(colnames(D_alt_MOFA),All.samples))]
D_alt_MOFA = as.matrix(D_alt_MOFA)
save(D_alt_MOFA, file=paste0(path.MOFA,"D_alt_MOFA.RData"))

D_loh_MOFA = matrix(NA,nrow(D_loh_red),length(All.samples[which(!All.samples%in%colnames(D_loh_red))]),dimnames = list(rownames(D_loh_red),All.samples[which(!All.samples%in%colnames(D_loh_red))]))
D_loh_MOFA = as.data.frame(D_loh_MOFA)
D_loh_MOFA = cbind(D_loh_red,D_loh_MOFA)
D_loh_MOFA = D_loh_MOFA[,order(match(colnames(D_loh_MOFA),All.samples))]
D_loh_MOFA = as.matrix(D_loh_MOFA)
save(D_loh_MOFA, file=paste0(path.MOFA,"D_loh_MOFA.RData"))

D_cnv_MOFA = matrix(NA,nrow(D_cnv_red),length(All.samples[which(!All.samples%in%colnames(D_cnv_red))]),dimnames = list(rownames(D_cnv_red),All.samples[which(!All.samples%in%colnames(D_cnv_red))]))
D_cnv_MOFA = as.data.frame(D_cnv_MOFA)
D_cnv_MOFA = cbind(D_cnv_red,D_cnv_MOFA)
D_cnv_MOFA = D_cnv_MOFA[,order(match(colnames(D_cnv_MOFA),All.samples))]
D_cnv_MOFA = as.matrix(D_cnv_MOFA)
save(D_cnv_MOFA, file=paste0(path.MOFA,"D_cnv_MOFA.RData"))

D_expr_MOFA = matrix(NA,nrow(D_expr_red),length(All.samples[which(!All.samples%in%colnames(D_expr_red))]),dimnames = list(rownames(D_expr_red ),All.samples[which(!All.samples%in%colnames(D_expr_red))]))
D_expr_MOFA = as.data.frame(D_expr_MOFA)
D_expr_MOFA = cbind(D_expr_red ,D_expr_MOFA)
D_expr_MOFA = D_expr_MOFA[,order(match(colnames(D_expr_MOFA),All.samples))]
D_expr_MOFA = as.matrix(D_expr_MOFA)
save(D_expr_MOFA, file = paste0(path.MOFA,"D_expr_MOFA.RData"))

D_met.proB_MOFA = matrix(NA,nrow(D_met.proB),length(All.samples[which(!All.samples%in%colnames(D_met.proB))]),dimnames = list(rownames(D_met.proB),All.samples[which(!All.samples%in%colnames(D_met.proB))]))
D_met.proB_MOFA = as.data.frame(D_met.proB_MOFA)
D_met.proB_MOFA = cbind(D_met.proB,D_met.proB_MOFA)
D_met.proB_MOFA = D_met.proB_MOFA[,order(match(colnames(D_met.proB_MOFA),All.samples))]
D_met.proB_MOFA = as.matrix(D_met.proB_MOFA)
save(D_met.proB_MOFA, file=paste0(path.MOFA,"D_met.proB_MOFA.RData"))

D_met.bodB_MOFA = matrix(NA,nrow(D_met.bodB),length(All.samples[which(!All.samples%in%colnames(D_met.bodB))]),dimnames = list(rownames(D_met.bodB),All.samples[which(!All.samples%in%colnames(D_met.bodB))]))
D_met.bodB_MOFA = as.data.frame(D_met.bodB_MOFA)
D_met.bodB_MOFA = cbind(D_met.bodB,D_met.bodB_MOFA)
D_met.bodB_MOFA = D_met.bodB_MOFA[,order(match(colnames(D_met.bodB_MOFA),All.samples))]
D_met.bodB_MOFA = as.matrix(D_met.bodB_MOFA)
save(D_met.bodB_MOFA, file=paste0(path.MOFA,"D_met.bodB_MOFA.RData"))

D_met.enhB_MOFA = matrix(NA,nrow(D_met.enhB),length(All.samples[which(!All.samples%in%colnames(D_met.enhB))]),dimnames = list(rownames(D_met.enhB),All.samples[which(!All.samples%in%colnames(D_met.enhB))]))
D_met.enhB_MOFA = as.data.frame(D_met.enhB_MOFA)
D_met.enhB_MOFA = cbind(D_met.enhB,D_met.enhB_MOFA)
D_met.enhB_MOFA = D_met.enhB_MOFA[,order(match(colnames(D_met.enhB_MOFA),All.samples))]
D_met.enhB_MOFA = as.matrix(D_met.enhB_MOFA)
save(D_met.enhB_MOFA, file=paste0(path.MOFA,"D_met.enhB_MOFA.RData"))


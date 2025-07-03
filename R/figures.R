library(ggplot2)
library(ggpubr)
library(ggrepel)
library(vegan)
library(gridExtra)
library(data.table)
library(devtools)

devtools::install_github("ssun6/plotmicrobiome")
library(plotmicrobiome)

setwd("/Users/ssun5/git/Error_Control_Cancer_Microbiome/")

#Figure 2
kraken_n=read.table(file="data/shuffle/reports/FT00441_Kraken2_report.txt",quote="",sep="\t")
for (n in 1:35){
  kraken_a=read.table(file=paste0("data/shuffle/reports/FT00441_",n,"_randomized_report.txt"),quote="",sep="\t")
  kraken_n=merge(kraken_n,kraken_a,all=T,by=1)
}
colnames(kraken_n)=c("taxa","real",paste("random",c(1:35),sep="_"))
kraken_n[is.na(kraken_n)]=0

tax_level=sapply(kraken_n[,1],function(i){length(strsplit(i,"\\|")[[1]])})
kraken_gen=kraken_n[grepl("g__",kraken_n[,1]) & !grepl("s__",kraken_n[,1]),]
kraken_gen_log=cbind(sapply(strsplit(kraken_gen[,1],"\\|s__"),"[[",1),log10(kraken_gen[,2:37]+1))
kraken_gen_log=kraken_gen_log[-grep("Homo$",kraken_gen_log[,1]),]

kraken_spe=kraken_n[grepl("s__",kraken_n[,1]) & !grepl("t__",kraken_n[,1]),]
kraken_spe_log=cbind(sapply(strsplit(kraken_spe[,1],"\\|t__"),"[[",1),log10(kraken_spe[,2:37]+1))
kraken_spe_log=kraken_spe_log[-grep("Homo sapiens",kraken_spe_log[,1]),]

cor_mat_gen=matrix(nrow=35,ncol=2)
for (n in 3:37){
  cor1=cor.test(kraken_gen_log[,2],kraken_gen_log[,n],method = "spearman")
  cor_mat_gen[n-2,1]=cor1$estimate
  cor_mat_gen[n-2,2]=cor1$p.value
}

cor_mat_spe=matrix(nrow=35,ncol=2)
for (n in 3:37){
  cor1=cor.test(kraken_spe_log[,2],kraken_spe_log[,n],method = "spearman")
  cor_mat_spe[n-2,1]=cor1$estimate
  cor_mat_spe[n-2,2]=cor1$p.value
}

pdf("results/Figure2a.pdf")
plot(cor_mat_gen[,1],xlab="kmer randomized",ylab="Spearman correlation with non-randomized",ylim=c(0,1))
points(cor_mat_spe[,1],col="red")
legend("topleft",c("genus","species"),col=c("black","red"),pch=1)
dev.off()


pdf("results/Figure2b.pdf",height=11,width=5)
par(mfrow=c(2,1))
plot(kraken_gen_log[,2],kraken_gen_log[,3],xlab="Unrandomized genus abundance",ylab="1-mer randomized genus abundance",pch=16)
plot(kraken_gen_log[,2],kraken_gen_log[,5],xlab="Unrandomized genus abundance",ylab="3-mer randomized genus abundance",pch=16)
dev.off()


#Figure1
#modified kraken scripts /users/ssun5/kraken2_Mar26/kraken2-2.1.3
#exported modified database at /projects/afodor_research3/ssun5/kraken_library

tax_id=fread(file="data/hash_profile/tax_id.txt")
tax_freq=read.table(file="data/hash_profile/hash_tax_freq_table.txt",sep=" ")
id_match=read.table(file="data/hash_profile/external_to_internal_map.txt",sep="\t",header = T)
id_match1=id_match[match(tax_freq[,2],id_match[,2]),1]
tax_id2=match(id_match1,unlist(tax_id[,1]))
tax_id_d=tax_id[tax_id2,]
tax_id_d1=data.frame(tax_id_d)
tax_id_d1$num_hash=tax_freq[,1]
tax_id_d1$freq_hash=tax_freq[,1]/sum(tax_freq[,1])
pdf("results/Figure1b.pdf")
hist(log10(tax_id_d1$num_hash),xlab="hash/taxa (log10)",main="")
dev.off()
write.csv(tax_id_d1,file="results/Figure1a.csv")


#Figure 3
#Fig. 3A: real pcoa
metadata1=meta_format(metadata="data/tissue_seqs/cancer_seq_meta.csv",metadata_sep=",",meta_sample_name_col=1)
taxa_table1=format_wgs(taxa_file ="data/tissue_seqs/merged_kraken.txt",method = "kraken",reads_cutoff = 0)
palette_group1=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400","darkgrey")
pdf("results/Figure3a.pdf",width=7,height=21)
mds_plot(taxa_table = taxa_table1, one_level = F, metadata = metadata1, test_metadata = "sample_type", log_norm = T, taxa_level = "genus",                
         method_mds = "pcoa", distance_type = "bray", palette_grou=palette_group1)
dev.off()

#Fig. 3B: 1mer random pcoa
pdf("results/Figure3b.pdf",width=7,height=21)
taxa_table2=format_wgs(taxa_file ="data/tissue_seqs/merged_kraken_rand1mer.txt",method = "kraken",reads_cutoff = 0)
mds_plot(taxa_table = taxa_table2, one_level = F, metadata = metadata1, test_metadata = "sample_type", log_norm = T, taxa_level = "genus",                
         method_mds = "pcoa", distance_type = "bray", palette_grou=palette_group1)
dev.off()

#Fig. 3C and D: reads and GC contents
meta_s=metadata1[intersect(colnames(taxa_table1),rownames(metadata1)),]
meta_s=data.frame(meta_s)
meta_s$reads=as.numeric(meta_s$Bases)/as.numeric(meta_s$AvgSpotLen)
aggregate(reads~sample_type,data=meta_s,mean)
aggregate(reads~sample_type,data=meta_s,sd)
aggregate(AvgSpotLen~sample_type,data=meta_s,mean)

gc=read.csv(file="data/tissue_seqs/gc_contents_summary_new.txt",header=F)

gc[,2]=as.numeric(gsub("%","",gc[,2]))
gc[,1]=gsub("_1.fastq","",gc[,1])
meta_s$GC=gc[match(rownames(meta_s),gc[,1]),2]

pdf("results/Figure3cd.pdf",height=5,width=6.5)
par(mar=c(5,15,5,5),mfrow=c(1,1))
boxplot(log10(reads)~sample_type,data=meta_s,xlab="Reads (log10)",las=1,ylab="",border=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"),col="white",horizontal = TRUE)
stripchart(log10(reads)~sample_type,data=meta_s,vertical = FALSE,  method = "jitter", add = TRUE, pch = 16,col=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"))

boxplot(GC~sample_type,data=meta_s,las=1,ylab="",xlab="GC (%)",border=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"),col="white",horizontal = TRUE)
stripchart(GC~sample_type,data=meta_s,vertical = FALSE,  method = "jitter", add = TRUE, pch = 16,col=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"))
dev.off()

#Figure 4
file1=read.table("data/tissue_seqs/merged_kraken.txt",sep="\t",quote="",header=T)
file2=read.table("data/tissue_seqs/merged_kraken_rand1mer.txt",sep="\t",quote="",header=T)

file12=merge(file1,file2,by=1,all=T)
file12[is.na(file12)]=0
file12_1=file12
#log10 transform
file12_1[,2:295]=log10(file12[,2:295]+1)
file12_1=file12_1[order(rowSums(file12[,2:295]),decreasing = T),]
#only keep the genus level
file12_1_g=file12_1[grepl("g__",file12_1[,1]) & !grepl("s__",file12_1[,1]),]
#remove human data
file12_1_g=file12_1_g[file12_1_g[,1]!="d__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates|f__Hominidae|g__Homo",]


#boxplots of correlation coefficients
cor_mat=matrix(ncol=2,nrow=147)
i=1
for (n in c(2:148)){
  cor1=cor.test(as.numeric(file12_1_g[,n]),as.numeric(file12_1_g[,n+147]),method = "spearman")
  cor_mat[i,1]=cor1$estimate
  cor_mat[i,2]=cor1$p.value
  i=i+1
}

rownames(cor_mat)=sapply(strsplit(colnames(file12_1_g[,2:148]),".x"),"[[",1)
colnames(cor_mat)=c("rho_1mer","p_1mer")
meta_s$rho_1mer=cor_mat[match(rownames(meta_s),rownames(cor_mat)),1]
meta_s$p_1mer=cor_mat[match(rownames(meta_s),rownames(cor_mat)),2]

pdf("results/Figure4a.pdf",height=7.5,width=4)
par(mar=c(15,5,5,5),mfrow=c(1,1))
boxplot(rho_1mer~sample_type,data=meta_s,las=2,xlab="",ylab="Spearman's rho",border=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"),col="white")
stripchart(rho_1mer~sample_type,data=meta_s,vertical = TRUE,  method = "jitter", add = TRUE, pch = 16,col=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"))
dev.off()


#example scatter plots
#only label the most abundant ones
tax_gen1=sapply(strsplit(file12_1_g[1:100,1],"g__"),"[[",2)
tax_gen=c(tax_gen1,rep("",nrow(file12_1_g)-100))
cor_mat=matrix(ncol=2,nrow=6)
p1_2=list()
i=1
for (n in c(2,41,49,73,105,121)){
  cor1=cor.test(as.numeric(file12_1_g[,n]),as.numeric(file12_1_g[,n+147]),method = "spearman")
  cor_mat[i,1]=cor1$estimate
  cor_mat[i,2]=cor1$p.value
  main1=paste0(gsub(".x","",colnames(file12_1_g)[n])," rho = ",round(cor1$estimate,3),"   P = ",formatC(cor1$p.value, format = "e", digits = 2))
  p1=ggplot(file12_1_g, mapping=aes_string(x=colnames(file12_1_g)[n], y=colnames(file12_1_g)[n+147])) +geom_point()+
    theme_classic(base_size = 15) + labs(title=main1,x ="real (log10)" , y ="1mer random (log10)")+xlim(0, 6) +ylim(0, 6)
  p1_2[[i]]=p1+geom_text_repel(aes(label =tax_gen),size = 5,force_pull=0,force=0.5,max.overlaps=20)+geom_abline(intercept = 0, slope = 1, linetype="dotted")
  i=i+1
}

pdf("results/Figure4b.pdf",onefile = T,width=10,height=10)
for (n in 1:6){
  print(p1_2[[n]])
}
dev.off()


#Figure 5
#correlation of classified abundance and hash abundance in database

#hash abudance table
tax_id_d1=read.csv("results/Figure1a.csv")

#classified abundance table

file1_2=file1[,2:148]
rownames(file1_2)=file1[,1]
file1_2=file1_2[grepl("g__",rownames(file1_2)) & !grepl("s__",rownames(file1_2)),]
rownames(file1_2)=sapply(strsplit(rownames(file1_2),"g__"),"[[",2)
file1_2=file1_2[rownames(file1_2)!="Homo",]

file2_2=file2[,2:148]
rownames(file2_2)=file2[,1]
file2_2=file2_2[grepl("g__",rownames(file2_2)) & !grepl("s__",rownames(file2_2)),]
rownames(file2_2)=sapply(strsplit(rownames(file2_2),"g__"),"[[",2)
file2_2=file2_2[rownames(file2_2)!="Homo",]

match1=intersect(rownames(file1_2),tax_id_d1[,4])
id1=tax_id_d1[match(match1,tax_id_d1[,4]),]
comp1=file1_2[match(match1,rownames(file1_2)),]

match2=intersect(rownames(file2_2),tax_id_d1[,4])
id2=tax_id_d1[match(match2,tax_id_d1[,4]),]
comp2=file2_2[match(match2,rownames(file2_2)),]

ids=list()
ids[[1]]=id1
ids[[2]]=id2

comps=list()
comps[[1]]=comp1
comps[[2]]=comp2


cors=list()
for (n in 1:2){
  id_n=ids[[n]]
  comp_n=comps[[n]]
  cor_n=matrix(nrow=147,ncol=2)
  colnames(cor_n)=c("rho","p")
  rownames(cor_n)=colnames(comp_n)
  for (m in 1:147){
    cor1=cor.test(id_n[,10],comp_n[,m],method="spearman")
    cor_n[m,1]=cor1$estimate
    cor_n[m,2]=cor1$p.value
  }
  cors[[n]]=cor_n
}

meta_s$rho_raw_hash=cors[[1]][,1]
meta_s$p_raw_hash=cors[[1]][,2]

meta_s$rho_1mer_hash=cors[[2]][,1]
meta_s$p_1mer_hash=cors[[2]][,2]

pdf("results/Figure5a.pdf",height=7,width=10,onefile = T)
par(mar=c(15,5,5,5),mfrow=c(1,2))
boxplot(rho_raw_hash~sample_type,data=meta_s,ylab="Spearman's rho",las=2,xlab="",main="raw",border=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"),col="white",horizontal = F)
stripchart(rho_raw_hash~sample_type,data=meta_s,vertical = T,  method = "jitter", add = TRUE, pch = 16,col=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"))

boxplot(rho_1mer_hash~sample_type,data=meta_s,ylab="Spearman's rho",las=2,xlab="",main="1mer",border=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"),col="white",horizontal = F)
stripchart(rho_1mer_hash~sample_type,data=meta_s,vertical = T,  method = "jitter", add = TRUE, pch = 16,col=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"))
dev.off()

#examples
cor1=cor.test(id1[,10],comp1[,1],method="spearman")
cor2=cor.test(id2[,10],comp2[,1],method="spearman")
pdf("results/Figure5b.pdf",height=10,width=6)
par(mfrow=c(2,1))
plot(log10(id1[,10]),log10(comp1[,1]),main=paste("kraken_default \ncor =",round(cor1$estimate,2),"p-value < 2.2e-16"),xlab="log10(num_hash)",ylab="log10(abd) SRR25046362")
plot(log10(id2[,10]),log10(comp2[,1]),main=paste("kraken_conf0.5 \ncor =",round(cor2$estimate,2),"p-value < 2.2e-16"),xlab="log10(num_hash)",ylab="log10(abd) SRR25046362")
dev.off()


#Figure 6
#add metadata groups to separate one cancer from the others
meta_s$brain_binary=meta_s$sample_type
meta_s$colon_binary=meta_s$sample_type
meta_s$crc_binary=meta_s$sample_type
meta_s$lung_binary=meta_s$sample_type
meta_s$skin_binary=meta_s$sample_type
meta_s$ovarian_binary=meta_s$sample_type

meta_s$brain_binary[meta_s$brain_binary!="Brain tumor"]="others"
meta_s$colon_binary[meta_s$colon_binary!="Colon"]="others"
meta_s$crc_binary[meta_s$crc_binary!="CRC Tumor"]="others"
meta_s$lung_binary[meta_s$lung_binary!="Lung cancer"]="others"
meta_s$skin_binary[!grepl("Mucoepidermoid carcinoma",meta_s$skin_binary)]="others"
meta_s$ovarian_binary[meta_s$ovarian_binary!="Ovarian"]="others"


#calculate difference between real and randomized data
file12=merge(file1_2,file2_2,by=0,all=T)
file12[is.na(file12)]=0
df12=file12[,2:148]-file12[,149:295]
rownames(df12)=file12[,1]
colnames(df12)=colnames(file1_2)

comp1=t(t(file1_2)/colSums(file1_2))*mean(colSums(file1_2))
comp2=t(t(file2_2)/colSums(file2_2))*mean(colSums(file2_2))
comp12=t(t(df12)/colSums(df12))*mean(colSums(df12))


list1=c("brain","colon","crc","lung","skin","ovarian")
Figure1=list()
for (i in 1:6){
  group1=paste(list1[i],"binary",sep="_")
  fdrs1=stat_test(taxa_table = comp1, metadata=meta_s,test_metadata=group1,log_norm = F,method = "wilcoxon")
  fdrs2=stat_test(taxa_table = comp2, metadata=meta_s,test_metadata=group1,log_norm = F,method = "wilcoxon")
  fdrs3=stat_test(taxa_table = comp12, metadata=meta_s,test_metadata=group1,log_norm = F,method = "wilcoxon")
  
  Figure1[[i*2-1]]=p_compare(fdrs1, fdrs2,p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=0.05,cor_method="spearman",x.reverse=F,y.reverse=F,exclude_unclassified=T,one_level=T,direction=T)
  Figure1[[i*2]]=p_compare(fdrs1, fdrs3,p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=0.05,cor_method="spearman",x.reverse=F,y.reverse=F,exclude_unclassified=T,one_level=T,direction=T)
}

#in the generated Figure, left is real vs random (Figure S3), right is real vs subtracted difference (Figure 6)
pdf("results/Figure6.pdf",height=30,width=20,onefile = T)
for (i in seq(1, length(Figure1), by = 6)) {
  do.call(grid.arrange, c(Figure1[i:min(i+5, length(Figure1))], ncol=2, nrow=3))
}
dev.off()

#high microbial biomass examples
#urban vs rural
file1=read.table("data/tissue_seqs/urban_rural/merged_kraken.txt",sep="\t",quote="",header=T,row.names = 1)
file2=read.table("data/tissue_seqs/urban_rural/merged_kraken_random1mer.txt",sep="\t",header=T,row.names = 1)

file1_2=file1[grepl("g__",rownames(file1)) & !grepl("s__",rownames(file1)),]
rownames(file1_2)=sapply(strsplit(rownames(file1_2),"g__"),"[[",2)
file1_2=file1_2[rownames(file1_2)!="Homo",]

file2_2=file2[grepl("g__",rownames(file2)) & !grepl("s__",rownames(file2)),]
rownames(file2_2)=sapply(strsplit(rownames(file2_2),"g__"),"[[",2)
file2_2=file2_2[rownames(file2_2)!="Homo",]

file12=merge(file1_2,file2_2,by=0,all=T)
file12[is.na(file12)]=0

df12=file12[,2:41]-file12[,42:81]
rownames(df12)=file12[,1]
colnames(df12)=colnames(file1)

meta1=data.frame(cbind(colnames(df12),c(rep("urban",16),"rural",rep("urban",3),"rural","urban",rep("rural",18))))
colnames(meta1)=c("sample","group")
rownames(meta1)=meta1[,1]

comp1=t(t(file1_2)/colSums(file1_2))*mean(colSums(file1_2))
comp2=t(t(file2_2)/colSums(file2_2))*mean(colSums(file2_2))
comp12=t(t(df12)/colSums(df12))*mean(colSums(df12))


fdrs1=stat_test(taxa_table = comp1, metadata=meta1,test_metadata="group",log_norm = F,method = "wilcoxon")
fdrs2=stat_test(taxa_table = comp2, metadata=meta1,test_metadata="group",log_norm = F,method = "wilcoxon")
fdrs3=stat_test(taxa_table = comp12, metadata=meta1,test_metadata="group",log_norm = F,method = "wilcoxon")
pdf("results/Figure6_urban_rural.pdf",height=10,width=10)
p_compare(fdrs1, fdrs2,p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=0.05,cor_method="spearman",x.reverse=F,y.reverse=F,exclude_unclassified=T,one_level=T,direction=T)
p_compare(fdrs1, fdrs3,p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=0.05,cor_method="spearman",x.reverse=F,y.reverse=F,exclude_unclassified=T,one_level=T,direction=T)
dev.off()

#stool vs swab
taxa_table1=read.table(file = "data/tissue_seqs/stool_swab/merged_kraken_1mer_stool_swab_default.txt",sep="\t",row.names = 1,header = T,quote="")
colnames(taxa_table1)=sapply(strsplit(colnames(taxa_table1),"_"),"[[",1) #fix sample names
taxa_table2=taxa_table1[grepl("g__",rownames(taxa_table1)) & !grepl("s__",rownames(taxa_table1)),]
rownames(taxa_table2)=sapply(strsplit(rownames(taxa_table2),"g__"),"[[",2)
taxa_table2=taxa_table2[rownames(taxa_table2)!="Homo",]

#sample type are in sample names
metadata2=data.frame(sapply(colnames(taxa_table2),function(i){substr(i,1,2)}))
colnames(metadata2)="sample_type"
metadata2$name=sapply(strsplit(colnames(taxa_table2),"_"),"[[",1)

tab1=read.csv(file="data/tissue_seqs/stool_swab/kraken2.csv",row.names = 1)
tab2=tab1[grepl("g__",rownames(tab1)) & !grepl("s__",rownames(tab1)),]
rownames(tab2)=sapply(strsplit(rownames(tab2),"g__"),"[[",2)
tab3=tab2[,match(metadata2$name,colnames(tab2))]

file12=merge(tab3,taxa_table2,by=0,all=T)
file12[is.na(file12)]=0

df12=file12[,2:191]-file12[,192:381]
rownames(df12)=file12[,1]
colnames(df12)=colnames(tab3)

comp1=t(t(tab3)/colSums(tab3))*mean(colSums(tab3))
comp2=t(t(taxa_table2)/colSums(taxa_table2))*mean(colSums(taxa_table2))
comp12=t(t(df12)/colSums(df12))*mean(colSums(df12))


fdrs1=stat_test(taxa_table = comp1, metadata=metadata2,test_metadata="sample_type",log_norm = F,method = "wilcoxon")
fdrs2=stat_test(taxa_table = comp2, metadata=metadata2,test_metadata="sample_type",log_norm = F,method = "wilcoxon")
fdrs3=stat_test(taxa_table = comp12, metadata=metadata2,test_metadata="sample_type",log_norm = F,method = "wilcoxon")
pdf("results/Figure6_stool_swab.pdf",height=10,width=10)
p_compare(fdrs1, fdrs2,p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=0.05,cor_method="spearman",x.reverse=F,y.reverse=F,exclude_unclassified=T,one_level=T,direction=T)
p_compare(fdrs1, fdrs3,p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=0.05,cor_method="spearman",x.reverse=F,y.reverse=F,exclude_unclassified=T,one_level=T,direction=T)
dev.off()


#Figure 7
meta=read.csv(file="data/conf/cancer_seq_meta_conf0.csv",row.names = 1)
conf=seq(0,0.6,0.05)

tab_ran=read.csv(file="data/conf/kraken_report_1mer.csv",row.names = 1)
genus_ran=tab_ran[grepl("g__",rownames(tab_ran)) & !grepl("s__",rownames(tab_ran)),]
genus_ran=genus_ran[,order(colnames(genus_ran))]

tab=read.csv(file="data/conf/kraken_report.csv",row.names = 1)
genus=tab[grepl("g__",rownames(tab)) & !grepl("s__",rownames(tab)),]
genus=genus[,order(colnames(genus))]

names1=c("Streptomyces","Pseudomonas","Bacillus","Bradyrhizobium","Vibrio","Paenibacillus","Enterobacteriaceae","Escherichia coli","Rhodococcus","Rhizobium")

gFigure1=list()
i=1
for (n in names1[1:4]){
  for (m in names(table(meta$sample_type_BioProject))){
    n1=sapply(strsplit(rownames(meta)[meta$sample_type_BioProject==m & meta$random_group == "shuffled"],"_1mer_report_conf_0.00.txt"),"[[",1)
    n2=paste0(rep(n1,13),"_1mer_report_conf_",rep(as.character(formatC(conf,digits=2,format = "f")),each=length(n1)),".txt")
    mat1=matrix(log10(as.numeric(genus_ran[grep(paste0("g__",n),rownames(genus_ran),fixed=T),na.omit(match(n2,colnames(genus_ran)))])+1),ncol=13,byrow = F)
    mat2=matrix(colnames(genus_ran)[match(n2,colnames(genus_ran))],ncol=13,byrow = F)
    colnames(mat1)=paste0(as.character(formatC(conf,digits=2,format = "f")),"_random")
    
    n3=paste0(rep(n1,13),"_report_conf_",rep(as.character(formatC(conf,digits=2,format = "f")),each=length(n1)),".txt")
    mat1_1=matrix(log10(as.numeric(genus[grep(paste0("g__",n),rownames(genus),fixed=T),na.omit(match(n3,colnames(genus)))])+1),ncol=13,byrow = F)
    mat2_1=matrix(colnames(genus)[match(n3,colnames(genus))],ncol=13,byrow = F)
    colnames(mat1_1)=paste0(as.character(formatC(conf,digits=2,format = "f")),"_real")
    
    mat3=cbind(mat1_1,mat1)
    mat3=mat3[,order(colnames(mat3))]
    
    mat_long=melt(mat3)
    mat_long[,2]=as.character(mat_long[,2])
    mat_long$conf=sapply(strsplit(mat_long[,2],"_"),"[[",1)
    mat_long$group=sapply(strsplit(mat_long[,2],"_"),"[[",2)
    mat_long$group=factor(mat_long$group,levels=c("real","random"))
    main1=paste0(n,"\n",m)
    gFigure1[[i]]=ggboxplot(mat_long, x = "conf", y = "value", color = "group",palette = "jco", add = "jitter")+
      ggtitle(main1)+xlab("Confidence score")+ylab("Reads (log10)")
    i=i+1
  }
}

pdf("results/Figure7.pdf",onefile = T,width = 36, height = 20)
ggarrange(plotlist = gFigure1, ncol = 6, nrow = 4)
dev.off()

#Figure S2
file1=read.table("data/tissue_seqs/merged_kraken.txt",sep="\t",quote="",header=T)
taxa_level=sapply(file1[,1],function(i){length(strsplit(i,"\\|")[[1]])})
human1=file1[grep("Homo",file1[,1])[3],]
kingdom1=file1[which(taxa_level==1),]
class1=(colSums(kingdom1[,-1])-human1[-1])/colSums(kingdom1[,-1])

file2=read.table("data/tissue_seqs/merged_kraken_conf0.5.txt",sep="\t",quote="",header=T)
taxa_level=sapply(file2[,1],function(i){length(strsplit(i,"\\|")[[1]])})
human2=file2[grep("Homo",file2[,1])[1],]
kingdom2=file2[which(taxa_level==1),]
class2=(colSums(kingdom2[,-1])-human2[-1])/colSums(kingdom2[,-1])

file3=read.table("data/tissue_seqs/merged_kraken_rand3mer.txt",sep="\t",quote="",header=T)
taxa_level=sapply(file3[,1],function(i){length(strsplit(i,"\\|")[[1]])})
human3=file3[grep("Homo",file3[,1])[3],]
kingdom3=file3[which(taxa_level==1),]
class3=(colSums(kingdom3[,-1])-human3[-1])/colSums(kingdom3[,-1])

file4=read.table("data/tissue_seqs/merged_kraken_rand3mer_conf0.5.txt",sep="\t",quote="",header=T)
taxa_level=sapply(file4[,1],function(i){length(strsplit(i,"\\|")[[1]])})
human4=file4[grep("Homo",file4[,1])[1],]
kingdom4=file4[which(taxa_level==1),]
class4=(colSums(kingdom4[,-1])-human4[-1])/colSums(kingdom4[,-1])

file5=read.table("data/tissue_seqs/merged_kraken_rand6mer.txt",sep="\t",quote="",header=T)
taxa_level=sapply(file5[,1],function(i){length(strsplit(i,"\\|")[[1]])})
human5=file5[grep("Homo",file5[,1])[3],]
kingdom5=file5[which(taxa_level==1),]
class5=(colSums(kingdom5[,-1])-human5[-1])/colSums(kingdom5[,-1])

file6=read.table("data/tissue_seqs/merged_kraken_rand6mer_conf0.5.txt",sep="\t",quote="",header=T)
taxa_level=sapply(file6[,1],function(i){length(strsplit(i,"\\|")[[1]])})
human6=file6[grep("Homo",file6[,1])[1],]
kingdom6=file6[which(taxa_level==1),]
class6=(colSums(kingdom6[,-1])-human6[-1])/colSums(kingdom6[,-1])

file7=read.table("data/tissue_seqs/merged_kraken_rand1mer.txt",sep="\t",quote="",header=T)
taxa_level=sapply(file7[,1],function(i){length(strsplit(i,"\\|")[[1]])})
human7=file7[grep("Homo",file7[,1])[3],]
kingdom7=file7[which(taxa_level==1),]
class7=(colSums(kingdom7[,-1])-human7[-1])/colSums(kingdom7[,-1])

file8=read.table("data/tissue_seqs/merged_kraken_rand1mer_conf0.5.txt",sep="\t",quote="",header=T)
taxa_level=sapply(file8[,1],function(i){length(strsplit(i,"\\|")[[1]])})
human8=file8[grep("Homo",file8[,1])[1],]
kingdom8=file8[which(taxa_level==1),]
class8=(colSums(kingdom8[,-1])-human8[-1])/colSums(kingdom8[,-1])

reads=rbind(human1[-1],colSums(kingdom1[,-1]),human2[-1],colSums(kingdom2[,-1]),human3[-1],colSums(kingdom3[,-1]),human4[,-1],colSums(kingdom4[,-1]),human5[,-1],colSums(kingdom5[,-1]),human6[,-1],colSums(kingdom6[,-1]),human7[,-1],colSums(kingdom7[,-1]),human8[,-1],colSums(kingdom8[,-1]))
rownames(reads)=paste(rep(c("kraken","kraken_conf0.5","rand3","rand3_conf0.5","rand6","rand6_conf0.5","rand1","rand1_conf0.5"),each=2),rep(c("human","classified"),8),sep="_")

all1=fread("data/tissue_seqs/reads_num.csv")
all2=t(all1)
all2=all2[,match(colnames(reads),all2[1,])]
colnames(all2)=all2[1,]
all2=all2[-1,]
reads2=rbind(reads,all2)

meta_s$classified_raw=as.numeric(reads2["kraken_classified",match(rownames(meta_s),colnames(reads2))])
meta_s$classified_human_raw=as.numeric(reads2["kraken_human",match(rownames(meta_s),colnames(reads2))])
meta_s$classified_bac_raw=as.numeric(kingdom1[2,match(rownames(meta_s),colnames(kingdom1))])
meta_s$classified_bac_perc_raw=meta_s$classified_bac_raw/meta_s$classified_raw*100
meta_s$classified_perc_raw=meta_s$classified_raw/meta_s$reads*100
meta_s$classified_human_perc_raw=meta_s$classified_human_raw/meta_s$reads*100


meta_s$classified_1mer=as.numeric(reads2["rand1_classified",match(rownames(meta_s),colnames(reads2))])
meta_s$classified_human_1mer=as.numeric(reads2["rand1_human",match(rownames(meta_s),colnames(reads2))])
meta_s$classified_bac_1mer=as.numeric(kingdom7[2,match(rownames(meta_s),colnames(kingdom7))])
meta_s$classified_bac_perc_1mer=meta_s$classified_bac_1mer/meta_s$classified_1mer*100
meta_s$classified_perc_1mer=meta_s$classified_1mer/meta_s$reads*100
meta_s$classified_human_perc_1mer=meta_s$classified_human_1mer/meta_s$reads*100


pdf("results/FigureS2.pdf",height=20,width=8)
par(mar=c(15,5,5,5),mfrow=c(4,2))
boxplot(log10(classified_raw)~sample_type,data=meta_s,las=2,xlab="",ylab="Reads (log10)",main="Classified reads (unshuffled)",border=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"),col="white")
stripchart(log10(classified_raw)~sample_type,data=meta_s,vertical = TRUE,  method = "jitter", add = TRUE, pch = 16,col=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"))

boxplot(log10(classified_1mer)~sample_type,data=meta_s,las=2,xlab="",ylab="Reads (log10)",main="Classified reads (1-mer shuffled)",border=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"),col="white")
stripchart(log10(classified_1mer)~sample_type,data=meta_s,vertical = TRUE,  method = "jitter", add = TRUE, pch = 16,col=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"))

boxplot(log10(classified_bac_raw)~sample_type,data=meta_s,las=2,xlab="",ylab="Classified bacteria reads (log10)",main="Classified bacteria reads (unshuffled)",border=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"),col="white")
stripchart(log10(classified_bac_raw)~sample_type,data=meta_s,vertical = TRUE,  method = "jitter", add = TRUE, pch = 16,col=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"))

boxplot(log10(classified_bac_1mer)~sample_type,data=meta_s,las=2,xlab="",ylab="Classified bacteria reads (log10)",main="Classified bacteria reads (1-mer shuffled)",border=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"),col="white")
stripchart(log10(classified_bac_1mer)~sample_type,data=meta_s,vertical = TRUE,  method = "jitter", add = TRUE, pch = 16,col=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"))

boxplot(classified_bac_perc_raw~sample_type,data=meta_s,las=2,xlab="",ylab="Percentage of reads classified as bacteria (%)",main="Percentage of reads classified as bacteria (unshuffled)",border=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"),col="white")
stripchart(classified_bac_perc_raw~sample_type,data=meta_s,vertical = TRUE,  method = "jitter", add = TRUE, pch = 16,col=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"))

boxplot(classified_bac_perc_1mer~sample_type,data=meta_s,las=2,xlab="",ylab="Percentage of reads classified as bacteria (%)",main="Percentage of reads classified as bacteria (1-mer shuffled)",border=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"),col="white")
stripchart(classified_bac_perc_1mer~sample_type,data=meta_s,vertical = TRUE,  method = "jitter", add = TRUE, pch = 16,col=c("#3498DB","#E74C3C","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#D35400"))
dev.off()



#Figure S3
#left columns of p_compare_cancer.pdf


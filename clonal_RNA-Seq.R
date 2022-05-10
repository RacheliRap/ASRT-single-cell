############################################################################################
# clonal RNA-seq analysis.
# The data retrived from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148348
############################################################################################


library(gprofiler2)
library(EnsDb.Mmusculus.v79)


setwd("C:/Users/rachelRa/")

ensdb_genes <- genes(EnsDb.Mmusculus.v79)
ensdb_genes.coding = ensdb_genes[ensdb_genes$gene_biotype == "protein_coding"]
X_names <- ensdb_genes[seqnames(ensdb_genes) == "X"]$gene_id

file_names <- dir("Documents/", full.names = T, pattern = "GSM44*") #where you have your files
file.names.NPC = file_names[grep("_Ch", file_names)]

table <- do.call(cbind,lapply(file.names.NPC,read.delim))
table$Gene = gsub("\\..*","", table$Gene)
table = table[!table$Gene %in% X_names,]

async.list = read.delim("Documents/async_master_list.protein_coding.gencode.vM1.bed", header = F)
async.list = async.list[async.list$V1!="chrX",]
gConvert = gconvert(query = async.list$V5, organism = "mmusculus",
         target="ENSG", mthreshold = Inf, filter_na = F)
async.list$ensbl = gConvert[match(async.list$V5, gConvert$input),"target"]

table.async = table[table$Gene %in% async.list$ensbl,]
#table.async = table

table.async.ratio = table.async[,grep("allelicRatio_CAST_EiJ", names(table.async))]
row.names(table.async.ratio) <- table.async$Gene
table.async.ratio.naOmit = table.async.ratio[complete.cases(table.async.ratio),]
keep = apply(table.async.ratio.naOmit,1, function(x) all(x>0.7|x<0.3)) 
table.async.ratio.naOmit.filter = table.async.ratio.naOmit[keep,]

prop.async = table.async.ratio.naOmit.filter
prop.async[prop.async<=0.5] <- -1
prop.async[prop.async>0.5] <- 1

prop.async.has.sd = prop.async[apply(prop.async,1,sd)!=0,]

sum(row.names(prop.async.has.sd) %in% async.list$ensbl)
#prop.async.has.sd.chr = async.list[which(rownames(prop.async.has.sd) %in% async.list$ensbl), 1]
prop.async.has.sd.chr = ensdb_genes[names(ensdb_genes.coding)%in% rownames(prop.async.has.sd)] 
prop.async.has.sd.chr.df = data.frame(prop.async.has.sd.chr)
prop.async.has.sd.index = which(names(ensdb_genes.coding)%in% rownames(prop.async.has.sd))

sum(apply(prop.async,2,sd)!=0)

c = cor((prop.async.has.sd), method = "spearman")
c.genes =  cor(t(prop.async.has.sd), method = "spearman")
write.csv(c, "Documents/NPC.clones.cor_cells.csv")
write.csv(c.genes, "Documents/NPC.clones.cor_genes.csv")
write.csv(prop.async.has.sd, "Documents/NPC.clones.expression.csv")


# Analysis of the data with corr to the early ASRT
keep = apply(table.async.ratio.naOmit,1, function(x) sum(x>0.9|x<0.1)>=6) 
table.async.ratio.naOmit.filter = table.async.ratio.naOmit[keep,]

prop.async = table.async.ratio.naOmit.filter
prop.async[prop.async<=0.5] <- -1
prop.async[prop.async>0.5] <- 1

prop.async.has.sd = prop.async[apply(prop.async,1,sd)!=0,]
match = match(row.names(prop.async.has.sd), async.list$ensbl)
prop.async.has.sd.early = async.list[match,"V4"] 

tmp = cbind.data.frame(prop.async.has.sd[,5], prop.async.has.sd.early)
chisq.test(table(tmp))


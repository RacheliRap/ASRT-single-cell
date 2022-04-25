############################################################################################
# clonal RNA-seq analysis.
# The data retrived from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148348
############################################################################################


library(gprofiler2)


setwd("C:/Users/rachelRa/")

file_names <- dir("Documents/", full.names = T, pattern = "GSM44*") #where you have your files
file.names.NPC = file_names[grep("_Ch", file_names)]

table <- do.call(cbind,lapply(file.names.NPC,read.delim))
table$Gene = gsub("\\..*","", table$Gene)

async.list = read.delim("Documents/async_master_list.protein_coding.bed", header = F)
gConvert = gconvert(query = async.list$V5, organism = "mmusculus",
         target="ENSG", mthreshold = Inf, filter_na = F)
async.list$ensbl = gConvert[match(async.list$V5, gConvert$input),"target"]

table.async = table[async.list$ensbl %in% table$Gene,]

table.async.ratio = table.async[,grep("allelicRatio_CAST_EiJ", names(table.async))]
row.names(table.async.ratio) <- table.async$Gene
table.async.ratio.naOmit = table.async.ratio[complete.cases(table.async.ratio),]
keep = apply(table.async.ratio.naOmit,1, function(x) all(x>0.7|x<0.3)) 
table.async.ratio.naOmit.filter = table.async.ratio.naOmit[keep,]

prop.async = table.async.ratio.naOmit.filter
prop.async[prop.async<0.3] <- -1
prop.async[prop.async>0.7] <- 1

prop.async.has.sd = prop.async[apply(prop.async,1,sd)!=0,]

sum(apply(prop.async,2,sd)!=0)

c = cor((prop.async.has.sd), method = "spearman")
c.genes =  cor(t(prop.async.has.sd), method = "spearman")
write.csv(c, "Documents/NPC.clones.cor_cells.csv")
write.csv(c.genes, "Documents/NPC.clones.cor_genes.csv")
write.csv(prop.async.has.sd, "Documents/NPC.clones.expression.csv")


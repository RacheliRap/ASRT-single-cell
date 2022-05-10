##################
# This scripts analysis scRepli-Seq done by Gilbert (accession GSE108556)
# The initial analysis was dowwnloading the data and excluded only Day7Diff_ESC_single_MidS
# The pre-processing code:
# for f in *cba_400k_100S_MAP_HMM2_eps0.01.binary.bedGraph; 
# do 
# msm=${f/cba/msm}; 
# name=`basename -s cba_400k_100S_MAP_HMM2_eps0.01.binary.bedGraph $f`;  
# paste -d "\t" $f $msm | awk '$4==(-1) && $8==(1) {print "msm"}; $4==(1) && $8==(-1) {print "cba"}; $4==1 && $8==1 {print "early"}; $4==(-1) && $8==(-1) {print "late"}; $4==0 || $8 == 0 {print "NA"}' > $name.asrt.bedgraph; done
# sed -i '1i '$name'' $name.asrt.bedgraph 
# done
#
###################

library(dplyr)
library(GenomicRanges)

table = read.delim("Documents/all.asrt.bedgraph")[,-1]
cord = read.delim("Documents/GSM2905041_P293_17_1_cba_400k_100S_MAP_HMM2_eps0.01.binary.bedGraph", header = F)
names(cord) <- c("seqnames", "start", "end", "early")
cord.gr <- makeGRangesFromDataFrame(cord)

async.list = read.delim("Documents/async_master_list.bed", header = T)
names(async.list) <- c("seqnames", "start", "end", "early")
async.list.gr = makeGRangesFromDataFrame(async.list)

table.async = table[cord.gr %over% async.list.gr,]
cord.async = cord[cord.gr %over% async.list.gr,]

table.noNA = table[apply(table, 1, function(x) sum(!is.na(x))>38),]
cord.noNA = cord[apply(table, 1, function(x) sum(!is.na(x))>38),]

table.noUniqueRows = table.noNA[apply(table.noNA, 1, function(x) length(unique(unlist(x)))!=1),]
table.asrt = table.noNA[apply(table.noNA, 1, function(x) sum(c("msm", "cba") %in% x)==2 & sum(x=="cba" | x=="msm")>=25) ,]

table.asrt.code = table.asrt
table.asrt.code[table.asrt.code == "msm"] <- 1
table.asrt.code[table.asrt.code == "cba"] <- (-1)
table.asrt.code[table.asrt.code == "late" | table.asrt.code == "early" | is.na(table.asrt.code)] <- 0

table.asrt.code[] <- lapply(table.asrt.code, as.numeric)

table.asrt.code = table.asrt.code[apply(table.asrt.code, 1, function(x) sd(x)!=0),]

corr = cor((table.asrt.code), method = "spearman")


table.asrt.code.async = table.asrt.code[row.names(table.asrt.code) %in% row.names(cord.async),]

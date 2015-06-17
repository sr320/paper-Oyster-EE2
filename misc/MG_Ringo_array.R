source("http://bioconductor.org/biocLite.R")
library(Ringo)
library(ChIPpeakAnno)
library(rtracklayer)

#1st step from Ryan's code is to generate annotation object --what is an annotation object?:
#this command takes a Nimblegen POS file a tab delim which 
#maps each probe to it's position in the reference genome
#we received these files on a CD with the slides - I've saved them here -> http://eagle.fish.washington.edu/bivalvia/index.php?dir=array%2Fnimblegen%2F
#ok, got it - going to run Ryan's code, but save into my own folder (added the verbose=F so it doesn't write out the whole file in the console)
#...note: once this file is made once it can always be reference going forward, if you want to use the file Ryan made and sent, it's here -> http://eagle.fish.washington.edu/bivalvia/index.php?dir=array%2Fnimblegen%2F
#probeAnno<-posToProbeAnno("121101_OysterV9_MG_meth.pos", genome="Pacific Oyster", microarrayPlatform="NimbleGen 121101 Oyster V9 MG meth",verbose = FALSE)
#save(probeAnno,file="121101_OysterV9_MG_meth.ProbeAnnotation_MG.rda")

#load annotation object - the file Ryan made
setwd("/Volumes/web/bivalvia/array/Rstuff")
load("121101_OysterV9_MG_meth.ProbeAnnotation.rda")

#load gff
gff<-read.delim("oyster.v9.glean.final.rename.mRNA.gff", header=F)
head (gff)
names(gff)<-c("chr","source","type","start","end","score","strand","phase","name")

#Ryan's annotation of this step: 'read in altered targets file such that the control samlpe is in the denominator slot (Cy3)'
#I'm not sure of the necessity of the date thing...
todaysDate<-substr(format(Sys.time(), "%Y.%m.%d."),1,11)
#The help for readNimblegen function says it's a function to read in Ng intensity text files into an RGlist
#Ryan provided a 'spottypes' file and 'targets' file in a folder called 'pair' , the intensity files also need to be in here (see path argument) - 
#spot types are made or read by limma, I think targets is too.  I need to ask him about this a little more..
#I took the pair files from the raw data from Cassie in the 'Oyster Array' folder, and put them in the 'pair' folder with the spot files and target types
#renamed them to be consistent with Ryan's naming convention in the Target file
RG<-readNimblegen("2013.09.24.targets.txt", "spottypes.txt", headerPattern="# software=DEVA", path="pair")
save(RG, file="2014.03.03.RG.rda")
#load a dataset written with the function save
load("2014.03.03.RG.rda")
#This step removes control probes (labeled as RANDOM (versus BLOCK1))
RG<-RG[RG$gene$Status != "RANDOM",]

#Preprocessing - taking raw probe intensities into (background corrected) normalized log ratios (M-values), using the loess normalization (note the first round of analysis used the 'nimblegen' method - should see why they changed but I do know there are fewer sites with the loess)
NG<-preprocess(RG,method="loess")
#renames the samples to be "targetIDnumerator_vs_targetIDdemoninator"
sampleNames(NG) <-with(RG$targets, paste(FileNameCy5,"vs",FileNameCy3,sep="_"))

#Ryan's annotation: ProbeIDs all have ":+:" in their names, which has been converted to "..." when the preprocess function is run.  This is in conflict
#with the probeAnno
row.names(exprs(NG))<-gsub("...",":+:",row.names(exprs(NG)), fixed=TRUE)

#smooth probe intensities
#computes running medians on tiling expression sets, arguements include the normalized file, the probe anno made earlier, the column that has the 'phenoData' 
#win half size says this: Half the size of the window centered at a probe position, in which all other probes contribute to the calculation of the median.need to ask about this one
#in Ryan's original code the winHalfSize was 500 (full size 1000).  When they did the reananlysis the window size was changed to 600, so I changed this to 300
smoothNG<-computeRunningMedians(NG,probeAnno=probeAnno, modColumn="FileNameCy5", winHalfSize=300,verbose = FALSE)


###########Not running this next section, but I'm going to leave it in since this is how the original analysis was done.  Ryan and Jeff then suggested more conservative threshold:
###########thresholds for the final analysis: 2 fold change (log2(2) =1) for A02 and A03, and a less stringent threshold of 1.4 fold change in A01 (log2(1.4)=0.485)
###########when I asked why different threshold this was response: A lower threshold was used for the input vs input (A01) sample in an attempt to identify the most robust DMRs.  When using the same threshold for all samples, there were a large number of DMRs in A02 and A03, where the same, but slightly weaker pattern of logFC values was found.  
#the apply means applying functions over error margine, 2 indicates columns, 
#upperBoundNull is a function from Ringo - This function tries to pinpoint the minimum of data values (log2 ratios) that are more likely to arise from the alternative distribution, i.e. an upper bound for values following the null distribution.
#in the final analysis these thresholds weren't used.  Instead, a more conservative pick of Use a threshold of 2 fold change (log2(2) =1) for A02 and A03, and a less stringent threshold of 1.4 fold change in A01 (log2(1.4)=0.485), 
#y0 <- apply(exprs(smoothNG),2,upperBoundNull)
# A01_EE2.input_532.pair_vs_A01_Ctrl.input_635.pair.sm     A02_EE2.MBD_532.pair_vs_A02_Ctrl.MBD_635.pair.sm 
# 0.6742440                                            0.6978518 
# A03_EE2.MBD_635.pair_vs_A03_Ctrl.MBD_532.pair.sm 
# 0.7904904 

#final threshold selection
y0 <- c(.485,1,1)

#'chers' are ChIP enriched regions, thisfirst one uses the smoothed probe intensities, the enviroment of mapping 'probe anno' and the thresholds for
#for which you call something enriched, in this case it'sthe vcetor y0, which is the value of the upper bound of the null distribution - values above this are enriched, distCutOff is the
#max #of bp at which enriched probes are condensed into 1 cher
chersNG <- findChersOnSmoothed(smoothNG, probeAnno=probeAnno, thresholds=c(y0[1],y0[2],y0[3]), distCutOff=600, verbose = FALSE)

#this function relates found chers to annotated genomic features, such as transcripts
#gff is the mRNA regions (see above when it was read in)
chersNG <- relateChers(chersNG, gff)
#weird, this says called internally by other ringo functions, not normalled called by user
chersNGD <- as.data.frame.cherList(chersNG)
#order them in decreasing order by max level
chersNGD<-chersNGD[order(chersNGD$maxLevel, decreasing=TRUE),]
nrow(chersNGD) # --number of chers #4012

#chers detected within each array
nrow(chersNGD[grep("A01",chersNGD$name),])   #3036
nrow(chersNGD[grep("A02",chersNGD$name),]) #448
nrow(chersNGD[grep("A03",chersNGD$name),]) #528

###note: I will come back to this later if I have time, but skipping it for now
##Ryan's annotation: try to write wig using the writeWIG function from Starr
#library(Starr)
#writeWIG(smoothNG[,1],probeAnno,file="2013.09.25.A01_smoothed.wig")
#writeWIG(smoothNG[,2],probeAnno,file="2013.09.25.A02_smoothed.wig")
#writeWIG(smoothNG[,3],probeAnno,file="2013.09.25.A03_smoothed.wig")
##had to change NAs to 0s in order to load into IGV

#find overlapping regions
#create ranged data objects
#looks like here it's repeating finding chers, but each assay (e.g. A01) is generated separately
chersA01 <- findChersOnSmoothed(smoothNG[,1], probeAnno=probeAnno, thresholds=c(y0[1]), distCutOff=600, verbose = FALSE)
#described above
chersA01 <- relateChers(chersA01, gff)
#exports cherList into a file of gff or bed
exportCherList(chersA01,filename="2014.03.03.chers_A01_HYPER.bed",format="bed",genome="oyster")
#read back in without header
bedA01<-read.delim("2014.03.03.chers_A01_HYPER.bed", header=F)
#generate a 6th column with '1' in every row -it's a default for strand column
bedA01$V6<-1
#this ranged data is some version of a reformatted bed
#it's used later for findOverlapping peaks (I wonder why bedtools wouldn't be simpler for this?)
#this is required for the later step
rangedA01<-BED2RangedData(bedA01)

#doing the same thing for assay A02 and A03
chersA02 <- findChersOnSmoothed(smoothNG[,2], probeAnno=probeAnno, thresholds=c(y0[2]), distCutOff=600, verbose = FALSE)
chersA02<- relateChers(chersA02, gff)
exportCherList(chersA02,filename="2014.03.03.chers_A02_HYPER.bed",format="bed",genome="oyster")
bedA02<-read.delim("2014.03.03.chers_A02_HYPER.bed", header=F)
bedA02$V6<-1
rangedA02<-BED2RangedData(bedA02)

chersA03 <- findChersOnSmoothed(smoothNG[,3], probeAnno=probeAnno, thresholds=c(y0[3]), distCutOff=600, verbose = FALSE)
chersA03<- relateChers(chersA03, gff)
exportCherList(chersA03,filename="2014.03.03.chers_A03_HYPER.bed",format="bed",genome="oyster")
bedA03<-read.delim("2014.03.03.chers_A03_HYPER.bed", header=F)
bedA03$V6<-1
rangedA03<-BED2RangedData(bedA03)

#generate objects that contain overlapping regions
#this requires the 'ranged' data generated above
overlapA01.A02<-findOverlappingPeaks(rangedA01,rangedA02)
overlapA01.A03<-findOverlappingPeaks(rangedA01,rangedA03)
overlapA02.A03<-findOverlappingPeaks(rangedA02,rangedA03)

#number of overlapping peaks in each pairwise comparison
nrow(overlapA01.A02$OverlappingPeaks) #397
nrow(overlapA01.A03$OverlappingPeaks) #462
nrow(overlapA02.A03$OverlappingPeaks) #353

#This shows how multiple peaks in A02 are covered by 1 peak in A03. Here 346 regions of A02 overlap with A03, but 344 regions of A03 overlap with A02.
#############
nrow(rangedA02) #448
A02peaksWithOverlapA03<-overlapA02.A03$Peaks1withOverlaps; nrow(A02peaksWithOverlapA03) #346
A02noOverlapWithA03<-rangedA02[!rownames(rangedA02) %in% rownames(A02peaksWithOverlapA03),]; nrow(A02noOverlapWithA03) #102

nrow(rangedA03) #528
A03peaksWithOverlapA02<-overlapA02.A03$Peaks2withOverlaps; nrow(A03peaksWithOverlapA02) #344
A03noOverlapWithA02<-rangedA03[!rownames(rangedA03) %in% rownames(A03peaksWithOverlapA02),]; nrow(A03noOverlapWithA02) #184
################

#Ryan's annotation: we now have overlapping and nonoverlapping chers between A02 and A03, what we want now, is to remove those chers that also overlap with A01 from each set
A02peaksWithOverlapA01<-overlapA01.A02$Peaks2withOverlaps; nrow(A02peaksWithOverlapA01) #394
A02peaksWithOverlapA03butNotA01<-A02peaksWithOverlapA03[!rownames(A02peaksWithOverlapA03) %in% rownames(A02peaksWithOverlapA01),]; nrow(A02peaksWithOverlapA03butNotA01) #30

A03peaksWithOverlapA01<-overlapA01.A03$Peaks2withOverlaps; nrow(A03peaksWithOverlapA01) #461
A03peaksWithOverlapA02butNotA01<-A03peaksWithOverlapA02[!rownames(A03peaksWithOverlapA02) %in% rownames(A03peaksWithOverlapA01),]; nrow(A03peaksWithOverlapA02butNotA01) #28

#write RangedData to bed
#Ryan's annotation: https://stat.ethz.ch/pipermail/bioconductor/2012-April/045050.html
write.table(as.data.frame(A02peaksWithOverlapA03butNotA01)[,c(1,2,3)], file="2014.03.03.A02overlapA03_butNotA01_HYPER.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(as.data.frame(A03peaksWithOverlapA02butNotA01)[,c(1,2,3)], file="2014.03.03.A03overlapA02_butNotA01_HYPER.bed",sep="\t",quote=F,row.names=F,col.names=F)

######################################################################################################################
# Hypomethylation

#flip the ratios on these such that we're looking for hypomethylated regions
exprs(smoothNG) = -exprs(smoothNG)


chersNG <- findChersOnSmoothed(smoothNG, probeAnno=probeAnno, thresholds=c(y0[1],y0[2],y0[3]), distCutOff=600, verbose = FALSE)


chersNG <- relateChers(chersNG, gff)
chersNGD <- as.data.frame.cherList(chersNG)
chersNGD<-chersNGD[order(chersNGD$maxLevel, decreasing=TRUE),]
nrow(chersNGD) # --number of chers #3518

#chers detected within each array
nrow(chersNGD[grep("A01",chersNGD$name),])   #2758
nrow(chersNGD[grep("A02",chersNGD$name),]) #348
nrow(chersNGD[grep("A03",chersNGD$name),]) #412


#repeating finding chers, but each assay (e.g. A01) is generated separately
chersA01 <- findChersOnSmoothed(smoothNG[,1], probeAnno=probeAnno, thresholds=c(y0[1]), distCutOff=600, verbose = FALSE)
chersA01 <- relateChers(chersA01, gff)
exportCherList(chersA01,filename="2014.03.03.chers_A01_HYPO.bed",format="bed",genome="oyster")
bedA01<-read.delim("2014.03.03.chers_A01_HYPO.bed", header=F)
bedA01$V6<-1
rangedA01<-BED2RangedData(bedA01)

#doing the same thing for assay A02 and A03
chersA02 <- findChersOnSmoothed(smoothNG[,2], probeAnno=probeAnno, thresholds=c(y0[2]), distCutOff=600, verbose = FALSE)
chersA02<- relateChers(chersA02, gff)
exportCherList(chersA02,filename="2014.03.03.chers_A02_HYPO.bed",format="bed",genome="oyster")
bedA02<-read.delim("2014.03.03.chers_A02_HYPO.bed", header=F)
bedA02$V6<-1
rangedA02<-BED2RangedData(bedA02)

chersA03 <- findChersOnSmoothed(smoothNG[,3], probeAnno=probeAnno, thresholds=c(y0[3]), distCutOff=600, verbose = FALSE)
chersA03<- relateChers(chersA03, gff)
exportCherList(chersA03,filename="2014.03.03.chers_A03_HYPO.bed",format="bed",genome="oyster")
bedA03<-read.delim("2014.03.03.chers_A03_HYPO.bed", header=F)
bedA03$V6<-1
rangedA03<-BED2RangedData(bedA03)

#generate objects that contain overlapping regions
#this requires the 'ranged' data generated above
overlapA01.A02<-findOverlappingPeaks(rangedA01,rangedA02)
overlapA01.A03<-findOverlappingPeaks(rangedA01,rangedA03)
overlapA02.A03<-findOverlappingPeaks(rangedA02,rangedA03)

#number of overlapping peaks in each pairwise comparison
nrow(overlapA01.A02$OverlappingPeaks) #317
nrow(overlapA01.A03$OverlappingPeaks) #369
nrow(overlapA02.A03$OverlappingPeaks) #269

#This shows how multiple peaks in A03 are covered by 1 peak in A02. Here 255 regions of A02 overlap with A03, but 265 regions of A03 overlap with A02.
#############
nrow(rangedA02) #348
A02peaksWithOverlapA03<-overlapA02.A03$Peaks1withOverlaps; nrow(A02peaksWithOverlapA03) #255
A02noOverlapWithA03<-rangedA02[!rownames(rangedA02) %in% rownames(A02peaksWithOverlapA03),]; nrow(A02noOverlapWithA03) #93

nrow(rangedA03) #412
A03peaksWithOverlapA02<-overlapA02.A03$Peaks2withOverlaps; nrow(A03peaksWithOverlapA02) #265
A03noOverlapWithA02<-rangedA03[!rownames(rangedA03) %in% rownames(A03peaksWithOverlapA02),]; nrow(A03noOverlapWithA02) #147
################

#Ryan's annotation: we now have overlapping and nonoverlapping chers between A02 and A03, what we want now, is to remove those chers that also overlap with A01 from each set
A02peaksWithOverlapA01<-overlapA01.A02$Peaks2withOverlaps; nrow(A02peaksWithOverlapA01) #315
A02peaksWithOverlapA03butNotA01<-A02peaksWithOverlapA03[!rownames(A02peaksWithOverlapA03) %in% rownames(A02peaksWithOverlapA01),]; nrow(A02peaksWithOverlapA03butNotA01) #19

A03peaksWithOverlapA01<-overlapA01.A03$Peaks2withOverlaps; nrow(A03peaksWithOverlapA01) #364
A03peaksWithOverlapA02butNotA01<-A03peaksWithOverlapA02[!rownames(A03peaksWithOverlapA02) %in% rownames(A03peaksWithOverlapA01),]; nrow(A03peaksWithOverlapA02butNotA01) #18

#write RangedData to bed
#Ryan's annotation: https://stat.ethz.ch/pipermail/bioconductor/2012-April/045050.html
write.table(as.data.frame(A02peaksWithOverlapA03butNotA01)[,c(1,2,3)], file="2014.03.03.A02overlapA03_butNotA01_HYPO.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(as.data.frame(A03peaksWithOverlapA02butNotA01)[,c(1,2,3)], file="2014.03.03.A03overlapA02_butNotA01_HYPO.bed",sep="\t",quote=F,row.names=F,col.names=F)


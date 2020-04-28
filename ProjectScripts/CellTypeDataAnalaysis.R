source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
plotMA = DESeq2::plotMA
packageF("org.Hs.eg.db")
packageF("tabulizer")
packageF("limma")
packageF("edgeR")
packageF("Rsubread")

if("CellTypeSpecificCounts.Rds" %in% list.files(path = "Data")){
  CountsCells <- readRDS("Data/CellTypeSpecificCounts.Rds")
} else {
  #Import the peakset from Marzi et al. and create a SAF file
  temp <- read.table("Data/GSE102538_H3K27ac_EntorhinalCortex.bed.gz", header = F, sep = " ") %>%
    select(V4, V1, V2, V3)
  
  names(temp) <- c("GeneID" , "Chr",  "Start", "End")
  temp$Strand <- "*"
  
  #Quantified the H3K27ac reads from the NeuN+/- cellls based on Girdahal et al 2018 inside Marzi et al peaks
  temp2 <- featureCounts(files = list.files(path = "Data/CellTypesH3K27ac/BAMfiles/", full.names = T),nthreads = 5,
                         annot.inbuilt = "hg19",annot.ext = temp, isPairedEnd = T, countMultiMappingReads = F, countChimericFragments=F)
  
  SampleNames <- list.files(path = "Data/CellTypesH3K27ac/BAMfiles/", full.names = F)
  SampleNames <- sapply(SampleNames, function(x){
    x <- paste0(strsplit(x, "_")[[1]][c(7, 11)], collapse = "_")
    x <- gsub(".bam", "", x)
    if(grepl("\\+", x)){
      gsub("\\+", "pos", x)
    } else {
      gsub("-", "neg", x)
    }
  })
  
  CountsCells <- temp2$counts
  colnames(CountsCells) <- SampleNames
  saveRDS(CountsCells, file = "Data/CellTypeSpecificCounts.Rds")
}
  
CellTypeMeta <- read.table("Data/CellTypesH3K27ac/MSSM_U01MH103392_EpiMap_Metadata_ChIPseq_August2016Release.csv",
                           header = T, sep = ",", comment.char = "!") %>% filter(HistoneMark == "H3K27ac", BrainRegion == "DLPFC")
CellTypeMeta2 <- read.table("Data/CellTypesH3K27ac/MSSM_U01MH103392_EpiMap_Metadata_clinical.csv",
                            header = T, sep = ",") 

CellTypeMeta <- merge(CellTypeMeta, CellTypeMeta2, by = "Individual_ID")

CellTypeMeta$Sample_ID2 <- sapply(as.character(CellTypeMeta$File_Name), function(x){
  x <- paste0(strsplit(x, "_")[[1]][c(7, 11)], collapse = "_")
  if(grepl("\\+", x)){
    gsub("\\+", "pos", x)
  } else {
    gsub("-", "neg", x)
  }
})

CellTypeMeta$CellType2 <- sapply(as.character(CellTypeMeta$CellType), function(x){
  if(grepl("\\+", x)){
    "Neuron"
  } else {
    "Glia"
  }
})

CellTypeMeta %<>% mutate(Hemisphere = tolower(Hemisphere))

#Run edgeR for the difference between the cell types, similar to analysis by Marzi et al.
#Keep peaks with cpm > 1 in at least two samples
keep <- rowSums(cpm(CountsCells) > 1) >= 2
CountsCells <- CountsCells[keep,]

#Arrance the order of samples in the metadata file
CellTypeMeta <- CellTypeMeta[match(colnames(CountsCells), CellTypeMeta$Sample_ID2),]

CellTypecountList<-DGEList(counts=CountsCells, group = CellTypeMeta$CellType2)
CellTypecountList <- calcNormFactors(CellTypecountList)


CellTypeModelMatrix <- model.matrix(as.formula("~Sex + Hemisphere + AgeDeath + pH + CellType2"), data = CellTypeMeta)

CellTypecountList <- estimateDisp(CellTypecountList, CellTypeModelMatrix)

fitTMM <- glmQLFit(CellTypecountList, CellTypeModelMatrix)
qlf_CellType <- glmQLFTest(fitTMM, coef = "CellType2Neuron")
CellType_results <- topTags(qlf_CellType, n = Inf) %>% data.frame() %>% mutate(PeakName = rownames(.))

CellTypeResult_Anno <- CellType_results %>%
  AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot %>% select(-matches("GSM")), by.x = "PeakName", by.y = "PeakName") %>% arrange(FDR)

CellType_resultsSignifMarziCETS <- CellType_results %>% filter(PeakName %in% SignifCETS$PeakName)
CellType_resultsSignifMarziCETS$MethodSignif = "CETs"
CellType_resultsSignifMarziCETS$DirectictionChange <- sapply(CellType_resultsSignifMarziCETS$PeakName, function(x){
  if(SignifCETS %>% filter(PeakName == x) %>% .$logFC < 0){
    "Hypoacetylated"
  } else if(SignifCETS %>% filter(PeakName == x) %>% .$logFC > 0){
    "Hyperacetylated"
  }
})

CellType_resultsSignifMarziMSP <- CellType_results %>% filter(PeakName %in% SignifMSP$PeakName)
CellType_resultsSignifMarziMSP$MethodSignif = "MSP"
CellType_resultsSignifMarziMSP$DirectictionChange <- sapply(CellType_resultsSignifMarziMSP$PeakName, function(x){
  if(SignifMSP %>% filter(PeakName == x) %>% .$logFC < 0){
    "Hypoacetylated"
  } else if(SignifMSP %>% filter(PeakName == x) %>% .$logFC > 0){
    "Hyperacetylated"
  }
})

CellType_resultsSignifMarzi <- rbind(CellType_resultsSignifMarziCETS, CellType_resultsSignifMarziMSP)
CellType_resultsSignifMarzi$DirectictionChange <- factor(CellType_resultsSignifMarzi$DirectictionChange, levels = c("Hypoacetylated, hyperacetylated"))

ggplot(CellType_resultsSignifMarzi, aes(DirectictionChange, logFC, color = DirectictionChange)) +
  theme_minimal() +
  labs(x = "", y = "logFC (Neuron/Glia)") +
  geom_violin(aes(fill = DirectictionChange)) +
  geom_boxplot(width = 0.2)+
  scale_color_manual(values =  c("orange", "chartreuse4")) +
  scale_fill_manual(values =  c("orange", "chartreuse4")) +
  geom_hline(yintercept = 0, color = "red") +
  facet_wrap(~MethodSignif)

ggsave(paste0("MethodComparisonBoxplot.tif"), device = "tiff", width = 8, height = 6, dpi = 300, path = ResultsPath)


CellType_Anno <- CellType_results %>%
  AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName") %>% arrange(FDR)


#Get the TMM psudo counts

CountCellsTMM <- cpm(CellTypecountList)
CountCellsTMM <- apply(CountCellsTMM, c(1,2), function(x) log2(x+1))

CountCellsTMM2 <- cpm(CellTypecountList, log = T)


CountCellsTMMDown <- CountCellsTMM[rownames(CountCellsTMM) %in% c(SignifCETS %>% filter(logFC < 0) %>% .$PeakName),]

#Get the significantly hyperacetylated peaks in Marzi et al
CountCellsTMMUp <- CountCellsTMM[rownames(CountCellsTMM) %in% c(SignifCETS %>% filter(logFC > 0) %>% .$PeakName),]

#plot
packageF("pheatmap")

#Plot randmly selected 40K peaks
pheatmap(CountCellsTMM[sample(1:nrow(CountCellsTMM), 40000, replace = F),],
         scale = "row", border_color = NA, show_rownames = F, main = "Random peaks (40K)",
         filename = paste0(ResultsPath, "HeatmapCellTMMrandomPeak.pdf"),  width = 8, height = 8)

pheatmap(CountCellsTMMDown, scale = "row", border_color = NA, show_rownames = F, main = "CETs, Hypoacetylated",
         filename = paste0(ResultsPath, "HeatmapCellTMM_CETs_Hypoacetylated.pdf"),  width = 8, height = 8)
pheatmap(CountCellsTMMUp, scale = "row", border_color = NA, show_rownames = F, main = "CETs, Hyperacetylated",
         filename = paste0(ResultsPath, "HeatmapCellTMM_CETs_Hyperacetylated.pdf"),  width = 8, height = 8)


#Get the significantly hypoacetylated peaks in after MSP normalization
CountCellsTMMDownMSP <- CountCellsTMM[rownames(CountCellsTMM) %in% c(SignifMSP %>% filter(logFC < 0) %>% .$PeakName),]

#Get the significantly hyperacetylated peaks after MSP normalization
CountCellsTMMUpMSP <- CountCellsTMM[rownames(CountCellsTMM) %in% c(SignifMSP %>% filter(logFC > 0) %>% .$PeakName),]

pheatmap(CountCellsTMMDownMSP, scale = "row", border_color = NA, show_rownames = F, main = "MSP, Hypoacetylated",
         filename = paste0(ResultsPath, "HeatmapCellTMM_MSPs_Hypoacetylated.pdf"),  width = 8, height = 8)
pheatmap(CountCellsTMMUpMSP, scale = "row", border_color = NA, show_rownames = F, main = "MSP, Hyperacetylated",
         filename = paste0(ResultsPath, "HeatmapCellTMM_MSPs_Hyperacetylated.pdf"),  width = 8, height = 8)


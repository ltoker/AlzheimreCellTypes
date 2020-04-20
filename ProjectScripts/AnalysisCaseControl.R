source("generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
plotMA = DESeq2::plotMA
packageF("org.Hs.eg.db")
packageF("tabulizer")
packageF("limma")
packageF("edgeR")

ResultsPath = "Results"
if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")


CellTypePeakCountLoc = "Data/marzi_counts.tsv"
CountMatrixLoc = "Data/MarziPaperCounts.tsv"
Cohort = "Marzi"

MetadataSup <- extract_tables("41593_2018_253_MOESM1_ESM.pdf", pages = c(27:30)) %>% lapply(data.frame) %>% rbindlist() #Reading the metadata from the supplementary table 
names(MetadataSup) <- c("SampleID", "Group", "BraakStage", "Age", "Sex", "NeuralProportion", "PMI", "Experiment")
MetadataSup <- MetadataSup[-c(1:3),]
MetadataSup$AgeNumeric <- as.numeric(as.character(MetadataSup$Age))
MetadataSup$PMINumeric <- as.numeric(as.character(MetadataSup$PMI))
MetadataSup$NeuralProportionNumeric <- sapply(MetadataSup$NeuralProportion, function(x){
  if(grepl("%", x)){
    x = gsub("%", "", x) %>% as.character() %>% as.numeric()
    x/100
  } else {
    NA
  }
})

#Getting additional metadata
softDown("GSE102538", file = "Data/GSE102538.soft")
Metadata <- ReadSoft("Data/GSE102538.soft") %>% select(-matches("antib|Sample_source|Sample_platform|orga")) %>% data.frame()
names(Metadata) <- c("GSM", "SampleName", "Group", "Age", "Sex", "CETS")
Metadata$Age <- as.numeric(as.character(Metadata$Age))
Metadata$CETS <- as.numeric(as.character(Metadata$CETS))

Metadata %<>% arrange(CETS)
Metadata$Eno2 <- c(0.65, 0.99, 0.56, 0.48, 0.56, 0.55, 0.69, 0.26, 1.27, 0.88,
                      0.42, 0.51, 2.07, 0.84, 0.52, 1.30, 0.90, 0.80, 0.19, 1.16,
                      0.34, 0.69, 0.34, 0.77, 1.43, 0.72, 0.33, 0.36, 0.29, 0.59,
                      0.42, 0.83, 0.49, 1.53, 0.98, 1.31, 0.85, 1.61, 0.65, 1.11,
                      0.72, 1.26, 0.80, 0.48,  2.28, 1.17, NA)
Metadata$deltaCTEno2 <- sapply(Metadata$Eno2, function(x){
  if(!is.na(x)){
    -log2(x) 
  } else {
    NA
  }
}) 

#Following the steps in Marzi et al.
Metadata$CETS[is.na(Metadata$CETS)] <- median(Metadata$CETS, na.rm = T)
Metadata$CETSif <-  cut(Metadata$CETS, breaks=5, ordered_result = T)
Metadata$Agef <-  cut(Metadata$Age, breaks=5, ordered_result = T)


Metadata %<>% arrange(GSM)
rownames(Metadata) <- Metadata$GSM %>% as.character()
Metadata$Group <- factor(Metadata$Group, levels = c("C", "AD"))
##### Add relative cell proportion for based on differential NeuN positive and negative cell H3K27ac peaks ##########  
Metadata <- GetCellularProportions(Metadata = Metadata)

MetaCellMelt <- gather(Metadata %>% select(matches("Group|CETS$|^Eno2|MSP"), -OligoPrecursors_MSP, -Microglia_deactivation_MSP, -Microglia_activation_MSP), key = "CountType", value = "Value", -Group)
MetaCellMelt$CountType <- factor(MetaCellMelt$CountType, levels = c("Astrocyte_MSP" , "Endothelial_MSP", "Microglia_MSP", "Microglia_activation_MSP",
                                                          "Oligo_MSP", "GabaVIPReln_MSP", "Pyramidal_MSP", "NeuNall_MSP", "CETS","Eno2"))

ggplot(MetaCellMelt, aes(Group, Value, color = Group)) +
  theme_minimal() +
  geom_boxplot() +
  geom_jitter(height = 0) +
  facet_wrap(~CountType, scales = "free_y", nrow = 3)



############################# HTseq counts ######################################################################

annoFileCollapsed <- GetGenomeAnno(genome = "hg19")

HTseqCounts <- read.table(CountMatrixLoc, header = T, sep = "\t")

names(HTseqCounts) <- sapply(names(HTseqCounts), function(x){
  x = gsub(".*bamfiles.0?", "", x)
  strsplit(x,"_")[[1]][1]
})


HTseqCounts %<>% mutate(Peak.Location = paste0("chr", Chr, ":", Start, "-", End))
names(HTseqCounts)[2:4] <- c("CHR", "START", "END")

#Arrange the samples to match Metadata order
HTseqCounts %<>% select(c("Geneid", "CHR", "START", "END", "Strand",  "Length", as.character(Metadata$GSM)))

#Keep peaks with cpm > 1 in at least two samples
keep <- rowSums(cpm(HTseqCounts %>% select(matches("GSM"))) > 1) >= 2
HTseqCounts <- HTseqCounts[keep,]

HTseqCounts$MaxCount = apply(HTseqCounts %>% select(matches("GSM")), 1, max)
HTseqCounts %<>% mutate(NormalizedMaxCount = MaxCount/Length)  


AllCalledData <- GetCountMatrixHTseq(countsDF = HTseqCounts, meta = Metadata, MetaCol = c("GSM", "SampleName", "Group", "Age", "Agef", "Sex", "CETS", "CETSif", grep("Eno2|MSP", names(Metadata), value = "T")),
                                     MetaSamleCol = "GSM", countSampleRegEx = "GSM", OtherNormRegEx = "^C1orf43_|^CHMP2A_|^EMC7_|^GPI_|^PSMB2_|^PSMB4_|^RAB7A_|^REEP5_|^SNRPD3_|^VCP_|^VPS29")

closeDev()

###########################################################################################################
############ RERUN USING PEAKS BASED ON ALL THE SAMPLES, without collapsing peaks #########################
###########################################################################################################
pdf(paste0(ResultsPath, "SampleCorAllPeaks", Cohort, ".pdf"), useDingbats = F, width = 10, height = 8)
countMatrixFullAllCalled <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)), collapseBy = "PeakName",MetaSamleCol = "GSM",
                                               CorMethod = "pearson", countSampleRegEx = "GSM",MetaSamleIDCol = "GSM", groupCol = "Group",
                                               FilterBy = "", meta = AllCalledData$SampleInfo, title = paste0("Sample correlation, ", Cohort))
closeDev()

#Get the pvalues for associasion of each covariate with the first 5 PCs
PCAsamples <- prcomp(t(countMatrixFullAllCalled$CPMdata), scale. = T)
countMatrixFullAllCalled$Metadata %<>% mutate(PC1 = PCAsamples$x[,1],
                                              PC2 = PCAsamples$x[,2],
                                              PC3 = PCAsamples$x[,3],
                                              PC4 = PCAsamples$x[,4],
                                              PC5 = PCAsamples$x[,5]) 
VarExplained <- PCAsamples %>% summary() %>% .$importance %>%
  .[2, 1:sum(grepl("^PC", names(countMatrixFullAllCalled$Metadata)))]*100 


CovarPvalues <- sapply(grep("^PC", names(countMatrixFullAllCalled$Metadata), value = T), function(PC){
  temp <- lm(as.formula(paste0(PC, "~ Group + Age + Sex + CETS + NeuNall_MSP + Astrocyte_MSP + Microglia_MSP + Oligo_MSP")),
             data = countMatrixFullAllCalled$Metadata) %>% summary
  temp$coefficients[-1,4]
}, simplify = F) %>% do.call(cbind, .) %>% data.frame()

names(CovarPvalues) <- paste0(names(CovarPvalues), "(", round(VarExplained, digits = 1), "%)")
CovarPvalues %<>% mutate(Variable = factor(rownames(CovarPvalues), levels = rownames(CovarPvalues)))
CovarPvaluesMelt <- gather(CovarPvalues, key = "PC", value = "pValue", -Variable)
CovarPvaluesMelt %<>% mutate(pPvalue = -log10(pValue)) 
CovarPvaluesMelt$Significant <- sapply(CovarPvaluesMelt$pValue, function(x){
  if(x < 0.05){
    "Yes"
  } else {
    "No"
  }
}) %>% factor(levels = c("No", "Yes"))

Plot  <- ggplot(CovarPvaluesMelt, aes(PC, Variable)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0)) +
  labs(x = "", y = "", title = "Association of variables with main PCs") +
  geom_tile(aes(fill = pPvalue), colour = "white") +
  scale_fill_gradient(high = "steelblue", low = "gray94", name = "-log10(p)") +
  geom_text(aes(label = signif(pValue, 2), colour = Significant), show.legend = F) +
  scale_color_manual(values = c("black", "red"))
ggsave(paste0("AssociationWithPCs", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 8, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()


#Look at the correlation betweem CETS and ChIP-seq based cell abundance estimates
CellData <- countMatrixFullAllCalled$Metadata %>% select(-matches("__|Total")) %>% gather(matches("MSP|CETS$"), key = "CellType", value = "MSP") 
CellData$CellType <- sapply(CellData$CellType, function(x) gsub("_MSP", "", x)) %>%
  factor(levels = c("Astrocyte", "Endothelial", "Microglia", "Microglia_activation",
                    "Microglia_deactivation","Oligo", "OligoPrecursors",
                    "GabaVIPReln", "Pyramidal", "NeuNall","CETS"))
CellData %<>% mutate(MSPnorm = MSP/MeanRatioOrg)


Plot <- ggplot(CellData, aes(Group, MSP, color = Group)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  geom_boxplot(outlier.shape = NA, aes(fill = Group), alpha = 0.4) +
  geom_jitter(width = 0.2, aes(color = Group)) +
  scale_color_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_fill_manual(values = c("dodgerblue4", "chocolate1"), name = "Group") +
  scale_shape_manual(values = c(16, 4)) +
  facet_wrap(~CellType, scales = "free")

ggsave(paste0("CellTypeDifferences", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 6, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()

CellTypeStats <- sapply(levels(CellData$CellType), function(cellType){
  Data = CellData %>% filter(CellType == cellType)
  lm(MSP~Group + Sex + Agef, data = Data) 
}, simplify = F)

CellTypeStatSummary <- lapply(CellTypeStats, function(x){
  temp <- summary(x) %>% .$coef %>% .[2, c(1,4)]
  temp2 <- confint(x) %>% .[2,]
  out <- c(temp, temp2) %>% t
  colnames(out) <- c("Coeficient", "pValue", "Low", "High")
  out
}) %>% do.call(rbind, .) %>% data.frame()
rownames(CellTypeStatSummary) <- names(CellTypeStats)
CellTypeStatSummary %<>% mutate(CellType = factor(rownames(.), levels = rownames(.)))


CellData2 <- countMatrixFullAllCalled$Metadata %>% select(-matches("__|Total")) %>% gather(matches("MSP|CETS$"), -NeuNall_MSP, key = "CellType", value = "MSP") 
CellData2$CellType <- sapply(CellData2$CellType, function(x) gsub("_MSP", "", x)) %>%
  factor(levels = c("Astrocyte", "Endothelial", "Microglia", "Microglia_activation",
                    "Microglia_deactivation","Oligo", "OligoPrecursors",
                    "GabaVIPReln", "Pyramidal","CETS"))
CellData2 %<>% mutate(MSPnorm = MSP/MeanRatioOrg)
CellTypeStats2 <- sapply(levels(CellData2$CellType), function(cellType){
  Data = CellData2 %>% filter(CellType == cellType)
  lm(MSP~Group + Sex + Agef + NeuNall_MSP, data = Data) 
}, simplify = F)

CellTypeStatSummary2 <- lapply(CellTypeStats2, function(x){
  temp <- summary(x) %>% .$coef %>% .[2, c(1,4)]
  temp2 <- confint(x) %>% .[2,]
  out <- c(temp, temp2) %>% t
  colnames(out) <- c("Coeficient", "pValue", "Low", "High")
  out
}) %>% do.call(rbind, .) %>% data.frame
rownames(CellTypeStatSummary2) <- names(CellTypeStats2)
CellTypeStatSummary2 %<>% mutate(CellType = factor(rownames(.), levels = rownames(.)))


ggplot(CellTypeStatSummary, aes(CellType, Coeficient)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", title = ) +
  geom_point() +
  geom_errorbar(aes(ymin = Low, ymax = High)) +
  geom_hline(yintercept = 0, color = "red", linetype = 2)


ggplot(CellTypeStatSummary2, aes(CellType, Coeficient)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", title = ) +
  geom_point() +
  geom_errorbar(aes(ymin = Low, ymax = High)) +
  geom_hline(yintercept = 0, color = "red", linetype = 2)


ggplot(CellData %>% filter(CellType == "NeuNall"), aes(Group, MSP, color = Group)) +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Neuronal MSP", title = "MSP analysis enthorinal cortex") +
  geom_boxplot(outlier.shape = NA, aes(fill = Group), alpha = 0.4) +
  geom_jitter(width = 0.2, aes(color = Group)) +
  scale_color_manual(values = c("darkgrey", "darkred"), name = "Group") +
  scale_fill_manual(values = c("darkgrey", "darkred"), name = "Group") +
  annotate("text", x = 1.5, y = 1, label = paste0("p = ", signif(CellTypeStats$NeuNall$coefficients[2,4], digits = 2)))

ggsave(paste0("NeuNDifferences", Cohort, ".pdf"), device = "pdf", width = 5, height = 5, dpi = 300, useDingbats = F, path = ResultsPath)
closeDev()

#Detect outliers
MedianCor <- apply(countMatrixFullAllCalled$SampleCor, 1, function(x) median(x, na.rm = TRUE))
Outlier <- MedianCor[MedianCor < (median(MedianCor) - 1.5*iqr(MedianCor))]


#Run analyses (repeating the steps in Marzi et al.)
ddsAD <- DESeqDataSetFromMatrix(countData = countMatrixFullAllCalled$countMatrix, colDat = countMatrixFullAllCalled$Metadata, design = ~Group)
ddsAD <- estimateSizeFactors(ddsAD)

ADcounts <- counts(ddsAD)

ADcounts<-as.data.frame(ADcounts)


ADcountList<-DGEList(counts=ADcounts, group=Metadata$Group)
ADcountList <- calcNormFactors(ADcountList)


MarziModelMatrix <- model.matrix(as.formula("~Agef + CETSif + Group"), data = Metadata)

ADcountList <- estimateDisp(ADcountList, MarziModelMatrix)

fitTMM <- glmQLFit(ADcountList, MarziModelMatrix)
qlf_group <- glmQLFTest(fitTMM, coef = "GroupAD")
group_results <- topTags(qlf_group, n = Inf) %>% data.frame() %>% mutate(PeakName = rownames(.))

groupResult_Anno <- group_results %>%
  AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName") %>% arrange(FDR)


#Heatmap of significant peaks
SignifCETsDF <- groupResult_Anno %>% filter(FDR < 0.05) %>% filter(!duplicated(PeakName)) %>% select(matches("PeakName|logFC|FDR|GSM"))
SignifCETsMatrix <- SignifCETsDF %>% select(matches("GSM"))


pheatmap()


#Repeat after accounting for Neuronal and Microglia MSP
MarziModelMatrixMSP <- model.matrix(as.formula("~Agef + NeuNall_MSP + Microglia_MSP + Group"), data = countMatrixFullAllCalled$Metadata)

ADcountList2 <- estimateDisp(ADcountList, MarziModelMatrixMSP)

fitTMM_MSP <- glmQLFit(ADcountList2, MarziModelMatrixMSP)
qlf_groupMSP <- glmQLFTest(fitTMM_MSP, coef = "GroupAD")
group_resultsMSP<- topTags(qlf_groupMSP, n = Inf) %>% data.frame() %>% mutate(PeakName = rownames(.)) 


groupResult_MSPAnno <- group_resultsMSP %>%
  AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName") %>% arrange(FDR)

SignifCETS <- group_results %>% filter(FDR < 0.05)
SignifMSP <- group_resultsMSP %>% filter(FDR < 0.05)

UniquePeakMSP <- SignifMSP %>% filter(!PeakName %in% SignifCETS$PeakName)
UniquePeakMSP <- groupResult_MSPAnno %>% filter(PeakName %in% UniquePeakMSP$PeakName) %>% select(PeakName, logFC, logCPM, PValue, FDR, symbol, Peak_Gene, GeneAnnoType, Peak.width, Peak.Location)

UniquePeakCETS <- SignifCETS %>% filter(!PeakName %in% SignifMSP$PeakName)


#Repeat after adjusting also for oligos
MarziModelMatrixMSPb <- model.matrix(as.formula("~Agef + NeuNall_MSP + Microglia_MSP + Oligo_MSP + Group"), data = countMatrixFullAllCalled$Metadata)

ADcountList2b <- estimateDisp(ADcountList, MarziModelMatrixMSP)

fitTMM_MSPb <- glmQLFit(ADcountList2b, MarziModelMatrixMSPb)
qlf_groupMSPb <- glmQLFTest(fitTMM_MSPb, coef = "GroupAD")
group_resultsMSPb <- topTags(qlf_groupMSPb, n = Inf) %>% data.frame() %>% mutate(PeakName = rownames(.)) 


groupResult_MSPAnno_b <- group_resultsMSPb %>%
  AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName") %>% arrange(FDR)



#Run the analysis similar to PD pipeline
Model = as.formula("~Agef + NeuNall_MSP + Microglia_MSP + Group")
DESeqOutAll_Full <- RunDESeq(data = countMatrixFullAllCalled$countMatrix, UseModelMatrix = T, MetaSamleCol = "GSM", SampleNameCol = "GSM", 
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg", sampleToFilter = "none",
                             FullModel = Model, test = "Wald", FitType = "local")
                                  

#DESegResultsSex_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "Age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "GroupAD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")

#Adjusting for oligos too
Modelb = as.formula("~Agef + NeuNall_MSP + Microglia_MSP + Oligo_MSP + Group")
DESeqOutAll_Fullb <- RunDESeq(data = countMatrixFullAllCalled$countMatrix, UseModelMatrix = T, MetaSamleCol = "GSM", SampleNameCol = "GSM", 
                             meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg", sampleToFilter = "none",
                             FullModel = Modelb, test = "Wald", FitType = "local")


#DESegResultsSex_FullAll <- GetDESeqResults(DESeqOutAll_Full, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge_FullAllb <- GetDESeqResults(DESeqOutAll_Fullb, coef = "Agef.L") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAllb <- GetDESeqResults(DESeqOutAll_Fullb, coef = "GroupAD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")


## Rerun with RLE normalization
DESeqOutAll_Full_RLE <- RunDESeq(data = countMatrixFullAllCalled$countMatrix, UseModelMatrix = T, MetaSamleCol = "GSM", SampleNameCol = "GSM", 
                                 meta = countMatrixFullAllCalled$Metadata, normFactor = NULL, sampleToFilter = "none",
                                 FullModel = Model, test = "Wald", FitType = "local")

#DESegResultsSex_FullAll_RLE <- GetDESeqResults(DESeqOutAll_Full_RLE, coef = "sexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge_FullAll_RLE <- GetDESeqResults(DESeqOutAll_Full_RLE, coef = "Age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAll_RLE <- GetDESeqResults(DESeqOutAll_Full_RLE, coef = "GroupAD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")



#### Rerun with CETS cell type correction
Model2 = as.formula("~Age + CETS + Group")

DESeqOutAll_Full_noCorrection <- RunDESeq(data = countMatrixFullAllCalled$countMatrix, UseModelMatrix = T, MetaSamleCol = "GSM", SampleNameCol = "GSM", 
                                           meta = countMatrixFullAllCalled$Metadata, normFactor = "MeanRatioOrg", sampleToFilter = "none",
                                           FullModel = Model2, test = "Wald", FitType = "local")


#DESegResultsSex_FullAll_noCorrection <- GetDESeqResults(DESeqOutAll_Full_noCorrection, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge_FullAll_noCorrection <- GetDESeqResults(DESeqOutAll_Full_noCorrection, coef = "Age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAll_noCorrection <- GetDESeqResults(DESeqOutAll_Full_noCorrection, coef = "GroupAD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")

#RLE normalization
DESeqOutAll_Full_RLE_noCorrection <- RunDESeq(data = countMatrixFullAllCalled$countMatrix, UseModelMatrix = T, MetaSamleCol = "GSM", SampleNameCol = "GSM", 
                                              meta = countMatrixFullAllCalled$Metadata, normFactor = NULL, sampleToFilter = "none",
                                              FullModel = Model2, test = "Wald", FitType = "local")


#DESegResultsSex_FullAll_RLE_noCorrection <- GetDESeqResults(DESeqOutAll_Full_RLE_noCorrection, coef = "SexM") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsAge_FullAll_RLE_noCorrection <- GetDESeqResults(DESeqOutAll_Full_RLE_noCorrection, coef = "Age") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")
DESegResultsGroup_FullAll_RLE_noCorrection <- GetDESeqResults(DESeqOutAll_Full_RLE_noCorrection, coef = "GroupAD") %>% AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName")




#Identify top genes
ModifiedDF <- DESegResultsGroup_FullAll
ModifiedDF$type <- sapply(ModifiedDF$type, function(x){
  rev(strsplit(x, "_")[[1]])[1]
}) %>% factor(levels = c("promoters", "exons", "5UTRs", "3UTR", "1to5kb","intronexonboundaries", "introns", "intergenic"))
ModifiedDF %<>% mutate(AnnoOrder = type)
levels(ModifiedDF$AnnoOrder) <- c(1:length(levels(ModifiedDF$AnnoOrder)))
ModifiedDF %<>% arrange(AnnoOrder) 
ModifiedDF %<>% filter(!duplicated(Peak_Gene))

SignifGenes <- ModifiedDF %>% filter(padj < 0.05) %>% .$symbol %>% unique
ModifiedDF %<>% filter(symbol %in% SignifGenes, !is.na(symbol), !is.na(padj))


TopGenes <- sapply(unique(ModifiedDF$symbol), function(gene){
  GenePeaks <- ModifiedDF %>% filter(symbol == gene)
  GenePeaks %<>% mutate(Adj2 = p.adjust(.$padj, method = "BH"))
  TotalPeaks <- nrow(GenePeaks)
  UpPeaks <- sum(GenePeaks$log2FoldChange > 0)
  DownPeaks <- sum(GenePeaks$log2FoldChange < 0)
  TopList <- if(UpPeaks/TotalPeaks == 1 | DownPeaks/TotalPeaks == 1 | "promoters" %in% c(GenePeaks %>%
                                                                                          filter(padj < 0.05) %>%
                                                                                          .$type %>% as.character())){
    "Yes"
  } else {
    "No"
  }
  NomSignif = sum(GenePeaks$pvalue < 0.05)
  AdjSignif = sum(GenePeaks$padj < 0.05)
  PromoterSignif <- if("promoters" %in% c(GenePeaks %>% filter(padj < 0.05) %>% .$type %>% as.character())){
    "Yes"
  } else if("promoters" %in% c(GenePeaks %>% filter(padj > 0.05) %>% .$type %>% as.character())){
    "No"
  } else {
    NA
  }
  Direction <- if(sum(GenePeaks$log2FoldChange) > 0){
    "Up"
  } else {
    "Down"
  }
  data.frame(GeneSymbol = gene,
             TotalPeaks = TotalPeaks,
             TopList = TopList,
             NomSignif = NomSignif,
             AdjSignif = AdjSignif,
             PromoterSignif = PromoterSignif,
             UpPeaks = UpPeaks,
             DownPeaks = DownPeaks,
             Direction = Direction,
             TopP = min(GenePeaks$Adj2)) %>% mutate(NomSignifProp = round(NomSignif/TotalPeaks, digits = 2),
                                               AdjSignifProp = round(AdjSignif/TotalPeaks, digits = 2))
}, simplify = F) %>% do.call(rbind, .) 


TopGenes %<>% filter(TopList == "Yes") %>% .[!grepl("^MIR|^SNOR", .$GeneSymbol),]

save.image(paste0(ResultsPath, Cohort, ".Rdata"))
saveRDS(DESeqOutAll_Full, file = paste0(ResultsPath, Cohort, "DEoutput.Rds"))
saveRDS(DESeqOutAll_Full_noCorrection, file = paste0(ResultsPath, Cohort, "DEoutputNoCorrection.Rds"))
saveRDS(DESeqOutAll_Full_RLE_noCorrection, file = paste0(ResultsPath, Cohort, "DEoutputRLENoCorrection.Rds"))
saveRDS(DESeqOutAll_Full_RLE, file = paste0(ResultsPath, Cohort, "DEoutputRLE.Rds"))
saveRDS(DESegResultsGroup_FullAll, file = paste0(ResultsPath, Cohort, "DEresults.Rds"))
saveRDS(DESegResultsGroup_FullAll_RLE, file = paste0(ResultsPath, Cohort, "DEresultsRLE.Rds"))
saveRDS(DESegResultsGroup_FullAll_noCorrection, file = paste0(ResultsPath, Cohort, "DEresultsNoCorrection.Rds"))
saveRDS(DESegResultsGroup_FullAll_RLE_noCorrection, file = paste0(ResultsPath, Cohort, "DEresultsRLENoCorrection.Rds"))

ADgene <- read.table("Data/GeneCards-SearchResults.csv", header = T, sep = ",")
PDgene <- read.table("../ChIPseqPD_reproduce/GeneralResults/PDgeneStat.tsv", header = T, sep = "\t")

library(minfi)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(dplyr)
library(wateRmelon)
library(meffil)
library(dplyr)
library(readxl)
library(data.table)
library(MASS) 
library(sandwich) 
library(lmtest) 
library(parallel) 
library(R.utils)
library(openxlsx)
library(dplyr)
library(qqman)

load("RG.set.rdata") # Load in your RG.set.rdata
load("pd.rdata")  # Load in your pd.rdata  

# Sample QC (remove bad quality and sex-mismatch samples)
# Function for removing bad quality samples 
record_failed_samples <- function(RG.set, Pheno.data, cutoff = 0.04) {
  det.p <- detectionP(RG.set)
  failed.probes <- det.p > 0.01
  failed.fraction <- colMeans(failed.probes)
  failed.fraction.df <- data.frame(Sample_Name = Pheno.data$Sample_Name, failed_fraction = failed.fraction)
  failed_samples <- failed.fraction.df[failed.fraction.df$failed_fraction >= cutoff, ]
  write.csv(failed_samples, file = "failed_samples_by_probe_frequency.csv")
  return(failed_samples$Sample_Name)
} 

fail.detp <- record_failed_samples(RG.set, Pheno.data)

# Function for removing sex-mismatch samples
check_gender_mismatch <- function(RG.set, Pheno.data) {
  MSet.raw <- preprocessRaw(RG.set)
  ratioSet <- ratioConvert(MSet.raw, what = "both", keepCN = TRUE)
  gset <- mapToGenome(ratioSet)
  predictedSex <- getSex(gset, cutoff = -2)
  gset <- addSex(gset, sex = predictedSex)
  sex.df <- data.frame(id = Pheno.data$Sample_Name, predict = predictedSex[,3], document = Pheno.data$Gender)
  sex.mismatch <- sex.df[sex.df$predict != sex.df$document, ] %>% pull(id)
  return(sex.mismatch)
  return(gset)
} 
sex.mismatch <- check_gender_mismatch(RG.set, Pheno.data)
ppi = 300
png("sexcheck.png", width = 5*ppi, height = 5*ppi, res = ppi)
plotSex(gset,id = gset$Sample_Name)
dev.off()

fail.sample <- unique(c(sex.mismatch, fail.detp)) 
filtered_RG.set <- RG.set[,-which(colnames(RG.set) %in% fail.sample)]
save(filtered_RG.set, file = "filtered_RG.set.rdata")


# Probes QC 
det.p <- detectionP(filtered_RG.set)
bad.probes <- rowMeans(det.p > 0.01) > 0.1
table(bad.probes)
bad.probe.names.detP <- rownames(det.p[bad.probes,]) 
length(bad.probe.names.detP) 

# get the annotation of the IlluminaHumanMethylationEPICv2manifest
anno = read.csv("EPIC-8v2-0_A1.csv", skip = 7)
probe_id_all_arrays = data.frame(IlmnID = anno$IlmnID, EPICv2_Name = anno$Name, EPICv1 = anno$EPICv1_Loci, M450 = anno$Methyl450_Loci)
# epic v1 bad probes
snp1 = read.csv("13059_2016_1066_MOESM4_ESM.csv") %>% filter(EUR_AF>0.05) %>% pull(PROBE)
snp2 = read.csv("13059_2016_1066_MOESM5_ESM.csv") %>% filter(EUR_AF>0.05) %>% pull(PROBE)
snp4 = read.csv("13059_2016_1066_MOESM1_ESM.csv") %>% pull(X)
crprobe = read.table("AppendixE_CrossReactiveProbes_EPICv1.txt", sep = "\t", head = T) %>% pull(Probe)
length(crprobe)
maskprobe = read.table("AppendixD_Zhou_et_al_MASKgeneral_list.txt", sep = "\t", head = F) %>% pull(V1)
length(maskprobe)
epicv1_remove_probes = unique(c(snp1, snp2, snp4, crprobe, maskprobe))
remove_id_g1 = probe_id_all_arrays %>% filter(EPICv1 %in% epicv1_remove_probes) %>% pull(IlmnID) 
length(remove_id_g1) 

bead_number = read.table("AppendixC_Probes_BeadNumber_outliers.txt", sep = "\t", head = T) %>% pull(Probe_Id)
length(bead_number) 
FlaggedProbes = read.table("AppendixG_FlaggedProbes_EPICv2.txt", sep = "\t", head = F) %>% pull(V1)
length(FlaggedProbes) 
MappingInaccuracy = read.table("AppendixH_MappingInaccuracy_EPICv2.txt", sep = "\t", head = F) %>% pull(V1)
length(MappingInaccuracy) 
epicv2_remove_probes <- c(bead_number, FlaggedProbes, MappingInaccuracy)
remove_id_g2 <- probe_id_all_arrays %>% filter((IlmnID %in% epicv2_remove_probes | EPICv2_Name %in% epicv2_remove_probes)) %>% pull(IlmnID)
length(remove_id_g2) 

remove_probe = unique(c(bad.probe.names.detP, remove_id_g1, remove_id_g2))

pd <- pData(filtered_RG.set)
MSetraw <- preprocessRaw(filtered_RG.set)
MSetrp <- MSetraw[!(rownames(MSetraw) %in% remove_probe),] 
gmsetrp <- mapToGenome(MSetrp)
MSet.sq <- preprocessQuantile(gmsetrp, sex = pd$Gender)

X <- anno$IlmnID[anno$CHR %in% c("chrX")]
Y <- anno$IlmnID[anno$CHR %in% c("chrY")]
remove.probexy <- unique(c(X,Y))
m.set.flt <- MSet.sq[!(rownames(MSet.sq) %in% remove.probexy),]
bval = getBeta(m.set.flt)
Mval <- getM(m.set.flt)

save(Mval, file = "Mval.rdata")
save(bval, file = "bval.rdata")











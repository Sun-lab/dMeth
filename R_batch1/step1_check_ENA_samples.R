

setwd("~/research/Deconvolution/original_data/")

# ---------------------------------------------------------------------------
# arrayexpress 450k
# ---------------------------------------------------------------------------

A450 = read.table("ArrayExpress_450k.txt", sep="\t", quote="", 
                  as.is=TRUE, header=TRUE)

dim(A450)
names(A450)
A450[1:2,]

sort(table(A450$Organism), decreasing=TRUE)
sort(table(substring(A450$Type,1,80)), decreasing=TRUE)

names(A450)

A450[grep("monocyte", A450$Title, ignore.case=TRUE),1:2]

A450[grep("neutrophils", A450$Title, ignore.case=TRUE),1:2]
# it seems E-GEOD-63499 is unpublished

A450[grep("neutrophils", A450$Title, ignore.case=TRUE),1:2]

A450[grep("T cell", A450$Title, ignore.case=TRUE),1:2]

A450[grep("Nature Killer", A450$Title, ignore.case=TRUE),1:2]
A450[grep("NK", A450$Title, ignore.case=TRUE),1:2]

# ---------------------------------------------------------------------------
# ENCODE
# ---------------------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/_ENA_by_cell_type")

E = read.table("ENCODE_report.tsv", sep="\t", header=TRUE)

dim(E)
names(E)
E[1:2,]

table(E$Assay.Type)
table(E$Biosample)

cellTypes = c("B cell", "CD14-positive monocyte")
E1 = E[which(E$Biosample %in% cellTypes),]
E1

# ---------------------------------------------------------------------------
# B cells
# ---------------------------------------------------------------------------

B = read.table("Bcell.txt", sep="\t", header=TRUE)

dim(B)
names(B)

nNA = colSums(is.na(B))
table(nNA)
B = B[,which(nNA < nrow(B))]

dim(B)
names(B)

B[1:2,]
table(B$cell_line)
B = B[which(B$cell_line == ""),]

dim(B)
names(B)
summary(B)

B1 = B[,c(1,4,9,5)]
B1 = B1[order(B1$accession),]
dim(B1)
B1

## manual inspection, most of them are not RNA-seq. the only ones that are 
## relevant are the last two rows

B[133:134,]

# GSM2400232
# GSM2400233
# GSM2400234
# GSM2400235

# ---------------------------------------------------------------------------
# T cells
# ---------------------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/_ENA_by_cell_type")

Tc = read.table("Tcell.txt", sep="\t", header=TRUE)

dim(Tc)
names(Tc)

nNA = colSums(is.na(Tc))
table(nNA)
Tc = Tc[,which(nNA < nrow(Tc))]

dim(Tc)
names(Tc)

Tc[1:2,]
table(Tc$cell_line)
Tc = Tc[which(Tc$cell_line == ""),]

dim(Tc)
names(Tc)
summary(Tc)

T1 = Tc[,c(1,4,10,6)]
T1 = T1[order(T1$accession),]
dim(T1)
T1

# the last 17 rows are relevant, inlcuding T cells from four healthy controls
# and 14 SLE patients
# T Cell Transcriptomes Describe Patient Subtypes in Systemic Lupus Erythematosus
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0141171

# ---------------------------------------------------------------------------
# dendritic cells
# ---------------------------------------------------------------------------

dendritic = read.table("dendritic.txt", sep="\t", header=TRUE)

dim(dendritic)
names(dendritic)

nNA = colSums(is.na(dendritic))
table(nNA)
dendritic = dendritic[,which(nNA < nrow(dendritic))]

dim(dendritic)
names(dendritic)
summary(dendritic)

dendritic

## manual inspection,
# if may be better to use the first 9 records which are the results of:
# 
# The role of antigen presenting cells in the induction of HIV-1 latency 
# in resting CD4+ T-cells, Retrovirology. 2015; 12: 76. 
# 
# the remaining sampels are in vitro differentiated dendritic cells
# from ENCODE
#
# myeloid  (m)DC  and  plasmacytoid  (p)DC 
# SLAN+ DC represent a subpopulation of monocytic cells with increased 
# potential to secrete pro-inflammatory cytokines and develop a DC phenotype, 
# however precise residence remains unknown. 
# CD14+ monocytes represent DC and macrophage precursors in blood

# ---------------------------------------------------------------------------
# Macrophage
# ---------------------------------------------------------------------------

macrophage = read.table("Macrophage.txt", sep="\t", header=TRUE)

dim(macrophage)
names(macrophage)

nNA = colSums(is.na(macrophage))
table(nNA)
macrophage = macrophage[,which(nNA < nrow(macrophage))]

dim(macrophage)
names(macrophage)

table(macrophage$cell_line)
macrophage = macrophage[which(macrophage$cell_line == ""),]
macrophage = macrophage[which(macrophage$tax_id == 9606),]

dim(macrophage)
names(macrophage)
summary(macrophage)

m1 = macrophage[,c(1,5,7)]
m1 = m1[order(m1$accession),]
m1$description = substring(m1$description, 1, 60)
dim(m1)
m1

# the 60 samples from this study may be good
# http://dx.doi.org/10.1371/journal.pgen.1006338
# Widespread Shortening of 3’ Untranslated Regions and Increased Exon Inclusion 
# Are Evolutionarily Conserved Features of Innate Immune Responses to Infection

# using primary human macrophages derived from whole blood samples from 60 
# individuals, we sequenced mRNA both before and after infection with two 
# live bacteria. 

# another study from this group
# Genetic ancestry and natural selection drive population differences in 
# immune responses to pathogens in humans
# Transcriptomic profiles of 503 infected (Listeria and Salmonella) and 
# non-infected samples at 2hr time point.
# GSE81046

mm1 = m1[which(m1$first_public=="2016-09-08" & grepl("mRNA", m1$description)),]
dim(mm1)

mm1 = mm1[grepl("_NI", mm1$description),]
dim(mm1)
mm1

# ---------------------------------------------------------------------------
# monocyte
# ---------------------------------------------------------------------------

monocyte = read.table("monocyte.txt", sep="\t", header=TRUE)

dim(monocyte)
names(monocyte)

nNA = colSums(is.na(monocyte))
table(nNA)
monocyte = monocyte[,which(nNA < nrow(monocyte))]

dim(monocyte)
names(monocyte)

table(monocyte$cell_line)
monocyte = monocyte[which(monocyte$cell_line == ""),]
monocyte = monocyte[which(monocyte$tax_id == 9606),]

dim(monocyte)
names(monocyte)
summary(monocyte)

m1 = monocyte[,c(1,5,7)]
m1 = m1[order(m1$accession),]
m1$description = substring(m1$description, 1, 60)
dim(m1)
m1

# the following three records are from the same paper for dendritic cells

# 142 SAMN03785020                    Donor 97 CD14+ mono   2015-12-25
# 143 SAMN03785024                    Donor 28 CD14+ Mono   2015-12-25
# 144 SAMN03785028                    Donor 36 CD14+ mono   2015-12-25


# the following three records are interesing since they study the interaction 
# of macrophage - cancer cell

# Characterization of macrophage - cancer cell crosstalk in estrogen receptor 
# positive and triple-negative breast cancer

# 170 SAMN04276593                                    M1C   2015-11-20
# 171 SAMN04276594                                    M1M   2015-11-20
# 172 SAMN04276595                                    M1T   2015-11-20
# 173 SAMN04276596                                    M2M   2015-11-20
# 174 SAMN04276597                                    M3C   2015-11-20
# 175 SAMN04276598                                    M3M   2015-11-20
# 176 SAMN04276599                                    M3T   2015-11-20

# another study that may have methylation data
# 
# "β-Glucan Reverses the Epigenetic State of LPS-Induced Immunological 
# Tolerance.", Cell, 2016 Nov 17;167(5):1354-1368.e14

# Two innate immune memory states can be induced in culture through an initial 
# exposure of primary human monocytes to either LPS or BG for 24 hours, 
# followed by removal of stimulus and differentiation to macrophages for an 
# additional 5 days. 
# ChIP-seq, RNA-seq, WGBS and ATAC-seq data were generated.


# ---------------------------------------------------------------------------
# neutrophils
# ---------------------------------------------------------------------------

neutrophils = read.table("neutrophils.txt", sep="\t", header=TRUE)

dim(neutrophils)
names(neutrophils)

nNA = colSums(is.na(neutrophils))
table(nNA)
neutrophils = neutrophils[,which(nNA < nrow(neutrophils))]

dim(neutrophils)
names(neutrophils)
summary(neutrophils)

n1 = neutrophils[,c(1,4,5)]
n1 = n1[order(n1$accession),]
n1$description = substring(n1$description, 1, 60)
dim(n1)
n1

# GSE60424
# Fresh blood samples were collected from healthy subjects and subjects 
# diagnosed type 1 diabetes, amyotrophic lateral sclerosis, and sepsis, 
# as well as multiple sclerosis patients. highly pure populations of 
# neutrophils, monocytes, B cells, CD4 T cells, CD8 T cells, and 
# natural killer cell

# GSE70068
# We report gene expression in human neutrophils isolated by two methods: 
# Polymorphprep (~95% purity) and negative selection (~99% purity) from two 
# healthy donors - one donor with low eosinophil contamination of neutrophils 
# and one donor with high eosinophil contamination of neutrophils. We report 
# the effect of the presence of contaminating leukocytes in neutrophil 
# preparations, and in reponse to inflammatory cytokines TNF-alpha and GM-CSF.

# ---------------------------------------------------------------------------
# Th and Tregs
# ---------------------------------------------------------------------------

Tr = read.table("Th_Treg.txt", sep="\t", header=TRUE, quote="")

dim(Tr)
names(Tr)

nNA = colSums(is.na(Tr))
table(nNA)
Tr = Tr[,which(nNA < nrow(Tr))]

dim(Tr)
names(Tr)
summary(Tr)

Tr1 = Tr[,c(1,3,4:6)]
Tr1 = Tr1[order(Tr1$accession),]
Tr1



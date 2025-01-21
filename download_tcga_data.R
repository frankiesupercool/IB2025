#All data used in this publication can be downloaded from
#https://tcga-data.nci.nih.gov/docs/publications/brca_2012/

dir.create("data", showWarnings = FALSE, recursive = TRUE)
#download gene expression data with subtype information
brca.rna.dest <- "data/BRCA.exp.547.med.txt"
brca.pam50 <- "data/BRCA.547.PAM50.SigClust.Subtypes.txt"
#download.file("http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.exp.547.med.txt", destfile=brca.rna.dest, method="wget")
download.file("https://api.gdc.cancer.gov/data/912fe02c-4f6f-4b42-a445-92dec4d4d09a", destfile=brca.rna.dest, method="wget")
#download.file("http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.547.PAM50.SigClust.Subtypes.txt", destfile=brca.pam50, method="wget")
download.file("https://api.gdc.cancer.gov/data/a2a299a9-6cad-4811-a74c-0bdfa61cbab2", destfile=brca.pam50, method="wget")

#download methylation data
temp <- tempfile()
#download.file("http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.methylation.27k.450k.zip", destfile=temp, method="wget")
download.file("https://api.gdc.cancer.gov/data/14bffe3a-4f07-401c-9178-d2d4c7aeb972", destfile=temp, method="wget")
unzip(temp, "BRCA.methylation.27k.450k.txt", exdir="data")
unlink(temp)

#download gene mapping data for methylation data
temp <- tempfile()
#download.file("https://tcga-data.nci.nih.gov/docs/publications/tcga/integration/adfs/tcga/jhu-usc.edu_TCGA_HumanMethylation27.v2.adf.txt.zip", destfile=temp, method="wget")
download.file("https://www.cancer.gov/ccg/research/structural-genomics/tcga/using-tcga-data/technology/illumina-humanmethylation27-adf", destfile=temp, method="wget")
unzip(temp, "jhu-usc.edu_TCGA_HumanMethylation27.v2.adf.txt", exdir="data")
unlink(temp)


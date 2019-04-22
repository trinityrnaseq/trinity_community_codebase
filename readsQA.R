#!/usr/bin/env Rscript
#Check if the packages required are alrady installed, and if it  is not, install them 
if (!require(fastqcr)) {
   install.packages("fastqcr")
   library(fastqcr)
}

if (!require(parallel)) {
   install.packages("parallel")
   library(parallel)
}

#Specify running variables
args <- commandArgs(TRUE)

if (length(args) < 2) {
   stop("Error: Incorrect number of arguments. It should be 2 or 3")
} else if (length(args)  > 3) {
   stop("Error: Incorrect number of arguments. It should be 2 or 3")
} else if (length(args)  == 2) {
   readsDir  <- normalizePath(args[1])
   outputDir <- normalizePath(args[2])
   nCoresRun <- c(4)
} else if (length(args)  == 3) {
   readsDir  <- normalizePath(args[1])
   outputDir <- normalizePath(args[2])
   nCoresRun <- as.numeric(args[3])
}

#Read the path to each reads file
readsPath <- list.files(path = readsDir, pattern = ".fastq.gz$", full.names = TRUE)

#Run fastqc for each reads file
fastQCfun     <- function(readsFile) {
   runCommand <- c("fastqc --noextract --nogroup --threads 1 --outdir")
   system(command = paste(runCommand, outputDir, readsFile))
}

mclapply(X = readsPath, FUN = fastQCfun, mc.cores = nCoresRun)

#Generate a summary of quality assestment
qcResult  <- qc_aggregate(qc.dir = outputDir)
qcSummary <- summary(qcResult)

#Export the final summary
write.table(x = qcSummary, file = paste0(outputDir, "/", "summary.tsv"), sep = "\t", col.names = TRUE, row.names = TRUE, append = FALSE)

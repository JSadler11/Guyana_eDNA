library(dada2)
library(ShortRead)
library(Biostrings)

path <- "/Volumes/ROSALIND/eDNA_Pilots/SCreek/raw/12SV5"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

FWD <- "ACTGGGATTAGATACCCC" 
REV <- "TAGAACAGGCTCCTCTAG"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients,
  primerHits, fn = fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
  fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

cutadapt <- "/Users/jacksadler/miniconda3/bin/cutadapt" 
system2(cutadapt, args = "--version") 

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

R1.flags <- paste("-g", FWD, "-a", REV.RC) 

R2.flags <- paste("-G", REV, "-A", FWD.RC) 

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients,primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

library(ShortRead)

for(f in cutFs[1:2]) {
  fs <- FastqStreamer(f, n = 1e5)  # read 100k reads at a time
  fq_chunk <- yield(fs)
  qs <- as(quality(fq_chunk), "matrix")
  avg_q <- colMeans(qs, na.rm = TRUE)
  plot(avg_q, type="l", xlab="Cycle", ylab="Mean Quality", main=basename(f))
  close(fs)
}

for(f in cutRs[1:2]) {
  fs <- FastqStreamer(f, n = 1e5)  # read 100k reads at a time
  fq_chunk <- yield(fs)
  qs <- as(quality(fq_chunk), "matrix")
  avg_q <- colMeans(qs, na.rm = TRUE)
  plot(avg_q, type = "l", xlab = "Cycle", ylab = "Mean Quality", 
       main = paste("Reverse:", basename(f)))
  close(fs)
}

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

table(nchar(getSequences(seqtab.nochim)))

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
               rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste0(">ASV_", seq_len(length(asv_seqs)))

# Write FASTA file
asv_fasta <- c(rbind(asv_headers, asv_seqs))
writeLines(asv_fasta, file.path(path, "ASVs.fa"))

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- paste0("ASV_", seq_len(nrow(asv_tab)))

write.csv(asv_tab, file.path(path, "ASV_counts.csv"))

#Run in Terminal
#blastn -query /Volumes/ROSALIND/eDNA_Pilots/SCreek/raw/12SV5/ASVs_for_BLAST.fa \
#-db /Volumes/ROSALIND/blastdb \
#-out /Volumes/ROSALIND/eDNA_Pilots/SCreek/raw/12SV5/ASV_BLAST_results.tsv \
#-evalue 1e-20 \
#-outfmt "6 qseqid sacc staxids pident length mismatch gapopen evalue bitscore stitle" \
#-max_target_seqs 5 \
#-num_threads 8

blast_res <- read.table(file.path(path, "ASV_BLAST_results.tsv"), sep="\t", header=FALSE,
                        col.names = c("ASV", "Accession", "TaxID", "PercentID", "Length", "Mismatches",
                                      "GapOpen", "Evalue", "Bitscore", "Title"))
head(blast_res)

asv_biom <- make_biom(data = asv_tab)
write_biom(asv_biom, biom_file = file.path(path, "ASV_table.biom"))

#qiime tools import \
#--input-path ASV_table.biom \
#--type 'FeatureTable[Frequency]' \
#--input-format BIOMV100Format \
#--output-path feature-table.qza

#qiime tools import \
#--input-path ASVs.fa \
#--output-path rep-seqs.qza \
#--type 'FeatureData[Sequence]'


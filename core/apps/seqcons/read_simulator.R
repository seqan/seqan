# ============================================================================
#                SeqAn - The Library for Sequence Analysis
#                          http://www.seqan.de 
# ============================================================================
#  Copyright (C) 2007
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 3 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  Lesser General Public License for more details.
# ============================================================================
# Author: Tobias Rausch <rausch@embl.de>
# ============================================================================

############################################################
# Variables
############################################################

# Path variables
directory = "./"
readFile = "reads.fasta"				# Multi-fasta file with all reads
fragmentFile = "reads.fastaF"				# All fragments where each fragment is a single read or 2 reads (for mate pairs)
libraryFile = "reads.fastaL"				# All libraries or a dummy library if no mate pairs were simulated
sourceFile = "reads.fastaS"				# The source sequence (or multiple sequences if haplotypes have been simulated)

# Basic parameters
alphabet = c('A','C','G','T') # The source sequence alphabet
seqLength = 100000		# Source sequence length
sourceSeq = ""			# Source sequence, if empty a random sequence is siumlated
numOfReads = 1250			# Number of reads to simulate	
fixedReadLength = 0		# 0 = false, avgReadLength applies; 1 otherwise

# Global error rate, fixedReadLength = 0
avgReadLength = 800		# Average read length
errorPerBaseCall = 0.02		# Sequencing error rate for each base in a read (substitution or indels)

# Or: Error distribution, fixedReadLength = 1
# Absolute probability for each base to mutate (only substitutions)
errorDist = c(0.001131099, 0.001613963, 0.002547571, 0.002702844, 0.003357026, 0.003673833, 0.004200899, 0.005458078, 0.007530188, 0.00973986, 0.01253226, 0.01396316, 0.01820574, 0.02151224, 0.02786259, 0.03494583, 0.04295635, 0.05040179, 0.06265487, 0.07137801, 0.08404933, 0.09705024, 0.1051062, 0.1209309, 0.1361337, 0.1483654, 0.1658058, 0.176571, 0.1886209, 0.2049304, 0.2115647, 0.2261853, 0.2333330, 0.2503056, 0.2637948, 0.2351953)
readsPerBucket = c(0.20, 0.30, 0.50)	# 20% with 0 errors, 30% with 1 error, and 50% with 2 errors
							# If empty "c()" all reads are put into their corresponding bucket, fractional content is shown at the end

# Source sequence repeat parameters (implanted repeats)
numOfRepeats = 0			# Repeat occurences
repeatLength = 20			# Repeat length of the implanted repeat, avgRepeatLength * numOfRepeats << seqLength
					# Repeats are simulated from a common root with numOfRepeats children and x SNPs and y Indels
repeatSNPRate = 0.01		# SNP rate
repeatIndelRate = 0.02		# Indel rate
repeatIndelBaseRange = 4:6	# Indel size range


# Haplotype parameters (overcompressed repeats)
					# Haplotypes are simulated from a common root with numOfHaplotypes children
numOfHaplotypes = 1		# Number of haplotypes
snpRate = 0.01			# SNP rate
indelRate = 0.01			# Indel rate
indelBaseRange = 4:6		# Indel size range


# Mate-pairs parameters
simulateMatePairs = 1		# 0 = no mate pairs are simulated, 1 = all reads are in a mate pair
librarySizes = c(10000,20000)	# Mean library sizes for mate pairs
librarySd = c(10,20)		# Standard deviations for above library sizes

# Path for reference file if generated.
referencePath = "reference.fasta"

############################################################
# Optional command line arguments override above settings
# Usage:  R --slave  --vanilla --args directory=~/matches/reads/test/ seqLength=100 numOfReads=10 avgReadLength=36 simulateMatePairs=0 fixedReadLength=0 errorPerBaseCall=0.4 < read_simulator.R
############################################################
meanLib = 0
meanSd = 0
args <- commandArgs()
for(arg in args) {
	if (length(grep("=", arg))) {
		keyvalue = strsplit(arg, "=")
		keyvalue = unlist(keyvalue)
		if (keyvalue[1] == "directory") directory = keyvalue[2]
		if (keyvalue[1] == "sourceSeq") sourceSeq = keyvalue[2]
		if (keyvalue[1] == "seqLength") seqLength = as.integer(keyvalue[2])
		if (keyvalue[1] == "numOfReads") numOfReads= as.integer(keyvalue[2])
		if (keyvalue[1] == "avgReadLength") avgReadLength = as.integer(keyvalue[2])
		if (keyvalue[1] == "errorPerBaseCall") errorPerBaseCall = as.double(keyvalue[2])
		if (keyvalue[1] == "fixedReadLength") fixedReadLength = as.integer(keyvalue[2])
		if (keyvalue[1] == "simulateMatePairs") simulateMatePairs = as.integer(keyvalue[2])
		if (keyvalue[1] == "meanLib") meanLib = as.integer(keyvalue[2])
		if (keyvalue[1] == "meanSd") meanSd = as.integer(keyvalue[2])
	}
}
rm(arg)
rm(args)

############################################################
# Print the current setting to the command line
############################################################
print(ls.str())


############################################################
# Functions
############################################################


####
# Creates a random sequence of a given length and over a given alphabet
####
createRandomSequence=function(alphabet = c('A','C','G','T'), seqLength = 100, p = rep(1/length(alphabet), length(alphabet))) {
	sample(alphabet, seqLength, replace = TRUE, prob=p)
}

####
# Creates a random sequence of a given length and over a given alphabet with repeat copies
####
createRandomSequenceWithRepeats=function(alphabet = c('A','C','G','T'), seqLength = 100, numOfRepeats = 2, repeatLength = 40, repeatSNPRate = 0.01, repeatIndelRate = 0.01, repeatIndelBaseRange = c(1)) {
	if (numOfRepeats < 2) {
		seq = createRandomSequence(alphabet, seqLength)
	} else {
		seq = createRandomSequence(alphabet, repeatLength)
		if (repeatLength * numOfRepeats > seqLength) numOfRepeats = floor(seqLength / repeatLength)
		
		# Create polymorphic copies (with x SNPs and y Indels)
		repeats = c()
		for(i in 1:numOfRepeats ) {
			seqPoly = createPolymorphicCopy(alphabet, seq, repeatSNPRate / 2, repeatIndelRate / 2, repeatIndelBaseRange)
			repeats = c(repeats, list(seqPoly))
		}

		# Plug repeats together
		interRepeatDist = floor((seqLength - (repeatLength * numOfRepeats)) / (numOfRepeats + 1))
		seq = createRandomSequence(alphabet, interRepeatDist)
		for(i in 1:numOfRepeats ) {
			seq = paste(c(seq, repeats[[i]]),sep = "")
			seq = paste(c(seq, createRandomSequence(alphabet, interRepeatDist)),sep = "")
		}
		seq
	}
}

####
# Creates a polymorphic sequence from a given sequence
####
createPolymorphicCopy = function(alphabet = c('A','C','G','T'), seqPoly = createRandomSequence(), snpRate = 0.01, indelRate = 0.01, indelBaseRange = c(1)) {
	seqLength = length(seqPoly)
	snp_pos = (1:seqLength)[runif(seqLength) <= snpRate]
	if (length(snp_pos) > 0) {
		for(j in snp_pos) seqPoly[j] = sample(alphabet[alphabet!= seqPoly[j]], 1)
	}
	indels = length((1:seqLength)[runif(seqLength) <= indelRate])
	if (indels > 0) {
		for(j in 1:indels) {
			seqLength = length(seqPoly)
			pos = sample((1:seqLength), 1)
			indelSize = sample(indelBaseRange, 1)
			if (runif(1) > 0.5) { # Insertion
				beginPiece = c()
				if (pos > 1) beginPiece = seqPoly[1:pos - 1]
				seqPoly= c(beginPiece, sample(alphabet, indelSize, replace = TRUE), seqPoly[pos:seqLength])
			} else { # Deletion
				beginPiece = c()
				if (pos > 1) beginPiece = seqPoly[1:pos - 1]
				endPiece = c()
				if (pos + indelSize <= seqLength) {
					end = pos + indelSize
					endPiece = seqPoly[end:seqLength]
				}
				seqPoly= c(beginPiece, endPiece)
			}
		}
	}
	seqPoly
}


####
# Reverse complements a given sequence
####
reverseComplement = function(seq, alphabet = c('A','C','G','T','N'), compAlphabet = c('T','G','C','A','N')) {
	seq = rev(seq)
	for(i in 1:length(alphabet)) {
		seq[seq == alphabet[i]] = i
	}
	seq = as.integer(seq)
	for(i in 1:length(alphabet)) {
		seq[seq == i] = compAlphabet[i]
	}
	seq
}

####
# Simulates a set of reads from a set of haplotypes
####
simulateReads=function(alphabet = c('A','C','G','T'), haplotypes = list(createRandomSequence(), createRandomSequence()), numOfReads = 100, avgReadLength = 35, errorPerBaseCall = 0.01, simulateMatePairs = 1, librarySizes = c(90, 80), librarySd = c(5, 2), readPath, fragmentPath) {
	tmpPath = paste(c(readPath, "tmp"),sep = "", collapse = "")
	unlink(readPath)
	unlink(tmpPath)
	unlink(fragmentPath)
	if (simulateMatePairs == 1) {
		if (numOfReads %% 2 != 0) numOfReads = numOfReads +1
		numOfReads = numOfReads / 2
	}
	for(readCounter in 1:numOfReads) {
		# Sample a read length
		readLength = round(rnorm(1, mean = avgReadLength, sd = avgReadLength / 10))

		# Pick a haplotype
		currentHaplotype = sample(1:(length(haplotypes)), 1)

		# Sample a read
		invalidMatePair = 1
		while(invalidMatePair == 1) {		
			sourceSeq = haplotypes[[currentHaplotype]]
			seqLength = length(sourceSeq)	
			start = sample(1:(seqLength - readLength + 1), 1)
			end = start + readLength - 1
			read = sourceSeq[start:end]

			# Reverse complement this read
			if (runif(1) > 0.5) {
				read = reverseComplement(read)
				tmp = start
				start = end
				end = tmp
			}
				
			if (simulateMatePairs == 0) {
				currentLibrary = 1
				invalidMatePair = 0
			} else {
				# Pick a library size
				currentLibrary = sample(1:(length(librarySizes)), 1)
				libSize = round(rnorm(1, mean=librarySizes[currentLibrary], sd=librarySd[currentLibrary]))
				if (start < end) {
					if (start + libSize <= seqLength) {
						startMatePair = start + libSize - readLength + 1 
						endMatePair = startMatePair + readLength - 1
						readMate = sourceSeq[startMatePair:endMatePair]
						readMate = reverseComplement(readMate)
						tmp = startMatePair 
						startMatePair = endMatePair 
						endMatePair = tmp
						invalidMatePair = 0
					}
				} else {
					if (start - libSize >= 1) {
						startMatePair = start - libSize  
						endMatePair = startMatePair + readLength - 1
						readMate = sourceSeq[startMatePair:endMatePair]
						invalidMatePair = 0
					}
				}
			}
		}
		
		# Sequence the reads
		read = createPolymorphicCopy(alphabet, read, errorPerBaseCall / 2, errorPerBaseCall / 2)
		if (simulateMatePairs == 1) readMate = createPolymorphicCopy(alphabet, readMate, errorPerBaseCall / 2, errorPerBaseCall / 2)

		# Simulate library noise
		if ((meanLib != 0) && (meanSd != 0)) {
			sourceSeqLen = length(haplotypes[[currentHaplotype]])	
			redoEverything = 1
			while(redoEverything == 1) {
				offsetRead = round(rnorm(1, meanLib , meanSd) - meanLib)
				if ((start + offsetRead >= 0) && (end + offsetRead >= 0) &&	(start + offsetRead < sourceSeqLen) && (end + offsetRead < sourceSeqLen)) {
					start = start + offsetRead 
					end = end + offsetRead 
					redoEverything = 0
				}
			}
		}

		# Add read to file
		if (start < end) header = paste(">",start - 1,sep="")
		else header = paste(">",end + length(read) - 1,sep="")
		header = paste(header,",", sep="")
		if (start < end) header = paste(header,start + length(read) - 1,sep="")
		else header = paste(header,end - 1,sep="")
		header = paste(header,"[id=", sep="")
		header = paste(header,readCounter - 1, sep="")
		header = paste(header,",fragId=", sep="")
		header = paste(header,readCounter - 1, sep="")
		# Add the haplotype to the header so it can be identified later
		header = paste(header,",repeatId=", sep="")
		header = paste(header,currentHaplotype - 1, sep="")
		header = paste(header,"]", sep="")
		read = paste(read,sep = "",collapse = "")
		write(header, file=readPath, sep="\n", append = TRUE)
		write(read, file=readPath, sep="\n", append = TRUE)

		# Add fragment
		header = paste(">",readCounter - 1,sep="")
		header = paste(header,"[libId=", sep="")
		header = paste(header,currentLibrary - 1, sep="")
		header = paste(header,"]", sep="")
		if (simulateMatePairs == 1) {
			read = paste((numOfReads + readCounter) - 1,readCounter - 1,sep=",",collapse = "")
		} else {
			read = paste("0", "0",sep=",",collapse = "")
		}
		write(header, file=fragmentPath, sep="\n", append = TRUE)
		write(read, file=fragmentPath, sep="\n", append = TRUE)

		if (simulateMatePairs == 1) {
			if (startMatePair < endMatePair) header = paste(">",startMatePair - 1,sep="")
			else header = paste(">",endMatePair + length(readMate) - 1,sep="")
			header = paste(header,",", sep="")
			if (startMatePair < endMatePair) header = paste(header,startMatePair + length(readMate) - 1,sep="")
			else header = paste(header,endMatePair - 1,sep="")
			header = paste(header,"[id=", sep="")
			header = paste(header,(numOfReads + readCounter - 1), sep="")
			header = paste(header,",fragId=", sep="")
			header = paste(header,readCounter - 1, sep="")
			# Add the haplotype to the header so it can be identified later
			header = paste(header,",repeatId=", sep="")
			header = paste(header,currentHaplotype - 1, sep="")
			header = paste(header,"]", sep="")
			read = paste(readMate,sep = "",collapse = "")
			write(header, file=tmpPath, sep="\n", append = TRUE)
			write(read, file=tmpPath, sep="\n", append = TRUE)
		}
	}
	if (simulateMatePairs == 1) {
		write(scan(tmpPath, what = 'character'), file=readPath, sep="\n", append = TRUE)
		unlink(tmpPath)
	}
}

####
# Simulates a set of reads from a set of haplotypes with a certain error distribution
####
simulateReadsFromErrorDist=function(alphabet = c('A','C','G','T'), haplotypes = list(createRandomSequence(), createRandomSequence()), numOfReads = 100, errorDist = c(0.01, 0.02, 0.03), readsPerBucket = c(0.50, 0.50), simulateMatePairs = 1, librarySizes = c(90, 80), librarySd = c(5, 2), readPath, fragmentPath) {
	tmpPath = paste(c(readPath, "tmp"),sep = "", collapse = "")
	unlink(readPath)
	unlink(tmpPath)
	unlink(fragmentPath)
	readCounter = 0
	readLength = length(errorDist)
	if (length(readsPerBucket) < 1) readsPerBucket = rep(1, length(errorDist) + 1)
	bucketCounter = rep(0, length(readsPerBucket))
	if (simulateMatePairs == 1) {
		if (numOfReads %% 2 != 0) numOfReads = numOfReads +1
		numOfReads  = numOfReads / 2
	}
	while (readCounter < numOfReads) {
		# Pick a haplotype
		currentHaplotype = sample(1:(length(haplotypes)), 1)

		# Sample a read
		invalidMatePair = 1
		while(invalidMatePair == 1) {		
			sourceSeq = haplotypes[[currentHaplotype]]
			seqLength = length(sourceSeq)	
			start = sample(1:(seqLength - readLength + 1), 1)
			end = start + readLength - 1
			read = sourceSeq[start:end]

			# Reverse complement this read
			if (runif(1) > 0.5) {
				read = reverseComplement(read)
				tmp = start
				start = end
				end = tmp
			}
				
			if (simulateMatePairs == 0) {
				currentLibrary = 1
				invalidMatePair = 0
			} else {
				# Pick a library size
				currentLibrary = sample(1:(length(librarySizes)), 1)
				libSize = round(rnorm(1, mean=librarySizes[currentLibrary], sd=librarySd[currentLibrary]))
				if (start < end) {
					if (start + libSize <= seqLength) {
						startMatePair = start + libSize - readLength + 1 
						endMatePair = startMatePair + readLength - 1
						readMate = sourceSeq[startMatePair:endMatePair]
						readMate = reverseComplement(readMate)
						tmp = startMatePair 
						startMatePair = endMatePair 
						endMatePair = tmp
						invalidMatePair = 0
					}
				} else {
					if (start - libSize >= 1) {
						startMatePair = start - libSize  
						endMatePair = startMatePair + readLength - 1
						readMate = sourceSeq[startMatePair:endMatePair]
						invalidMatePair = 0
					}
				}
			}
		}


		# Sequence the reads
		countErrors = 0
		for(pos in 1:length(read)) {
			if (runif(1) <= errorDist[pos]) {
				read[pos] = sample(alphabet[alphabet!= read[pos]], 1)
				countErrors = countErrors + 1
			}
		}
		countMateErrors = 0
		if (simulateMatePairs == 1) {
			for(pos in 1:length(readMate)) {
				if (runif(1) <= errorDist[pos]) {
					readMate[pos] = sample(alphabet[alphabet!= readMate[pos]], 1)
					countMateErrors = countMateErrors + 1
				}
			}
		}
	
		if ((countErrors < length(readsPerBucket)) &&
		    (countMateErrors < length(readsPerBucket))) {
			if ((simulateMatePairs == 1) || ((bucketCounter[countErrors+1] / numOfReads) < readsPerBucket[countErrors+1])) {
				bucketCounter[countErrors+1] = bucketCounter[countErrors+1] + 1
				readCounter = readCounter + 1
				
				# Add read to file
				if (start < end) header = paste(">",start - 1,sep="")
				else header = paste(">",end + length(read) - 1,sep="")
				header = paste(header,",", sep="")
				if (start < end) header = paste(header,start + length(read) - 1,sep="")
				else header = paste(header,end - 1,sep="")
				header = paste(header,"[id=", sep="")
				header = paste(header,readCounter - 1, sep="")
				header = paste(header,",fragId=", sep="")
				header = paste(header,readCounter - 1, sep="")
				# Add the haplotype to the header so it can be identified later
				header = paste(header,",repeatId=", sep="")
				header = paste(header,currentHaplotype - 1, sep="")
				header = paste(header,",errors=", sep="")
				header = paste(header,countErrors, sep="")
				header = paste(header,"]", sep="")
				read = paste(read,sep = "",collapse = "")
				write(header, file=readPath, sep="\n", append = TRUE)
				write(read, file=readPath, sep="\n", append = TRUE)

				# Add fragment
				header = paste(">",readCounter - 1,sep="")
				header = paste(header,"[libId=", sep="")
				header = paste(header,currentLibrary - 1, sep="")
				header = paste(header,"]", sep="")
				if (simulateMatePairs == 1) {
					read = paste((numOfReads + readCounter - 1),readCounter - 1,sep=",",collapse = "")
				} else {
					read = paste("0", "0",sep=",",collapse = "")
				}
				write(header, file=fragmentPath, sep="\n", append = TRUE)
				write(read, file=fragmentPath, sep="\n", append = TRUE)

				if (simulateMatePairs == 1) {
					bucketCounter[countMateErrors+1] = bucketCounter[countMateErrors+1] + 1

					if (startMatePair < endMatePair) header = paste(">",startMatePair - 1,sep="")
					else header = paste(">",endMatePair + length(readMate) - 1,sep="")
					header = paste(header,",", sep="")
					if (startMatePair < endMatePair) header = paste(header,startMatePair + length(readMate) - 1,sep="")
					else header = paste(header,endMatePair - 1,sep="")
					header = paste(header,"[id=", sep="")
					header = paste(header,(numOfReads + readCounter - 1), sep="")
					header = paste(header,",fragId=", sep="")
					header = paste(header,readCounter - 1, sep="")
					# Add the haplotype to the header so it can be identified later
					header = paste(header,",repeatId=", sep="")
					header = paste(header,currentHaplotype - 1, sep="")
					header = paste(header,",errors=", sep="")
					header = paste(header,countMateErrors, sep="")
					header = paste(header,"]", sep="")
					read = paste(readMate,sep = "",collapse = "")
					write(header, file=tmpPath, sep="\n", append = TRUE)
					write(read, file=tmpPath, sep="\n", append = TRUE)				}
			}
		}
	}
	if (simulateMatePairs == 1) {
		print(bucketCounter / (2* numOfReads))
	} else {
		print(bucketCounter / numOfReads)
	}

	if (simulateMatePairs == 1) {
		write(scan(tmpPath, what = 'character'), file=readPath, sep="\n", append = TRUE)
		unlink(tmpPath)
	}
}



############################################################
# Main
############################################################


# Set all path variables
readPath = paste(c(directory, readFile),sep = "", collapse = "")
fragmentPath = paste(c(directory, fragmentFile),sep = "", collapse = "")
sourcePath = paste(c(directory, sourceFile),sep = "", collapse = "")
libraryPath = paste(c(directory, libraryFile),sep = "", collapse = "")


# Create random source sequence with or without repeats
if (nchar(sourceSeq) != 0) {
	data = scan(sourceSeq, skip = 1, what = 'character', strip.white = TRUE)
	seq = paste(data ,sep = "", collapse = "")
	if (any(grep('>', seq ))) {
		print("Multiple Sequences in this file!!!")
	}
	seq = unlist(strsplit(seq,""))
} else {
	if (numOfRepeats < 2) {
		seq = createRandomSequence(alphabet, seqLength)
	} else {
		seq = createRandomSequenceWithRepeats(alphabet, seqLength, numOfRepeats, repeatLength, repeatSNPRate, repeatIndelRate, repeatIndelBaseRange)
	}
    write(">reference", file=referencePath, append = FALSE)
    write(paste(seq, sep = "", collapse = ""), file=referencePath, append = TRUE)
}


# Create polymorphic copies (with x SNPs and y Indels)
haplotypes = c()
if (numOfHaplotypes < 2) {
	haplotypes = c(list(seq))
	numOfHaplotypes = 1
} else {
	for(i in 1:numOfHaplotypes) {
		seqPoly = createPolymorphicCopy(alphabet, seq, snpRate / 2, indelRate / 2, indelBaseRange)
		haplotypes = c(haplotypes, list(seqPoly))
	}
}

# Simulate reads
if (fixedReadLength == 0) {
	simulateReads(alphabet, haplotypes, numOfReads, avgReadLength, errorPerBaseCall, simulateMatePairs, librarySizes, librarySd, readPath, fragmentPath)
} else {
	simulateReadsFromErrorDist(alphabet, haplotypes, numOfReads, errorDist, readsPerBucket, simulateMatePairs, librarySizes, librarySd, readPath, fragmentPath)
}

# Write libraries
if (simulateMatePairs == 1) {
	writeString = ""
	for(i in 1:length(librarySizes)) {
		writeString = c(writeString, ">")
		writeString = c(writeString, i - 1)
		writeString = c(writeString, "\n")
		writeString = c(writeString, librarySizes[i])
		writeString = c(writeString, ",")
		writeString = c(writeString, librarySd[i])
		writeString = c(writeString, "\n")
		writeString = paste(writeString,sep = "", collapse = "")
		write(writeString, file=libraryPath, append = FALSE)
	}
} else {
	# Write a dummy lib
	writeString = ""
	writeString = c(writeString, ">")
	writeString = c(writeString, 0)
	writeString = c(writeString, "\n")
	writeString = c(writeString, 0)
	writeString = c(writeString, ",")
	writeString = c(writeString, 0)
	writeString = c(writeString, "\n")
	writeString = paste(writeString,sep = "", collapse = "")
	write(writeString, file=libraryPath, append = FALSE)
}


# Write haplotypes
writeString = ""
for(i in 1:numOfHaplotypes) {
	writeString = c(writeString, ">Sequence")
	writeString = c(writeString, i - 1)
	writeString = c(writeString, "\n")
	writeString = c(writeString, haplotypes[[i]])
	writeString = c(writeString, "\n")
}
writeString = paste(writeString,sep = "", collapse = "")
write(writeString, file=sourcePath, append = FALSE)


rm(list=ls())


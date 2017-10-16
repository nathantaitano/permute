source("fstPairwisePermutationFun.R")

# tempDir <- "../temp/"                                                # For running at home
tempDir <- "/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/"  # For running on Zcluster
outDir <- "../../gbsOutputData/"

tempPairs <- c(paste(tempDir, "tus1Tus2Fst", sep=""),
               paste(tempDir, "tus1CagFst", sep=""),
               paste(tempDir, "tus1CosFst", sep=""),
               paste(tempDir, "tus1FruFst", sep=""),
               paste(tempDir, "tus2CagFst", sep=""),
               paste(tempDir, "tus2CosFst", sep=""),
               paste(tempDir, "tus2FruFst", sep=""),
               paste(tempDir, "cagCosFst", sep=""),
               paste(tempDir, "cagFruFst", sep=""),
               paste(tempDir, "cosFruFst", sep=""))

outPairs <- c(paste(outDir, "tus1Tus2Fst", sep=""),
              paste(outDir, "tus1CagFst", sep=""),
              paste(outDir, "tus1CosFst", sep=""),
              paste(outDir, "tus1FruFst", sep=""),
              paste(outDir, "tus2CagFst", sep=""),
              paste(outDir, "tus2CosFst", sep=""),
              paste(outDir, "tus2FruFst", sep=""),
              paste(outDir, "cagCosFst", sep=""),
              paste(outDir, "cagFruFst", sep=""),
              paste(outDir, "cosFruFst", sep=""))

for(i in 1:length(tempPairs)){
  critFst <- getCritFsts(tempPairs[i])
  write(critFst, paste(outPairs[i], "CritVals.txt", sep=""))
}
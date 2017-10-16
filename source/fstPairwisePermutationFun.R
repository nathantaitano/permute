makePairPermuteSamples <- function(popFile1, popFile2, prefix, n=1000){
  pop1 <- readLines(popFile1)
  pop2 <- readLines(popFile2)
  combinedPop <- c(pop1, pop2)
  
  for(p in 1:n){
  permutPop1 <- sample(combinedPop, length(pop1)) # Uses existing population sizes to make the permutation populations
  permutPop2 <- sample(combinedPop, length(pop2))
  writeLines(permutPop1,
             paste(prefix,p,"pop1.txt",sep = "."))
  writeLines(permutPop2,
             paste(prefix,p,"pop2.txt",sep = "."))
  }
}

getCritFsts <- function(prefix){
  require(stringr)
  permFiles <- system(paste("ls ",
                            paste(prefix, ".perm*.sum", sep="")), intern = T)
  weightFst <- c()
  meanFst <- c()
  for(i in 1:length(permFiles)){
    permFileLines <- readLines(permFiles[i])
    
    weightFstLine <- permFileLines[grepl("weighted",permFileLines)]
    weightFst <- c(weightFst, str_extract(weightFstLine, "(?<=: )[-\\.\\d]+$"))
  }
  weightFst <- sort(as.numeric(weightFst))
  
  quantiles <- quantile(weightFst, c(.95, .99, .999))
  ps <- sapply(quantiles, function(x) head(weightFst[weightFst>x], n=1))
  
  return(c(p.05=ps[1], p.01=ps[2], p001=ps[3]))
}
#!/bin/bash
# Uses 16 threads for alignment, and 20 gb mem. qsub command:
# qsub -q rcc-mc-30d -pe thread 16 -l mem_total=20g AddRefs.sh

tempDir="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/"
# tempDir="../temp/"
inDir="../../gbsInputData/"

outDir="../../gbsOutputData/"


refVcfHeadLine="##SAMPLE=<ID=XRef:9997,Barcode=NNNNNN,Column=12,Flowcell=NNNNNNNNNN,Flowcell_Lane=NNNNNNNNN_both,Lane=both,LibraryPrepID=9997,PlateName=A,Row=A,Sample=XRef,Status=private>"
alt1VcfHeadLine="##SAMPLE=<ID=chiltepin:9998,Barcode=NNNNNN,Column=0,Flowcell=NNNNNNNNNN,Flowcell_Lane=NNNNNNNNN_both,Lane=both,LibraryPrepID=9998,PlateName=A,Row=A,Sample=chiltepin,Status=private>"
alt2VcfHeadLine="##SAMPLE=<ID=zunla:9999,Barcode=NNNNNN,Column=0,Flowcell=NNNNNNNNNN,Flowcell_Lane=NNNNNNNNN_both,Lane=both,LibraryPrepID=9999,PlateName=A,Row=A,Sample=zunla,Status=private>"
refVcfHeadCol="XRef:9997"
alt1VcfHeadCol="chiltepin:9998"
alt2VcfHeadCol="zunla:9999"

snpContextRFile="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/addRefSnpContext.Rfasta"
snpContextFile="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/addRefSnpContext.fasta"
zunlaSaiFile=${tempDir}zunla.sai
chiltSaiFile=${tempDir}chilt.sai
zunlaSamFile=${tempDir}zunla.sam
chiltSamFile=${tempDir}chilt.sam
zunlaGenomeFile="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/addRef/zunla/Capsicum.annuum.L_Zunla-1_Release_2.0.fasta"
zunlaBt2Prefix="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/addRef/zunla/zunla"
chiltGenomeFile="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/MUMmerXRef/Capsicum.annuum.var.glabriusculum_Chiltepin_Release_2.0.fasta"
chiltBt2Prefix="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/addRef/chiltepin/chiltepin"
cm334GenomeFile="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/refgenome/Pepper_1.55_reform.fasta"

origVcf=${outDir}mikeyLmiss20Allhet05lmiss20imiss20.recode.vcf
vcfHeader="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/vcfHeader.txt"
rvcfFileAllGenomesAdded="/escratch4/nkt58510/nkt58510_Jun_03/Peppers/gitTemp/mikeyLmiss20Allhet05lmiss20imiss20Cm334ChiltZunla.Rvcf"
vcfFileAllGenomesAdded=${outDir}mikeyLmiss20Allhet05lmiss20imiss20Cm334ChiltZunla.vcf

export PATH=${PATH}:/usr/local/R/3.3.2/bin/
export PATH=${PATH}:/usr/local/vcftools/latest/bin/


# arg 1:population pair prefix arg2:vcf arg3:number of permutations
function calcPermFst {
  p=1
  while (( $p <= $3 ))
  do
    vcftools --vcf $2 \
    --weir-fst-pop $1.${p}.pop1.txt --weir-fst-pop $1.${p}.pop2.txt \
    --out $1.perm${p}.out 2> $1.perm${p}.sum
    ((p++))
  done
}

Rscript makePairwisePermPops.R
#
calcPermFst ${tempDir}tus1Tus2Fst $vcfFileAllGenomesAdded 10000
calcPermFst ${tempDir}tus1CagFst $vcfFileAllGenomesAdded 10000
calcPermFst ${tempDir}tus1CosFst $vcfFileAllGenomesAdded 10000
calcPermFst ${tempDir}tus1FruFst $vcfFileAllGenomesAdded 10000
#
calcPermFst ${tempDir}tus2CagFst $vcfFileAllGenomesAdded 10000
calcPermFst ${tempDir}tus2CosFst $vcfFileAllGenomesAdded 10000
calcPermFst ${tempDir}tus2FruFst $vcfFileAllGenomesAdded 10000
#
calcPermFst ${tempDir}cagCosFst $vcfFileAllGenomesAdded 10000
calcPermFst ${tempDir}cagFruFst $vcfFileAllGenomesAdded 10000
#
calcPermFst ${tempDir}cosFruFst $vcfFileAllGenomesAdded 10000

# Extracts fst values for each permutation from vcftools output, gets the .05, .01, .001 percentiles for fst p-values
Rscript getCritFstTusTest.R

Rscript makeSpidPopFileAllAccs.R

vcftools --vcf $vcfFileAllGenomesAdded --keep ../../gbsOutputData/allPassingKeepfile.txt \
  --maf .0001 --max-missing .9 --recode --out ${vcfFileAllGenomesAdded}SegOnly 
Removes monomorphic loci (30061 remain). Necessary for conversion to structure format in PGDSpider
# Also removes individual accessions with low coverage, keeping only Oaxaca individuals (not the reference indvs) with high coverage and SNPs with high cov as well.
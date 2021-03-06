# Generate fst and heterozygosity
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase31 \
  --fst Population \
  --het \
  --nonfounders \
  --out ${WORKDIR}/1000G/all_phase31

# Generate population specific allele frequencies
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase31 \
  --freq \
  --nonfounders \
  --loop-cats Population \
  --out ${WORKDIR}/1000G/all_phase31

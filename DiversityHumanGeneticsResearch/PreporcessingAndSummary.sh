# The following code was run on University of Florida Cluster Computer HiPerGator

# Decompress p files
./plink2 \
  --zst-decompress ${WORKDIR}/1000G/allphase3.pgen.zst \
  > ${WORKDIR}/1000G/allphase3.pgen
  
./plink2 \
  --zst-decompress ${WORKDIR}/1000G/allphase3.pvar.zst \
  > ${WORKDIR}/1000G/allphase3.pvar

# Remove multiallelic variants, and variants and samples with over 0.02 missing genotyping rate
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase3 \
  --max-alleles 2 \
  --geno 0.02 \
  --mind 0.02 \
  --make-pgen \
  --out ${WORKDIR}/1000G/all_phase31
   
# Remove nonautosomes
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase31 \
  --autosome \
  --snps-only \
  --make-pgen \
  --out ${WORKDIR}/1000G/all_phase31
  
# Assign IDs to variants with missing marker IDs
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase31 \
  --set-missing-var-ids @:# \
  --make-pgen \
  --out ${WORKDIR}/1000G/all_phase31
  
# Remove duplicate marker IDs
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase31 \
  --rm-dup force-first \
  --make-pgen \
  --out ${WORKDIR}/1000G/all_phase31
  
# Prune Variants
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase31 \
  --indep-pairwise 200 50 0.25 \
  --out ${WORKDIR}/1000G/all_phase31
  
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase31 \
  --exclude ${WORKDIR}/1000G/all_phase31.prune.out \
  --make-pgen \
  --out ${WORKDIR}/1000G/all_phase31
  
# Remove first degree individuals
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase31 \
  --make-king \
  --king-cutoff 0.177\
  --make-pgen \
  --out ${WORKDIR}/1000G/all_phase31
  
# Generate missingness report
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase31 \
  --missing \
  --out ${WORKDIR}/1000G/all_phase31

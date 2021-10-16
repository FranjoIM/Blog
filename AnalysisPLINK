# Generate fst and heterozygosity
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase31 \
  --fst Population \
  --out ${WORKDIR}/1000G/all_phase31

# Generate population specific allele frequencies
./plink2 \
  --pfile ${WORKDIR}/1000G/all_phase31 \
  --freq \
  --nonfounders \
  --loop-cats Population \
  --out ${WORKDIR}/1000G/all_phase31

# Calculate gene diversity for each population
SUM="${WORKDIR}/1000G/SUM"
SUMMARY="${WORKDIR}/1000G/SUMMARY"

mkdir -p ${SUM} S{SUMMARY}

cat ${WORKDIR}/1000G/all_phase31.psam | tr ' ' ',' | tr '\t' ',' |awk 'NR > 1' | awk -F ',' '{print $6}' | sort | uniq > ${WORKDIR}/1000G/populations.txt

cat ${WORKDIR}/1000G/populations.txt | while read item || [[ -n $line ]]
do
  awk '{print 1-(1-$5)**2-$5**2}' ${WORKDIR}/1000G/all_phase31.${item}.afreq > ${WORKDIR}/1000G/all_phase31.${item}.freq 
  paste -sd+ ${WORKDIR}/1000G/all_phase31.${item}.freq | bc > ${SUM}/all_phase31.${item}.freq.sum
  cat ${WORKDIR}/1000G/all_phase31.${item}.freq | datamash min 1 q1 1 median 1 q3 1 max 1 iqr 1 mean 1 sstdev 1 > ${SUMMARY}/all_phase31.${item}.freq.summary
done;

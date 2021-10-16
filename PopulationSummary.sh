# Summarize number of samples per superpopulation
cat ${WORKDIR}/1000G/phase3_corrected.psam | tr ' ' ',' | tr '\t' ',' | awk -F ',' '{print $5}' | sort | uniq -c

# Summarize number of samples per population
cat ${WORKDIR}/1000G/phase3_corrected.psam | tr ' ' ',' | tr '\t' ',' | awk -F ',' '{print $6}' | sort | uniq -c

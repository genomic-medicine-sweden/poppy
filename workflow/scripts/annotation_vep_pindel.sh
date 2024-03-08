#!/usr/bin/env bash

if ( zcat "${snakemake_input[vcf]}" | grep -q -v '^#' );
 then
 echo 'lines found' &&
 (vep --vcf --no_stats -o "${snakemake_output[vcf]}" -i "${snakemake_input[vcf]}" --dir_cache "${snakemake_input[cache]}" \
 --fork "${snakemake[threads]}" ${snakemake_params[mode]} --fasta "${snakemake_input[fasta]}" \
 ${snakemake_params[extra]} ) &> "${snakemake_log[0]}";
else
 echo 'no lines found' &&
 cp "${snakemake_input[vcf]}" "${snakemake_output[vcf]}" &&
 gzip -d "${snakemake_output[vcf]}".gz ;
fi

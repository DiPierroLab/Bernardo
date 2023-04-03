#!/bin/bash

c=7
joblim=64
cell_type=xx
region=95.4-96.5
nkeep=(5400 5250 5000 4500 4000 3000)
output_dir='/scratch/b.zubillagaherrera/dna_lig/amit_data/var_cutting_scheme/output/'
new_output_dir='/scratch/b.zubillagaherrera/dna_lig/amit_data/var_cutting_scheme/output/rc1.5'
rc=1.5
gamma_vals=(1)
nstart=1
nfiles=20
size=5500
for keep in ${nkeep[@]}
do
        for gamma in ${gamma_vals[@]}
        do
     		fn_base=C_${c}'_region_'${region}'_nkeep_'${keep}'_from_manystr'
     		runfile=mean_first_time_ligation_analysis_${fn_base}.py
     		cat combine_ft_lig_distoflig_data_source.py \
     		    | sed  -e "s|@chr|${c}|g" \
     		      -e "s|@region|${region}|g" -e "s|@nkeep|${keep}|g" -e "s|@nstart|${nstart}|g" -e "s|@size|${size}|g" -e "s|@new_output_dir|${new_output_dir}|g"\
     		      -e "s|@output_dir|${output_dir}|g" -e "s|@nfiles|${nfiles}|g" -e "s|@gamma|${gamma}|g" -e "s|@rc|${rc}|g" \
     		 	> ${runfile}
     			jobname=submit_mean_ligation_analysis_${fn_base}
     			batchfile=submit_${jobname}.sbatch
     			cp submit_analyze_source ${batchfile}
     			sed -i "s|@jobname|${jobname}|g" ${batchfile}
     			sed -i "s|@runfile|${runfile}|g" ${batchfile}
     			#sbatch ${batchfile}
	done
done

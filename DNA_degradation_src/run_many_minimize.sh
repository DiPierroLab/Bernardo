#!/bin/bash

c=7
joblim=20
cell_type=xx
region=95.4-96.5
output_dir='/scratch/amit.das/dna_lig/chr7_95.4-96.5/data'
wallpos=300
nstart=1
nfiles=20
k=0
ll=251
dt_for_fire=0.0001
dt_for_langvd=0.0001
maxirun=1000
njobs=4007
jobname=chr_${c}'_region_'${region}'_minimize_list_'${ll}
joblist_file=${jobname}.txt
batchfile=submit_${jobname}.sbatch
rm ${joblist_file}
echo " " > ${joblist_file}
ijob=0
for n in $(seq $nstart 1 $nfiles)
do
	for j in $(seq 2 5 999)   # (0 20 999)
	do
		fn_base=C_${c}'_region_'${region}'_file_'${n}'_'str_${j}
        	input_name=minimized_${fn_base}
        	finaloutfile=${output_dir}/${input_name}'_'minimized.gsd
        	if [ -f $finaloutfile ]; then
        	        echo run already done
        	        continue
        	fi

		runfile=minimize_energy_${fn_base}.py
		seed1=$RANDOM
		seed2=$RANDOM
        	input_traj=input/REP${n}/eq_chr7_95.4-96.5_conf_0.pdb
		cat minimize_chromosome_configuration_source.py \
			| sed -e "s|@data_dir|${output_dir}|g" -e "s|@wallpos|${wallpos}|g" -e "s|@snap_index|${j}|g" \
		 	-e "s|@seed1|${seed1}|g" -e "s|@seed2|${seed2}|g" -e "s|@dt_for_fire|${dt_for_fire}|g" -e "s|@dt_for_langvd|${dt_for_langvd}|g"\
		 	-e "s|@input_traj|${input_traj}|g" -e "s|@input_name|${input_name}|g" -e "s|@maxirun|${maxirun}|g" \
		 	> ${runfile}
        	let k=k+1
        	echo "python" ${runfile} "> out_"${fn_base} >> ${joblist_file}
                let ijob=ijob+1
        	if [ $k -eq ${joblim} ] || [ $ijob -eq $njobs ]; then
        		cp submit_chr_all_source ${batchfile}
        		sed -i "s|@jobname|${jobname}|g" ${batchfile}
        		sed -i "s|@joblist|${joblist_file}|g" ${batchfile}
        		sed -i "s|@ntasks|${joblim}|g" ${batchfile}
        		sbatch ${batchfile}
        		let k=0
        		let ll=ll+1
        		jobname=chr_${c}'_region_'${region}'_minimize_list_'${ll}
			joblist_file=${jobname}.txt
        		batchfile=submit_${jobname}.sbatch
        		rm ${joblist_file}
        		echo " " > ${joblist_file}
        	fi
	done
done


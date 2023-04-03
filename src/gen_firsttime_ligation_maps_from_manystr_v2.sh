#!/bin/bash

c=7
joblim=64
cell_type=xx
region=95.4-96.5
rc=1.5
nkeep=(5000) #(5400 5250 5000 4500 4000 3000)
gamma_vals=(1)
input_dir='/work/dipierrolab/amit.das/dna_lig/chr7_95.4-96.5/minimized_str_files/data'
output_dir='/scratch/b.zubillagaherrera/Fr_ij/output'
new_output_dir='/scratch/b.zubillagaherrera/Contact_Distance_Cycle_Ensemble/rc'${rc} #'_v2'
nstart=1
nsim=1 #5000
nfiles=20
njobs=$((5000*${nsim}))  #$((${nsim}))   #$((5620*${nsim}))
k=0
ll=1

ijob=0

for keep in ${nkeep[@]}
do
	jobname=ligan_chr_${c}'_region_'${region}'_nkeep_'${keep}'_from_manystr_list_'${ll}
	joblist_file=${jobname}.txt
	batchfile=submit_${jobname}.sbatch
	rm ${joblist_file}
	echo " " > ${joblist_file}
        for gamma in ${gamma_vals[@]}
        do
        	for n in $(seq $nstart 1 $nfiles)
        	do
			    for sim in $(seq 1 1 $nsim)
                do
        	        for j in $(seq 0 1 999)
        	        do
					    fn_base=C_${c}'_region_'${region}'_nkeep_'${keep}'_str_'${j}'_n_'${n}'_sim_'${sim}
        	            traj_file=${output_dir}/chromosome_crosslinked_fixed_chopped_${fn_base}.gsd

					           #echo $traj_file

        	            if [ -f ${traj_file} ]; then
						echo ${traj_file}
						runfile=first_time_ligation_analysis_${fn_base}.py
 						seed1=$RANDOM
						cat ft_lig_distoflig_analysis_source_v2.py \
					 		| sed -e "s|@seed_num|${seed1}|g" -e "s|@fn_base|${fn_base}|g" -e "s|@chr|${c}|g" \
                                                 	-e "s|@filenumber|${n}|g" -e "s|@sim|${sim}|g" \
						 	-e "s|@region|${region}|g" -e "s|@output_dir|${output_dir}|g" -e "s|@str|${j}|g" \
						 	-e "s|@input_dir|${input_dir}|g" -e "s|@rc|${rc}|g" -e "s|@new_output_dir|${new_output_dir}|g" \
					 		> ${runfile}
        		        		let k=k+1
        	                		let ijob=ijob+1
        	                		echo $ijob $njobs
        	                		echo "python" ${runfile} "> out_"${fn_base} >> ${joblist_file}
        	                		if [ $k -eq ${joblim} ] || [ $ijob -eq $njobs ]; then
        	                	        	echo $ijob
        		        			cp submit_chr_all_source ${batchfile}
        		        			sed -i "s|@jobname|${jobname}|g" ${batchfile}
        		        			sed -i "s|@joblist|${joblist_file}|g" ${batchfile}
        		        			sed -i "s|@ntasks|${joblim}|g" ${batchfile}
        		        			#sbatch ${batchfile}
        		        			let k=0
        		        			let ll=ll+1
        		        			jobname=ligan_chr_${c}'_region_'${region}'_nkeep_'${keep}'_from_manystr_list_'${ll}
        		        			joblist_file=${jobname}.txt
        		        			batchfile=submit_${jobname}.sbatch
        		        			rm ${joblist_file}
        		        			echo " " > ${joblist_file}
						fi
					fi
				done
			done
		done
	done
done

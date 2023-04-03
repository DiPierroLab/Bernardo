#!/bin/bash

c=7
dt=0.001
nsteps=1000000
nsteps_eq=0
seed_number=$RANDOM
freq=1000
joblim=80
cell_type=xx
region=95.4-96.5
nkeep=(4000 4500 3000) #(5400 5250 5000 4500 4000 3000) # chopped bonds = N_bonds - nkeep
input_dir='/work/dipierrolab/amit.das/dna_lig/chr7_95.4-96.5/minimized_str_files/data'
output_dir='/scratch/b.zubillagaherrera/Fr_ij_2/output'
wallpos=300
nevolve=5000  # n_crosslink = N - nevolve
#suffix=rela_k562_hg19
nstart=1
nfiles=20
runs_per_str=1
prev_run_per_str=0
njobs=$((16860*${runs_per_str}))
k=0
ll=1
jobname=chr_${c}'_region_'${region}'_chop_evolve_multirun_list_'${ll}
joblist_file=${jobname}.txt
batchfile=submit_${jobname}.sbatch
rm ${joblist_file}
echo " " > ${joblist_file}
ijob=0
for keep in ${nkeep[@]}
do
	for n in $(seq $nstart 1 $nfiles)
	do
        	for j in $(seq 0 1 999)   # (0 20 999)
        	do
                        run_start=$((${prev_run_per_str}+1))
                        run_end=$((${run_start}+${runs_per_str}-1))
                        for rr in $(seq ${run_start} 1 ${run_end})
			do
                		fn_base=C_${c}'_region_'${region}'_file_'${n}'_'str_${j}
       	        		input_name=minimized_${fn_base}
                		minimized_str_file=${input_dir}/${input_name}'_'minimized.gsd
        			if [ -f ${minimized_str_file} ]; then
                			#echo input file ${minimized_str_file} found
					new_fn_base=C_${c}'_region_'${region}'_nkeep_'${keep}'_file_'${n}'_str_'${j}'_run_'${rr}
					runfile=run_chr_chop_evolve_${new_fn_base}.py
					seed1=$RANDOM
					seed2=$RANDOM
					cat chromosome_fix_chop_evolve_from_singlestr_source.py \
						| sed -e "s|@dt|${dt}|g" -e "s|@nsteps|${nsteps}|g" \
						-e "s|@output_dir|${output_dir}|g" -e "s|@wallpos|${wallpos}|g" \
						-e "s|@nevolve|${nevolve}|g" -e "s|@minimized_str_file|${minimized_str_file}|g"\
					 	-e "s|@seed1|${seed1}|g" -e "s|@seed2|${seed2}|g" -e "s|@fn_base|${new_fn_base}|g" \
					 	-e "s|@freq|${freq}|g" -e "s|@chr|${c}|g" -e "s|@filenumber|${n}|g" -e "s|@str|${j}|g"\
        	        		        -e "s|@teq|${nsteps_eq}|g" -e "s|@region|${region}|g" -e "s|@nkeep|${keep}|g" -e "s|@run|${rr}|g" \
					 	> ${runfile}
        	        		let k=k+1
					let ijob=ijob+1
					#echo $ijob $njobs
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
						jobname=chr_${c}'_region_'${region}'_chop_evolve_multirun_list_'${ll}
        	        			joblist_file=${jobname}.txt
        	        			batchfile=submit_${jobname}.sbatch
        	        			rm ${joblist_file}
        	        			echo " " > ${joblist_file}
		        		fi
				else
					continue
				fi
			done
		done
	done
done

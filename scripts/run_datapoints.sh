while read line; do 
	echo "#!/bin/bash">run_${line}.sh;
        echo "#SBATCH --nodes=1">>run_${line}.sh;
        echo "#SBATCH --ntasks=1">>run_${line}.sh;
        echo "#SBATCH --cpus-per-task=4">>run_${line}.sh;
	echo "#SBATCH --mem=20GB">>run_${line}.sh;
	echo "#SBATCH --time=00:01:00">>run_${line}.sh;
	echo "python3 ../code/efficient.py ../data/datapoints/${line} ../outputs/${line}_output.txt" >> run_${line}.sh
done <../data/datapoints/inputs.txt

#PBS -N gTFIv2
#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu
#PBS -e /Users/azofeifa/correlation_networks/
#PBS -o /Users/azofeifa/correlation_networks/
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=64
#PBS -l mem=10gb

<<<<<<< HEAD
./CORR join -bed_dir /Users/joazofeifa/Lab/new_motif_distances/EMG_out_files_human/ -o /Users/joazofeifa/Desktop/test.bed -fast 0 -win\
dow 1000 -distance 10000
=======
#./CORR join -bed_dir /Users/joazofeifa/Lab/new_motif_distances/EMG_out_files_human/ -o /Users/joazofeifa/Desktop/test.bed -fast 0 -window 1000 -distance 10000

>>>>>>> b1d03cd208e9e161240b3ff215a1da4515a0859b

#./CORR insert -bed /Users/joazofeifa/Lab/correlation_networks/files/joined_bidir.bed  -o /Users/joazofeifa/Desktop/coverage_stats_2.bed  -bg_dir /Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files_2/

./CORR correlate -bed /Users/joazofeifa/Lab/correlation_networks/files/coverage_stats_3000.tsv -o ~/Desktop/correlations.txt

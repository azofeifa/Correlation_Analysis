#PBS -N gTFIv2
#PBS -m ae
#PBS -M joseph.azofeifa@colorado.edu
#PBS -e /Users/azofeifa/correlation_networks/
#PBS -o /Users/azofeifa/correlation_networks/
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=64
#PBS -l mem=10gb

#./CORR join -bed_dir /Users/joazofeifa/Lab/new_motif_distances/EMG_out_files_human/ -o /Users/joazofeifa/Desktop/test.bed -fast 0 -window 1000 -distance 10000


#./CORR insert -bed /Users/joazofeifa/Desktop/joined_bed_prelim_10000_5000  -o /Users/joazofeifa/Desktop/coverage_stats_2.bed  -bg_dir /Users/joazofeifa/Lab/gro_seq_files/HCT116/bed_graph_files_2/ -test 1

./CORR correlate -bed /Users/joazofeifa/Desktop/SRR942423_coverage_all.tsv -o ~/Desktop/correlations.txt -test 1 -threshold -1

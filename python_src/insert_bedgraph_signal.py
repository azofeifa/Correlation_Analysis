def load_bed_file(FILE):
	F 	= open(FILE, 'r')
	for line in FH:
		if line[0]!="#":
			chrom,start, stop 	= line.strip("\n").split("\t")[:3]
			if chrom not in G:
				G[chrom]=list()
			G[chrom].append([ int(start), int(stop),{} ]))
	return G
def insert_signal(FILE, G, ID,SUMS):
	SUMS[ID] 	= 0.0
	for chrom in G:
		for i in range(len(G[chrom])):
			G[chrom][i][2][ID] 	= 0.0
	SUM,prevchrom 	= 0.0, ""
	FH 	= open(FILE, 'r')
	for line in FH:
		chrom,start, stop, cov 	= line.strip("\n").split("\t")
		if chrom!=prevchrom and chrom in G:
			j,N 				= 0, len(G[chrom])
		while j < N and G[chrom][j][1] < int(start):
			j+=1
		y 	= abs(int(stop) - int(start))*abs(int(cov))
		if j < N and G[chrom][j][0] < int(stop):
			G[chrom][j][2][ID]+=y
		SUMS[ID]+=y

		prevchrom 				= chrom
	FH.close()
	return G
def iterate_through_bedgraph(root,G,OUT):
	SUMS 	= {}
	for f in os.listdir(root):
		ID 	= f.split(".fastq")[0]
		G 	= insert_signal(root+f, G, ID,SUMS)
	FHW 	= open(OUT, "w")
	SRRS 	= G.keys()
	FHW.write(",".join(SRRSs)+ "\t" + ",".join([ str(SUMS[S]) for S in SRRSs]) + "\n")
	for chrom in G:
		for j in range(len(G[chrom])):
			FHW.write(chrom+"\t" + str(G[chrom][j][0]) +"\t"+ str(G[chrom][j][1]) +"\t" + ',',join([ str(G[chrom][j][2][SRR]) for SRR in SRRSs]) +"\n"  )
	return G






if __name__ == "__main__":
	FILE 			= "/Users/azofeifa/bidirectional_hits_all/joined_bed_files/"
	bedgraph_root 	= "/scratch/Shares/dowell/processed_pubgro/human_SRR_rawcoverage/"
	G 				= load_bed_file(FILE)
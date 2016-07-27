import numpy as np
import matplotlib.pyplot as plt
import time
import os
import multiprocessing as mp
import seaborn as sns
def joint(points, distance):
	joined=True
	points.sort()
	while joined and len(points)>1 :
		joined 	= False
		d,i 	= min([(points[i+1] - points[i],i) for i in range(len(points)-1)])
		if d < distance:
			joined 	= True
			c 		= (points[i] + points[i+1])/2.
			points 	= points[:i] + [c] + points[i+2:]
	points 	= [int(p) for p in points]
	return points
def joint2(points, distance):
	points.sort()
	new_points 	= list()
	i 			= 0
	while i < len(points):
		new_points.append(points[i])
		while i < len(points) and (points[i] - new_points[-1]) < distance:
			new_points[-1] 	= (new_points[-1] + points[i])/2.
			i+=1
#		i+=1
	new_points 	= [int(n) for n in new_points]
	return new_points



def insert_bed_file(F, G):
	FH 	= open(F,'r')
	for line in FH:
		if line[0]!="#":
			chrom,start, stop 	= line.strip("\n").split("\t")[:3]
			if chrom not in G:
				G[chrom]=list()
			x 	= (float(start) + float(stop)) / 2.
			G[chrom].append(x)
	return G 

def load_all(DIR):
	G 	= {}
	for k,f in enumerate(os.listdir(DIR)):
		G 	= insert_bed_file(DIR+f, G)
		if k > 10:
			break
	return G
def wrapper(OUT, distance, G):
	FHW 	= open(OUT+"joined_"+str(distance)+".bed", 'w')
	chroms 	= G.keys()[:4]
	for chrom in chroms:
		points 	= joint(G[chrom],distance)	
		for x in points:
			FHW.write(chrom + "\t" + str(int(x)-1500) + "\t" + str(int(x)+1500) + "\n")
	FHW.close()
	pass


def join_all(G,OUT):
	jobs 	= list()
	for distance in np.linspace(1,500,3):
		p = mp.Process(target=wrapper, args=(OUT, int(distance), G))
		jobs.append(p)
		p.start()
	for j in jobs:
		j.join()
	

	return G


if __name__ == "__main__":
	y1,y2 		= list(),list()
	for i in range(300):
		print i
		l,n 	= 1000,200
		distance = 10
		npoints 	= list(np.random.randint(0,l,n))
		
		n1 			= len(joint(npoints, distance))
		
		n2 			= len(joint2(npoints, distance))

		y1.append(n1)
		y2.append(n2)
	F 	= plt.figure()
	plt.hist(y1,bins=30,alpha=0.5,label="MIN")
	plt.hist(y2,bins=30,alpha=0.5,label="GREED")
	plt.legend(loc="best")
	plt.show()
	# DIR 	= "/Users/joazofeifa/Lab/new_motif_distances/EMG_out_files_human/"
	# OUT 	= "/Users/joazofeifa/Desktop/test_stuff/"
	# print "loading...",
	# G 		= load_all(DIR)
	# print "done"
	# join_all(G,OUT)


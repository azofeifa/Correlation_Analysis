import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def load(FILE,test=True):
	header 	= True
	X 		= list()
	t 		= 0
	with open(FILE) as FH:
		for line in FH:
			if header:
				header  				= False
				RAW 					= [float(x.split("|")[1]) for x in line.strip("\n").split(",")]

			elif not header:
				chrom,start, stop,info 	= line.strip("\n").split("\t")
				info 					= [ sum([float(y) for y in x.split("|")]) / RAW[i] for i,x in enumerate(info.split(","))]
				X.append(info)
				if test and t > 1000:
					break
				t+=1
	return X
def compute_distance(X,OUT):
	FHW 	= open(OUT, "w")
	for i in range(len(X)):
		line 	= "1.0,"
		for j in range(i+1,len(X)):
			corr 	= np.corrcoef(X[i], y=X[j])
			line+=str(corr[0,1])+ ","
		line = line.strip(",")
		FHW.write(line+"\n")
	FHW.close()
def load_out(FILE):
	Y 			= list()
	with open(FILE) as FH:
		for line in FH:
			X 	= [float(x) for x in line.strip("\n").split(",")]
			Y.append(X)
	return Y
def get_histogram(Y):
	return [x for X in Y for x in X]
def make_network(Y,threshold=0.8,weight=0.05,add=False):
	G 	= nx.Graph()

	for i in range(len(Y)):
		if add:
			G.add_node(i)
		for j in range(1, len(Y[i])):
			if Y[i][j] > threshold:
				G.add_edge(i,j+i,weight=weight)
				G.add_edge(j+i,i,weight=weight)
	return G
def display(Y):
	weight 			= 0.03

	node_size 		= 50
	node_alpha 		= 0.5
	node_color 		= "blue"

	edge_tickness 	= 0.5
	edge_alpha 		= 0.5
	edge_color 		= "black"

	scores 	= get_histogram(Y)
	F 	= plt.figure()
	ax1 = F.add_subplot(2,2,1)
	ax1.hist(scores,bins=30, color="green", edgecolor="white")
	ax1.set_xlabel("Pearons Correlation Coefficient")
	ax1.set_ylabel("Frequency")

	ax2 = F.add_subplot(2,2,2)
	G 	= make_network(Y,threshold=0.95,weight=weight)
	pos = nx.spring_layout(G)	
	nx.draw(G, pos,ax=ax2,node_size=15)
	nx.draw_networkx_nodes(G,pos,node_size=node_size, 
		alpha=node_alpha, node_color=node_color)
	nx.draw_networkx_edges(G,pos,width=edge_tickness,
		alpha=edge_alpha,edge_color=edge_color)

	ax2.set_title("Network Thresholded, > 0.95")



	ax3 = F.add_subplot(2,2,3)
	ax3.set_title("Network Thresholded, > 0.55")
	G 	= make_network(Y,threshold=0.55,weight=weight)
	pos = nx.spring_layout(G)	
	nx.draw(G, pos,ax=ax3,node_size=15)
	nx.draw_networkx_nodes(G,pos,node_size=node_size, 
		alpha=node_alpha, node_color=node_color)
	nx.draw_networkx_edges(G,pos,width=edge_tickness,
		alpha=edge_alpha,edge_color=edge_color)



	ax4 = F.add_subplot(2,2,4)

	penalties 	= np.linspace(0.0, 1.0,20)

	counts 		= [len([x for x in nx.connected_components(make_network(Y,threshold=p, weight=1,add=True))]) for p in penalties]
	ax4.plot(penalties, counts)
	ax4.scatter(penalties, counts)
	ax4.set_xlabel("Pearon's Treshold")
	ax4.set_ylabel("Number of Connected Components")

	plt.tight_layout()
	plt.show()


if __name__ == "__main__":
	REMAKE 	= False
	OUT 	= "/Users/joazofeifa/Lab/correlation_networks/files/pearons.tsv"
		
	if REMAKE:
		FILE 	= "/Users/joazofeifa/Lab/correlation_networks/files/counts_2_current.tsv"
		X 		= load(FILE)
		compute_distance(X,OUT)
	X 			= load_out(OUT)
 	display(X)


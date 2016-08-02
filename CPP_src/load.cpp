#include "load.h"

#include "split.h"
#include <stdio.h>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <string>
#include <sstream>
using namespace std;
//========================================================================
//The very very important segment class

segment::segment(string chr, int st, int sp,vector<string> IDS){
	chrom 	= chr, start = st, stop = sp;
	for (int i = 0 ; i < IDS.size();i++){
		G[IDS[i]] 	= 0.0;
	}


}

segment::segment(){
}


//================================================================================================
//interval tree code

node::node(){};

node::node(vector<segment * > segments ){
	center 	= (double(segments[0]->start)  + double(segments[segments.size()-1]->stop)) / 2.;
	vector<segment * > Left;
	vector<segment * > Right;
	left=NULL, right=NULL;
	for (int i = 0 ; i < segments.size(); i++){
		if (segments[i]->stop < center){
			Left.push_back(segments[i]);
		}
		else if (segments[i]->start > center){
			Right.push_back(segments[i]);
		}
		else{
			current.push_back(segments[i]);
		}
	}
	if (Left.size() > 0){
		left 	= new node(Left);
	}
	if (Right.size() > 0){
		right 	= new node(Right);
	}
}
void node::insert_coverage(double x, double y, string ID){
	for (int i = 0 ; i < current.size(); i++){
		if (x > current[i]->start and  x < current[i]->stop  ){
			current[i]->G[ID]+=y;
		}
	}	

	if (x >= center and right != NULL ){
		right->insert_coverage(x, y, ID);
	}
	if (x <= center and left !=NULL){
		left->insert_coverage(x, y,  ID);
	}
}
void node::searchInterval(int start, int stop, vector<int>& finds ){
	for (int i = 0 ; i < current.size(); i++){
		if (stop > current[i]->start and  start < current[i]->stop  ){
			finds.push_back(1);
		}
	}	
	if (start >= center and right != NULL ){
		right->searchInterval(start, stop, finds);
	}
	if (stop <= center and left !=NULL){
		left->searchInterval(start, stop, finds);
	}	
}

void node::retrieve_nodes(vector<segment*> & saves){
	for (int i = 0; i < current.size(); i++){
		saves.push_back(current[i]);
	}
	if (right!= NULL){
		right->retrieve_nodes(saves);
	}
	if (left != NULL){
		left->retrieve_nodes(saves);		
	}
}

//================================================================================================

vector<string> split_under_tab(string str, char delimiter) {
	vector<string> internal;
	stringstream ss(str); // Turn the string into a stream.
	string tok;

	while(getline(ss, tok, delimiter)) {
		internal.push_back(tok);
	}
	return internal;
}

vector<segment *> sort_segments(vector<segment *> x){
	bool changed 	= true;
	while (changed){
		changed = false;
		for (int i = 1 ; i < x.size(); i++){
			if (x[i-1]->start> x[i]->start){
				changed 		= true;
				segment * copy 	= x[i-1];
				x[i-1] 	= x[i];
				x[i] 	= copy;
			}
		}
	}
	return x;

}


map<string, int> load::insert_bedgraph_data(map<string, node> A, vector<string> bedgraph_files, vector<string> IDS){
	map<string, int> NS;
	for (int i = 0 ; i < bedgraph_files.size(); i++){
		NS[IDS[i]] 	= 0.0;
	}
	#pragma omp parallel for
	for (int i = 0 ; i < bedgraph_files.size(); i++){
		string file 	= bedgraph_files[i], ID 	= IDS[i];
		printf("%s\n",file.c_str() );
		ifstream 	FH(file);
		string line,chrom; 
		vector<string> line_array;
		double x,y;
		int t 			= 0;
		string prevchrom 	= "";
		int j 			= 0, k =0, N 	= 0;
		vector<segment *> intervals;

		if (FH){
			while (getline(FH, line)){
				line_array 	= split_under_tab(line, '\t');
				if (line_array.size()==4){
					chrom  		= line_array[0];
					if (chrom!=prevchrom){
						j=0,N 	= 0;
						if (A.find(chrom)!=A.end()){
							N 	= 1;
						}

					}
					if (N>0){
						x 			= (stoi(line_array[2]) + stoi(line_array[1])) / 2.;
						y 			= (stoi(line_array[2]) - stoi(line_array[1]))*abs(stoi(line_array[3]));
						
						NS[IDS[i]] +=y;
						A[chrom].insert_coverage(x,y, ID);
					}
					//if (t > 100000){
					//	break;
					//}
					t+=1;
					prevchrom=chrom;

				
				}
			}			


		}else{
			printf("Couldn't Open: %s\n", file.c_str() );

		}
	}
	return NS;
}


void load::write_out_inserted_bedgaph_data(map<string, node> A, map<string, int> NS, vector<string> IDS, string out){
	ofstream FHW(out);
	typedef map<string, node>::iterator it_type;
	typedef vector<string>::iterator string_it;
	string line 	= "#";
	string line2 	= "#";
	for (string_it ID =IDS.begin(); ID!=IDS.end(); ID++){
		line+=(*ID)+",";
		line2+=to_string(NS[*ID])+",";
		
	}
	line=line.substr(0,line.size()-1) + "\n";
	line2=line2.substr(0,line2.size()-1) + "\n";


	FHW<<line;
	FHW<<line2;
	for (it_type a = A.begin(); a!=A.end(); a++){
		vector<segment *> vals;
		a->second.retrieve_nodes(vals);
		for (int i = 0 ; i < vals.size(); i++){
			int S 	= 0;
			string info 	= "";
			for (string_it ID =IDS.begin(); ID!=IDS.end(); ID++){
				S+=vals[i]->G[*ID];
				info+=to_string(vals[i]->G[*ID])+",";
			}
			if (S){
				line 	= "";
				line 	= vals[i]->chrom+"\t" + to_string(vals[i]->start)  + "\t" + to_string(vals[i]->stop) + "\t";
				info 	= info.substr(0,info.size()-1);
				FHW<<line<<info<<endl;
			}
		}
	}

}



map<string, node> load::make_interval_tree(string FILE, vector<string> IDS){
	ifstream FH(FILE);
	string line,chrom;
	int start, stop;

	vector<string> line_array;
	map<string, vector<segment*>> A;
	if (FH){
		while (getline(FH, line)){
			line_array 	= split_under_tab(line, '\t');
			if (line_array.size() > 2){
				chrom 	= line_array[0] , start = stoi(line_array[1]) , stop = stoi(line_array[2]);
				segment * S 	= new segment(chrom, start, stop, IDS);
				A[S->chrom].push_back(S);

			}	
		}
	}else{
		printf("Couldn't Open: %s\n", FILE.c_str() );
	}
	map<string, node> NT;
	typedef map<string, vector<segment *> >::iterator it_type_5;
	for(it_type_5 c = A.begin(); c != A.end(); c++) {
		NT[c->first] 	= node(c->second);
	}
	return NT;
}


void load::get_bg_files(string dir, vector<string> & bg_files, vector<string> & IDS){
	struct dirent *entry;
    DIR *dp;

    dp = opendir(dir.c_str());
    if (dp == NULL) {
        perror("opendir: Path does not exist or could not be read.");
    }
    int t = 0;
    cout.flush();
    vector<string> file_array;
    int i 	= 1;
    while ((entry = readdir(dp))){
    	string ID 	= "";
    	file_array 	= split_under_tab(entry->d_name, '.');
    	if (file_array.size() > 1){
    		ID 		= file_array[0];
    	}
    	if (!ID.empty()){
    		bg_files.push_back(dir+entry->d_name);
    		IDS.push_back(ID+"_"+to_string(i));
    		i++;
    	}
    }
}
vector<double> transform_line_array_to_double(vector<string> line_array){
	vector<double> x;
	for (int i = 0 ; i < line_array.size(); i++){
		x.push_back(stod(line_array[i]));
	}
	return x;

}


vector<vector<double>> load::coverage_stats_file(string FILE, vector<double> & counts2,vector<int> & IDS, vector<double> & centers, vector<string> & chroms ){
	ifstream FH(FILE);
	string line;
	vector<string> line_array;
	vector<string> tab_array;
	vector<double> counts;
	int ct =1;
	int ID =0;
	vector<vector<double>> D;
	while (getline(FH, line)){
		if (line.substr(0,1)=="#"){
			if (ct ==2){
				counts 	= transform_line_array_to_double(split_under_tab(line.substr(1, line.size()), ','));
			}
			ct+=1;
		}else{
			tab_array 			= split_under_tab(line , '\t');

			double center 		= (stod(tab_array[1]) + stod(tab_array[2])) / 2.;


			vector<double> x    = transform_line_array_to_double(split_under_tab(tab_array[3], ','));
			D.push_back(x);
			IDS.push_back(ID);
			chroms.push_back(tab_array[0]);
			centers.push_back(center);
			if (D.size()>1500){
				break;
			}
			ID+=1;

		}
	}
	counts2 	= counts;

	FH.close();
	return D;
}








vector<double> sort(vector<double> x){
	bool changed 	= true;
	while (changed){
		changed = false;
		for (int i = 1 ; i < x.size(); i++){
			if (x[i-1]> x[i]){
				changed 		= true;
				double copy 	= x[i-1];
				x[i-1] 	= x[i];
				x[i] 	= copy;
			}
		}
	}
	return x;
}

void intert_bed(string FILE,map<string, vector<double> > & G){
	ifstream FH(FILE);
	string line;
	vector<string> line_array;
	string chrom;
	double x;
	if (FH){
		while (getline(FH, line)){
			if (line.substr(0,1)!="#"){
				line_array 	= splitter2(line, "\t");
				if (line_array.size() > 2){
					chrom 		= line_array[0];
					x 			= (stod(line_array[1]) + stod(line_array[2]))/2.;
					G[chrom].push_back(x);
				}
			}
		}
	}else{
		printf("Couldn't Open: %s\n", FILE.c_str() );
	}
}

bool check(vector<double>x ){
	for (int i = 1; i <x.size(); i++){
		if (x[i-1] > x[i]){
			return false;
		}
	}
	return true;
}

int get_min(vector<double> x, int & idx, double & MIN,double threshold){
	MIN=pow(10,1000);
	for (int i = 0; i < x.size()-1; i++){
		if ((x[i+1] - x[i]) < MIN) {
			MIN 	= (x[i+1] - x[i]), idx=i;
		}
	}
	return MIN;
}

vector<double> join(vector<double> x,double distance,vector<vector<int>> & stats){
	bool joined =true;
	int idx;
	double MIN;
	double c;
	int s=1;
	for (int i = 0; i < x.size();i++){
		vector<int> 	row = {int(x[i]),int(x[i]),s};
		stats.push_back(row);
	}


	while (joined and x.size()>1){
		joined=false;
		get_min(x, idx, MIN, distance);
		if (MIN < distance){
			c 		= (x[idx] + 	x[idx+1])/2.;
			x[idx] 	= c;	
			int start 	= min(stats[idx][0],stats[idx+1][0]);
			int stop 	= max(stats[idx][1],stats[idx+1][1]);
			stats[idx][0]=start,stats[idx][1]=stop, stats[idx][2]+=stats[idx+1][2];
			
			x.erase(x.begin()+idx+1);
			stats.erase(stats.begin()+idx+1);
			
			joined=true;
		}
	}
	return x;

}

vector<double> join_fast(vector<double> x,double distance, vector<vector<int>> & stats){
	bool joined =true;
	int i 	= 0;
	vector<double> nx;
	int s 	= 0;
	while (i < x.size()){
		double c = x[i];

		vector<int> 	row = {int(x[i]),int(x[i]),s};
		while (i < x.size() and (x[i]-c)<distance ){
			c 		= (c+x[i])/2.;
			row[0] 	= min(int(x[i]),row[0] ), row[1] 	= max(int(x[i]),row[1]);
			row[2]++;

			i++;
		}
		stats.push_back(row);
		nx.push_back(c);
	}
	return nx;

}

map<string, vector<double> > load::get_bed_positions(string dir,string out,int distance,int fast, int window){
	map<string, vector<double> > G;
	struct dirent *entry;
    DIR *dp;

    dp = opendir(dir.c_str());
    if (dp == NULL) {
        perror("opendir: Path does not exist or could not be read.");
    }
    int t = 0;
    printf("loading....................");
    cout.flush();
    while ((entry = readdir(dp))){
    	string file_name 	= dir + entry->d_name;
    	intert_bed(file_name, G);
		// if (t > 10){
		// 	break;
		// }
    	t++;
    }
    printf("done\n");
    cout.flush();
    printf("sorting and joining........");
    cout.flush();
    typedef map<string, vector<double> >::iterator it_type;
    vector<string> chroms;
    for (it_type c = G.begin(); c!=G.end();c++){
    	chroms.push_back(c->first);
    }
    map<string, vector<vector<int> >> S;
    #pragma omp parallel for
    for (int c = 0; c<chroms.size();c++){

    	string chrom 	= chroms[c];
    	vector<vector<int>> stats;
    	S[chrom] 		= stats;
    	sort (G[chrom].begin(), G[chrom].end());
    	if (fast){
    		G[chrom] 		= join_fast(G[chrom],distance, S[chrom]);
    	}else{
	    	G[chrom] 		= join(G[chrom],distance,S[chrom]);
	    }
    	
    }
    printf("done\n");
    cout.flush();
    ofstream FHW(out);
    for (it_type c=G.begin(); c!=G.end();c++){
    	for (int i = 0 ; i < c->second.size();i++){
    		int start 	= S[c->first][i][0], stop 	= S[c->first][i][1];
    		if (start > 0){
    			if (window){
    	    		double x 	= (start + stop	)/2.;
    	    		start 		= x-window, stop 	= x+window;
		    	}
	    		FHW<<(c->first+"\t" + to_string(start) + "\t" + to_string(stop) + "\t" + to_string(S[c->first][i][2])+ "\n" );

	    	}
    	}
    }


    return G;
}












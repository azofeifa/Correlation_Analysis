#ifndef load_H
#define load_H
#include <string>
#include <vector>
#include <map>
#include "read_in_parameters.h"

class segment{
public:
	string chrom; 
	int start, stop;
	map<string, int> G;
	segment();
	segment(string, int , int,vector<string>);

};

class node{
public:
	double center;
	int start, stop;
	node * left;
	node * right;
	vector<segment * > current;
	void retrieve_nodes(vector<segment * >&);
	void insert_coverage(double, double, string);
	node();
	node(vector<segment *>);
	void searchInterval(int, int, vector<int> &) ;
};

class segment_fits{
public:
	string chrom;
	int start, stop, TSS;
	double N;
	double N_pos, N_neg;
	map<int, double> M;
	map<int, string> parameters;
	int model;
	double BIC_ratio;
	string ID;
	segment_fits();
	segment_fits(string, int, int, double, double, string);
	void get_model(double);
	string write();
};
namespace load{
	void get_bg_files(string , vector<string> & , vector<string> & );

	map<string, node> make_interval_tree(string,vector<string>);
	map<string, int> insert_bedgraph_data(map<string, node> , vector<string> , vector<string> );
	void write_out_inserted_bedgaph_data(map<string, node> , map<string, int>, vector<string>  , string );

	map<string, vector<double> > get_bed_positions(string,string,int,int,int);
}

#endif

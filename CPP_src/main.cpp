#include "read_in_parameters.h"
#include "load.h"
#include <omp.h>
using namespace std;

int main(int argc,char* argv[]){
    int threads  	= omp_get_max_threads();

	params * P = new params();
	fill_in_options(argv, P, 0);
	if (P->EXIT){
		printf("exiting...\n");
		return 0;
	}
	P->display();
	if (P->module=="join"){
		string bed_directory 	= P->p["-bed_dir"];
		string out_file		 	= P->p["-o"];
		int distance 			= stoi(P->p["-distance"]);
		int fast 				= stoi(P->p["-fast"]);
		int window 				= stoi(P->p["-window"]);
		
		load::get_bed_positions(bed_directory,out_file,distance,fast, window);
	}

	if (P->module=="insert"){
		string bed_file 		= P->p["-bed"];
		string bg_dir 			= P->p["-bg_dir"];
		string out_file 		= P->p["-o"];
		vector<string> bg_files, IDS;
		load::get_bg_files(bg_dir, bg_files, IDS);
		map<string, node> 	IT 	= load::make_interval_tree(bed_file,IDS);

		map<string, int> NS 	= load::insert_bedgraph_data(IT,  bg_files,  IDS);
		load::write_out_inserted_bedgaph_data(IT, NS, IDS,  out_file);


	}





	return 0;
}
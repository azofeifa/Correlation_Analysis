#include "read_in_parameters.h"
#include <iostream>
using namespace std;
params::params(){
	p["-bed_dir"] 	= "";
	p["-distance"] 	= "1000";
	p["-o"] 		= "";
	p["-bed"] 		= "";
	p["-bg_dir"] 	= "";
	p["-fast"] 		= "0";
	p["-norm"] 		= "0";
	p["-window"] 	= "0";
	p["-test"] 		= "0";
	p["-threshold"] = "0.7";
	module 		= "";
	EXIT 		= false;

}
bool params::check(){
	typedef map<string, string >::iterator it_type; 
	bool EXIT = false;
	for (it_type i = p.begin(); i!= p.end(); i++){
		if (i->second.empty()){
			printf("%s was not specified\n",i->first.c_str() );
			EXIT=true;
		}
	}
	return EXIT;
}
void params::help(){
}
const std::string currentDateTime2() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%m/%d/%Y %X", &tstruct);

    return buf;
}

void params::display(){
	string header 	= "";
	header  += "===============================================\n";

	if (module=="join"){
		header 	+= "...joining bed file intervals...\n";
		header  += "Date Time    : "+currentDateTime2()+"\n";
		header 	+= "-bed_dir     : "+p["-bed_dir"]+"\n";
		header 	+= "-distance    : "+p["-distance"]+"\n";
		header 	+= "-fast        : "+p["-fast"]+"\n";
		header 	+= "-window      : "+p["-window"]+"\n";
		
	}
	if (module=="insert"){
		header 	+= "....inserting coverage data....\n";
		header  += "Date Time    : "+currentDateTime2()+"\n";
		header 	+= "-bed_dir     : "+p["-bg_dir"]+"\n";
		header 	+= "-bed         : "+p["-bed"]+"\n";

	}
	if (module=="correlate"){
		header 	+= "....computing correlation coefficients....\n";
		header  += "Date Time    : "+currentDateTime2()+"\n";
		header 	+= "-bed         : "+p["-bed"]+"\n";
		header	+= "-norm        : "+p["-norm"]+"\n";		
	}

	header 	+= "-o           : "+p["-o"]+"\n";
	header  += "===============================================\n";
	printf("%s\n", header.c_str() );	

}




void fill_in_options(char* argv[],params * P, int rank){
	string F 		= "";
	char * COM 		= "-";
	string current 	= "";
	argv++;
	if (*argv){
		P->module  = string(*argv);
		if (not (P->module== "join" or P->module=="insert"  or P->module=="correlate"  )){
			printf("did not specify module, either DB or EVAL\ntry -h or --help for quick reference software usage\n");
			P->EXIT 	= true;
		}
	}else{
		printf("No module specified\ntry -h or --help for quick reference software usage\n");
		P->EXIT 		= true;
	}
	while (*argv){
		F 	= string(*argv); 
		if(!current.empty()){
			P->p[current] 	= F;
			current 		= "";
		}
		if ((*argv)[0] == COM[0]){
			current 	= F;
			if ( P->p.find(F) ==P->p.end() ){
				if (rank == 0 and F.substr(0,2)!="-h" and F.substr(0,7)!="--help" ){
					printf("Unknown user option: %s\n", F.c_str() );
				}
				current 	= "";
			}
			else if(P->p.find(F) !=P->p.end() ){
				current 	= F;
			}
		}
		if (F.substr(0,2) == "-h" or F.substr(0,7)=="--help" ){
			P->help();
			P->EXIT;
		}
		
		argv++;
	}

	
	if (!current.empty()){
		P->p[current] 	= F;
	}

}





















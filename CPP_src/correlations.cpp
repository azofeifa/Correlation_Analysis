#include "correlations.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <omp.h>
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
using namespace std;
void make_weigts(double * weights, vector<double> counts, int m){
	double S=0;
	for (int i = 0 ; i < m; i++ ){
		S+=counts[i];
	}
	for (int i = 0 ; i < m; i++ ){
		weights[i] 	= counts[i]/S;
	}
}
void fill_mean_std(double ** D, double * mean, double * std, double * weights, int n, int m){
	for (int i = 0 ; i < n; i++){
		double S 	= 0;
		for (int j = 0 ; j < m; j++){
			S+=D[i][j]*weights[j];
		}
		mean[i] 	= S;
		S 			= 0;
		for (int j = 0 ; j < m; j++){
			S+=(pow(D[i][j]-mean[i], 2))*weights[j];
		}
		std[i] 		= sqrt(S);
	}
}

double  covariance(double * x, double * y, int m, double * weights, double mx, double my){
	double S 	 = 0;
	for (int j =0; j < m; j++){
		S+=weights[j]*(x[j]-mx)*(y[j] -my  ) ;
	}
	return S;
}
double LOG(double x){
	if (not x){
		return 0.0;
	}
	return log10(x);
}


void correlate::compute(vector<vector<double>> A, vector<double> counts, vector<int> IDS, vector<string> chroms, vector<double> centers, string out_file){
	int n 	= A.size(), m = counts.size();
	double ** D 		= new double *[n];
	double ** P 		= new double *[n];
	double * mean 		= new double[n];
	double * std 		= new double[n];
	double * weights 	= new double[m];

	make_weigts(weights, counts, m);

	for (int i = 0 ; i < n;i++){
		D[i] 	= new double[m];
		P[i] 	= new double[n];
		for (int j = 0 ; j < n;j++){
			P[i][j] 	= 0.0;
		}
		for (int j = 0 ; j < m; j++){
			D[i][j] 	= LOG(A[i][j]/counts[j]) ;
		}
	}
	double pearons;
	fill_mean_std(D, mean, std, weights, n, m);

	int threads  	= omp_get_max_threads();
	int cts 		= n / threads;
//	threads 		= 1;
	double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
	#pragma omp parallel num_threads(threads)
	{
		int tid 	= omp_get_thread_num();
		int start 	= tid*cts;
		int stop 	= (tid+1)*(cts);
		if (tid+1 == threads){
			stop 	= n;
		}
		for (int i = start ; i < stop;i++){
			P[i][i] 		= 1.0;
			for (int j = i +1; j < n ; j++){
				pearons 	= covariance(D[i], D[j], m, weights,  mean[i], mean[j] )/(std[i]*std[j]);
				P[i][j] 	= pearons, P[j][i] 	= pearons;
			}
		}
	}
	double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    // cout<<endl;
    // cout<< "Wall Time = " << wall1 - wall0 << endl;
    // cout<< "CPU Time  = " << cpu1  - cpu0  << endl;

    
	ofstream FHW(out_file);
	string line 	= "#Locations\t";
	string line2 	= "#MeanExpression\t";
	string line3 	= "#VarExpression\t";
		
	for (int i = 0 ; i < n ; i++){
		line += chroms[i]+":" + to_string(int(centers[i]));
		line2 += to_string(mean[i]);
		line3 += to_string(std[i]);
		if (i+1 < n){
			line+=",";
			line2+=",";
			line3+=",";
		}
		else{
			line+="\n";
			line2+="\n";
			line3+="n";
		}
	}
	FHW<<line;
	FHW<<line2;
	FHW<<line3;
	for (int i = 0 ; i < n; i++){
		line 	= "";
		for (int j = i+1; j < n ; j++){
			int dist 	= int((abs(centers[i]-centers[j]))/1000.0);
			if (chroms[i]!=chroms[j]){
				dist 	= -1;
			}

			if (P[i][j]>0.7){
				line+=  to_string(IDS[i]) + ","+ to_string(IDS[j]) + "\t" + to_string(P[i][j]) + "," + to_string(dist) + "\n";
			}
		}
		FHW<<line;			
		
	}

	FHW.close();


	delete [] mean, std, weights;

	for (int i = 0 ; i < n;i++){
		delete [] D[i];
		delete [] P[i];
	}
	delete [] D;
	delete [] P;



}



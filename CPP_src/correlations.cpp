#include "correlations.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <omp.h>
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
		S 	= 0;
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
			D[i][j] 	= (A[i][j] / counts[j])*pow(10,3);
		}
	}
	double pearons;
	fill_mean_std(D, mean, std, weights, n, m);

//	#pragma omp parallel for
	for (int i = 0 ; i < n;i++){
		P[i][i] 		= 1.0;
		for (int j = i +1; j < n ; j++){
			pearons 	= covariance(D[i], D[j], m, weights,  mean[i], mean[j] )/(std[i]*std[j]);
			P[i][j] 	= pearons, P[j][i] 	= pearons;
		}
	}

    
	ofstream FHW(out_file);
	string line 	= "";
		
	for (int i = 0 ; i < n ; i++){
		line 	= chroms[i]+":" + to_string(int(centers[i]));
		if (i+1 < n){
			FHW<<line+",";
		}
		else{
			FHW<<line<<endl;
		}

		
	}

	for (int i = 0 ; i < n; i++){
		line 	= "";
		for (int j = i+1; j < n ; j++){
			int dist 	= int(abs(centers[i]-centers[j]));
			if (chroms[i]!=chroms[j]){
				dist 	= -1;
			}
			line+=  to_string(IDS[i]) + ","+ to_string(IDS[j]) + "\t" + to_string(P[i][j]) + "," + to_string(dist) + "\n";
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



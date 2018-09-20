#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <cmath>

using namespace std;

int bipartition(int dim, int i0, int iN, double **data, int cluster_start1, int cluster_start2, int cluster_size1, int cluster_size2, 
	double *cluster_bdry1_min, double *cluster_bdry1_max, double *cluster_bdry2_min, double *cluster_bdry2_max, double *cluster_centroid1, 
	double *cluster_centroid2, int *cluster_assign, int cluster1, int cluster2) {

	int max_iterations = ((iN - i0) * 2);
	int i = i0;
	int count = 0;
	cout << "New call " << cluster1 << " " << cluster2 << " " << i0 << " " << iN << endl;
	while (i < max_iterations) {
		if (count == 0 && i == iN) {
			iN = (iN - i0) + iN;
		}
		cout << i << endl;
		i++;
	}
	cout << "----End of call----" << endl;
	return 0;
}

int kdtree(int dim, int ndata, double **data, int kk, int *cluster_start, int *cluster_size, double **cluster_bdry, 
	double **cluster_centroid, int *cluster_assign) {
	
	int i = 0;
	int j = -1;
	int skip = 0;
	int temp1, temp2 = 0;
	int num = 0;
	int *temp = (int*)calloc(1, sizeof(int));
	int depth = pow(2, skip);
	int kd_size = (kk - 1) * 2;
	int* kdtree = (int*)calloc(kd_size, sizeof(int));
	int *size = (int*)calloc(kd_size, sizeof(int));
	int x = ndata;
	double remainder = 0;

	while (i < kd_size) {
		if (j == kk) {
			skip++;
			depth = pow(2, skip);
			j = 0 - depth;
		}
		kdtree[i] = j + depth;
		i++;

		cout << kdtree[i] << endl;
	}
/*
	i = 0;
	j = pow(2, skip);

	while (i < kd_size) {
		if (kdtree[i] == 0 and i != 0) {
			skip--;
			j = pow(2, skip);
		}
		if (x % j != 0) {
			remainder = ndata % 2;
			x = x - remainder;
		}
		if (remainder != 0) {
			size[i] = (ndata / j) + 1;
			remainder--;
		}
		else {
			size[i] = (ndata) / j;
		}

		cout << size[i] << endl;

		i++;
	}
*/
	delete [] kdtree;

	return 0;
}

int main()
{
	random_device rd;

	mt19937 randomGenerator(rd());

	uniform_int_distribution<int> distributionDim(2, 3);
	uniform_int_distribution<int> distributionNdata(2, 10);
	uniform_real_distribution<double> distributionData(1, 10);
	uniform_int_distribution<int> distributionKk(1, 4);

	int i, j = 0;
	int dim = 4;
	int ndata = 24;
	int kk = 8;
	int *cluster_start = (int*)calloc(kk, sizeof(int));
	int *cluster_size = (int*)calloc(kk, sizeof(int));
	int *cluster_assign = (int*)calloc(ndata, sizeof(int));
	double *sum = (double*)calloc(dim, sizeof(double));
	double **data = (double**)calloc(ndata, sizeof(double));		//suggestion: have ndata*dim
	double **cluster_bdry = (double**)calloc(kk, sizeof(double));
	double **cluster_centroid = (double**)calloc(kk, sizeof(double));

	for (i = 0; i < ndata; i++) {
		data[i] = (double*)calloc(dim, sizeof(double));
	}

	for (i = 0; i < kk; i++) {
		cluster_bdry[i] = (double*)calloc(dim*2, sizeof(double));
		cluster_centroid[i] = (double *)calloc(dim, sizeof(double));
	}

	cluster_size[0] = ndata;
	cout << "------------------------------" << endl;
	cout << "Data Set" << endl;
	for (i = 0; i < ndata; i++) {
		cout << i << "- ";
		for (j = 0; j < dim; j++) {
			data[i][j] = distributionData(randomGenerator);
			cout << data[i][j] << " ";
			
			sum[j] = sum[j] + data[i][j];

			if (i == 0) {
				cluster_bdry[0][j * 2] = *&data[i][j];
				cluster_bdry[0][j * 2 + 1] = *&data[i][j];
			}
			else if (cluster_bdry[0][j * 2] > *&data[i][j]) {
				cluster_bdry[0][j * 2] = *&data[i][j];
			}
			else if (cluster_bdry[0][j * 2 + 1] < *&data[i][j]) {
				cluster_bdry[0][j * 2 + 1] = *&data[i][j];
			}
		}
		cout << endl;
	}

	cout << "------------------------------" << endl;

	cout << "Cluster Boundaries" << endl;
	for (i = 0; i < dim*2; i++) {
		cout << cluster_bdry[0][i] << " ";
	}
	cout << endl;
	cout << "------------------------------" << endl;

	cout << "Cluster Centroids" << endl;
	for (i = 0; i < dim; i++) {
		cluster_centroid[0][i] = sum[i] / ndata;
		cout << cluster_centroid[0][i] << " ";
	}
	cout << endl;
	cout << "------------------------------" << endl;

	cout << "Dimensions: " << dim << endl;
	cout << "Ndata: " << ndata << endl;
	cout << "Kk: " << kk << endl;

	kdtree(dim, ndata, data, kk, cluster_start, cluster_size, cluster_bdry, cluster_centroid, cluster_assign);

	getchar();

}

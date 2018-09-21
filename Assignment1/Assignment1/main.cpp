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
	int j = 0;
	int skip = 0;
	int temp1, temp2 = 0;
	int grind1, grind2 = 0;
	int num = ndata;
	int depth = pow(2, skip);
	int kd_size = (kk - 1) * 2;
	int *kdtree = (int*)calloc(kd_size, sizeof(int));
	int *size = (int*)calloc(kd_size, sizeof(int));
	int x = ndata;
	double remainder = 0;

	while (i < kd_size) {
		if (j == kk) {
			j = 0;
			skip++;
			depth = pow(2, skip);
		}
		kdtree[i] = j;
		cout << kdtree[i] << endl;
		i++;
		j = j + depth;
	}

	i = 0;
	skip = skip + 2;

	while (i < kd_size) {
		if (kdtree[i] == 0) {
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

	for (i = kd_size - 1; i >= 0; i = i - 2) {
		temp1 = kdtree[i];
		temp2 = kdtree[i - 1];
		cluster_size[temp1] = size[i];
		cluster_size[temp2] = size[i - 1];

		if (num == 0 && temp2 == 0){
			grind1 = 0;
			grind2 = cluster_size[temp2] + cluster_size[temp1];
		}
		else {
			grind1 = num;
			num = num - cluster_size[temp2] - cluster_size[temp1];
			grind2 = num;
		}

		cout << grind1 << grind2 << endl;

		for (j = 0; j < dim; j++) {
			bipartition(j, grind2, grind1, data, cluster_start[temp2], cluster_start[temp1], cluster_size[temp2], cluster_size[temp1],
				&cluster_bdry[temp2][j], &cluster_bdry[temp2][j + 1], &cluster_bdry[temp1][j], &cluster_bdry[temp1][j + 1], &cluster_centroid[temp2][j],
				&cluster_centroid[temp1][j], cluster_assign, temp2, temp1);
		}

		if (num == 0) {
			num = ndata;
		}

	}

	delete[] kdtree;
	delete[] size;

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
	double **data = (double**)calloc(ndata, sizeof(double));
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
				cluster_bdry[0][j * 2] = data[i][j];
				cluster_bdry[0][j * 2 + 1] = data[i][j];
			}
			else if (cluster_bdry[0][j * 2] > data[i][j]) {
				cluster_bdry[0][j * 2] = data[i][j];
			}
			else if (cluster_bdry[0][j * 2 + 1] < data[i][j]) {
				cluster_bdry[0][j * 2 + 1] = data[i][j];
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

	delete[] cluster_start;
	delete[] cluster_size;
	delete[] sum;
	delete[] cluster_bdry;
	delete[] cluster_assign;
	delete[] data;

	getchar();

}

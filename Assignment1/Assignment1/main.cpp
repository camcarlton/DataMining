#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>

/* Let data point  (x1, x2, ... xd)
   A cluster's bdry is (m1, M1, m2, M2, . . ., md, Md) */
 
using namespace std;

void print(int dim, int ndata, double **data, int kk, int *cluster_start, int *cluster_size, double **cluster_bdry,
	double **cluster_centroid, int *cluster_assign) {

  	int i, j = 0;

	for (i = 0; i < kk; i++) {
		cout << "--Cluster " << i << "--" << endl;
		cout << "Start " << cluster_start[i] << endl;
		cout << "Size " << cluster_size[i] << endl;

		for (j = 0; j < dim; j++) {
			cout << "----Dim " << j << "----" << endl;
			cout << "Min " << cluster_bdry[i][j * 2] << endl;
			cout << "Max " << cluster_bdry[i][j * 2 + 1] << endl;
			cout << "Centroid " << cluster_centroid[i][j] << endl;
		}
	}

	/*int i = 0;

	for (i = 0; i < ndata; i++) {
		cout << "cluster_assign slot " << i << " holds " << cluster_assign[i] << endl;
	}*/
}

double calculateMean(int dim, int i0, int iN, double **data, double *mean) {
	int i, j = 0;

	for (i = i0; i < iN; i++) {
		for (j = 0; j < dim; j++) {
			mean[j] = mean[j] + data[i][j];
		}
	}

	for (i = 0; i < dim; i++) {
		mean[i] = mean[i] / iN;
	}

	return *mean;
}

double calculateVariance(int dim, int i0, int iN, double *mean, double **data, double *variance) {
	int i, j = 0;
	double new_data = 0;

	for (i = i0; i < iN; i++) {
		for (j = 0; j < dim; j++) {
			variance[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
		}
	
	}

	for (i = 0; i < dim; i++) {
		variance[i] = variance[i] / iN;
	}

	return *variance;
}

int maxVariance(int dim, double *variance) {
	int i = 0;
	double max = 0;
	int max_dim = -1;

	for (i = 0; i < dim; i++) {
		if (i == 0) {
			max = variance[i];
			max_dim = i;
		}
		else if (max < variance[i]) {
			max = variance[i];
			max_dim = i;
		}
	}

	return max_dim;
}

double sortAcsending(int dim, int i0, int iN, int max_dim, double **data) {
	int i = 0;
	int j = 0;
	int count = 0;

	while (i < iN-1) {
		if (data[i][max_dim] > data[i + 1][max_dim]) {
			for (j = 0; j < dim; j++) {
				int temp = data[i][j];
				data[i][j] = data[i + 1][j];
				data[i + 1][j] = temp;
			}
			count++;
		}
		if (i == iN-2 && count != 0) {
			i = 0;
			count = 0;
		}
		else {
			i++;
		}
	}

	return **data;
}

int assignCluster(int i0, int iN, int max_dim, double mean, double **data, int cluster1, int cluster2, int *cluster_assign) {
	int i = 0;

	for (i = i0; i < iN; i++) {
		if (data[i][max_dim] > mean) {
			cluster_assign[i] = cluster2;
		}
		else if (data[i][max_dim] < mean) {
			cluster_assign[i] = cluster1;
		}

		//cout << "i = " << cluster_assign[i] << endl;
		
	}
	cout << endl;
	return *cluster_assign;
}

int defineCluster(int i0, int iN, int *cluster_assign, int cluster1, int cluster2, int *cluster_size, int *cluster_start) {
	int i = 0;
	cluster_size[cluster1] = 0;
	cluster_size[cluster2] = 0;

	for (i = i0; i < iN; i++) {
		if (cluster_assign[i] == cluster1 and cluster_size[cluster1] == 0) {
			cluster_size[cluster1]++;
			cluster_start[cluster1] = i;
		}
		else if (cluster_assign[i] == cluster1) {
			cluster_size[cluster1]++;
		}
		else if (cluster_assign[i] == cluster2 and cluster_size[cluster2] == 0) {
			cluster_size[cluster2]++;
			cluster_start[cluster2] = i;
		}
		else if (cluster_assign[i] == cluster2) {
			cluster_size[cluster2]++;
		}
	}

	return *cluster_size, *cluster_start;
}

double clusterBoundry(int dim, int i0, int iN, int *cluster_assign, int cluster1, int cluster2, double **data, double **cluster_bdry) {
	int i = 0;
	int j = 0;

	for (i = 0; i < dim; i++) {
		cluster_bdry[cluster1][i] = data[i0][i];
		cluster_bdry[cluster1][i + 1] = data[i0][i];

		for (j = i0; j < iN; j++) {
			if (cluster_assign[j] == cluster1) {
				if (cluster_bdry[cluster1][i] > data[j][i]) {
					cluster_bdry[cluster1][i] = data[j][i];
				}
				else if (cluster_bdry[cluster1][i + 1] < data[j][i]) {
					cluster_bdry[cluster1][i + 1] = data[j][i];
				}
			}
			else if (cluster_assign[j] == cluster2) {
				if (cluster_bdry[cluster2][i] > data[j][i]) {
					cluster_bdry[cluster2][i] = data[j][i];
				}
				else if (cluster_bdry[cluster2][i + 1] < data[j][i]) {
					cluster_bdry[cluster2][i + 1] = data[j][i];
				}
			}
		}
	}

	return **cluster_bdry;
}

int kdtree(int dim, int ndata, double **data, int kk, int *cluster_start, int *cluster_size, double **cluster_bdry, 
	double **cluster_centroid, int *cluster_assign) {
	
	int i = 0;
	int j = 0;
	int max_dim = 0;
	int depth = 0;
	int temp1, temp2 = 0;
	int start, end = 0;
	int sizeTotal = 0;
	int stride = pow(2, depth);
	int kd_size = (kk - 1) * 2;
	int *kdtree = (int*)calloc(kd_size, sizeof(int));
	double *mean = (double*)calloc(dim, sizeof(int));
	double *variance = (double*)calloc(dim, sizeof(int));

	while (i < kd_size) {
		if (j == kk) {
			j = 0;
			depth++;
			stride = pow(2, depth);
		}
		kdtree[i] = j;
		i++;
		j = j + stride;
	}

	for (i = kd_size - 1; i >= 0; i = i - 2) {
		temp1 = kdtree[i-1];
		temp2 = kdtree[i];
		sizeTotal = cluster_size[temp1];
		start = cluster_start[temp1];
		end = start + sizeTotal;

		*cluster_centroid[temp1] = calculateMean(dim, start, end, data, mean);
		*variance = calculateVariance(dim, start, end, cluster_centroid[temp1], data, variance);
		max_dim = maxVariance(dim, variance);
		**data = sortAcsending(dim, start, end, max_dim, data);
		*cluster_assign = assignCluster(start, end, max_dim, mean[max_dim], data, temp1, temp2, cluster_assign);
		*cluster_size, *cluster_start = defineCluster(start, end, cluster_assign, temp1, temp2, cluster_size, cluster_start);
		**cluster_bdry = clusterBoundry(dim, start, end, cluster_assign, temp1, temp2, data, cluster_bdry);
	}

	print(dim, ndata, data, kk, cluster_start, cluster_size, cluster_bdry, cluster_centroid, cluster_assign);

	delete[] mean;
	delete[] variance;

	return 0;
}

int main()
{
	random_device rd;

	mt19937 randomGenerator(rd());

	uniform_int_distribution<int> distributionDim(2, 3);
	uniform_int_distribution<int> distributionNdata(2, 10);
	uniform_real_distribution<double> distributionData(1, 100);
	uniform_int_distribution<int> distributionKk(1, 4);

	int i, j = 0;
	int dim = 4;
	int ndata = 100;
	int kk = 8;
	int *cluster_start = (int*)calloc(kk, sizeof(int));
	int *cluster_size = (int*)calloc(kk, sizeof(int));
	int *cluster_assign = (int*)calloc(ndata, sizeof(int));
	double *variance = (double*)calloc(dim, sizeof(double));
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

	for (i = 0; i < ndata; i++) {
		for (j = 0; j < dim; j++) {
			data[i][j] = distributionData(randomGenerator);

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
	}

	kdtree(dim, ndata, data, kk, cluster_start, cluster_size, cluster_bdry, cluster_centroid, cluster_assign);

	delete[] cluster_start;
	delete[] cluster_size;
	delete[] variance;
	delete[] cluster_bdry;
	delete[] cluster_assign;
	delete[] data;

	getchar();

}

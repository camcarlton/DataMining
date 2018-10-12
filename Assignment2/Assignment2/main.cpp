#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>

using namespace std;

//Possibly need to treat h differently
//Sort data before it is sent to LSH
//Cluster_hashval has the possibly of assign the same hashvalue for different clusters
//Assuming main with pass in the cluster number into LSH
//Vector math


int lsh(int dim, int ndata, int m, int *ncluster_ptr, int *cluster_start, int *cluster_size, int **cluster_hashval, double W, double **data, double **h){
	int x, y, z;
	double unitVector = 0;
	double min = -1;

	for (x = 0; x < m; x++) {
		for (y = 0; y < ndata; y++) {
			for (z = 0; z < dim; z++) {
				unitVector += pow(data[y][z] + h[x][z], 2);
			}
			unitVector = floor(sqrt(unitVector) / W);
			if (cluster_hashval[*ncluster_ptr][x] == 0 and cluster_start[*ncluster_ptr] == -1) {
				cluster_hashval[*ncluster_ptr][x] = unitVector;
				cluster_start[*ncluster_ptr] = y;
				*cluster_size += 1;
			}
			else if (cluster_hashval[*ncluster_ptr][x] == unitVector) {
				*cluster_size += 1;
			}
		}
	}

	return *cluster_size, *cluster_start, **cluster_hashval;
}

int main()
{
	random_device rd;

	mt19937 randomGenerator(rd());

	uniform_int_distribution<int> distributionDim(3, 8);
	uniform_int_distribution<int> distributionNdata(10, 10000);
	uniform_int_distribution<int> distributionClusters(1, 16);
	uniform_real_distribution<double> distributionData(1, 10);
	uniform_real_distribution<double> distributionW(0, 1);

	int i, j;
	int dim = distributionDim(randomGenerator);
	int ndata = distributionNdata(randomGenerator);
	int m = dim - 1;
	int ncluster = distributionClusters(randomGenerator);
	int *cluster_start = (int*)calloc(ncluster, sizeof(int));
	int *cluster_size = (int*)calloc(ncluster, sizeof(int));
	int **cluster_hashval = (int**)calloc(ncluster, sizeof(int));
	double W = distributionW(randomGenerator);
	double **data = (double**)calloc(ndata, sizeof(double));
	double **h = (double**)calloc(m, sizeof(double));

	for (i = 0; i < dim; i++) {
		data[i] = (double*)calloc(m, sizeof(double));
		h[i] = (double*)calloc(m, sizeof(double));
	}

	for (i = 0; i < ncluster; i++) {
		cluster_hashval[i] = (int*)calloc(m, sizeof(int));
	}

	for (i = 0; i < m; i++) {
		for (j = 0; j < dim; j++) {
			h[i][j] = distributionData(randomGenerator);
		}
	}

	for (i = 0; i < ndata; i++) {
		for (j = 0; j < dim; j++) {
			h[i][j] = distributionData(randomGenerator);
		}
	}

	for (i = 0; i < ncluster; i++) {
		lsh(dim, ndata, m, &i, cluster_start, cluster_size, cluster_hashval, W, data, h);
	}

	delete[] cluster_start;
	delete[] cluster_size;
	delete[] cluster_hashval;
	delete[] data;
	delete[] h;

	getchar();

}

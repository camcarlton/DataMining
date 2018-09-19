#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>

using namespace std;

int main()
{
	random_device rd;

	mt19937 randomGenerator(rd());

	uniform_int_distribution<int> distributionDim(2, 8);
	uniform_int_distribution<int> distributionNdata(10, 100);
	uniform_real_distribution<double> distributionData(1, 100);
	uniform_int_distribution<int> distributionKk(1, 4);

	int i, j = 0;
	int dim = distributionDim(randomGenerator);
	int ndata = distributionNdata(randomGenerator);
	int kk = 2 ^ distributionKk(randomGenerator);
	int *cluster_start = (int*)calloc(kk, sizeof(int));
	int *cluster_size = (int*)calloc(kk, sizeof(int));
	int *cluster_assign = (int*)calloc(ndata, sizeof(int));
	double *sum = (double*)calloc(dim, sizeof(double));
	double **data = (double**)calloc(ndata, sizeof(double));
	double **cluster_belry = (double**)calloc(kk, sizeof(double));
	double **cluster_centroid = (double**)calloc(kk, sizeof(double));

	for (i = 0; i < ndata; i++) {
		data[i] = (double*)calloc(dim, sizeof(double));
	}

	for (i = 0; i < kk; i++) {
		cluster_belry[i] = (double*)calloc(dim*2, sizeof(double));
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
				cluster_belry[0][j * 2] = *&data[i][j];
				cluster_belry[0][j * 2 + 1] = *&data[i][j];
			}
			else if (cluster_belry[0][j * 2] > data[i][j]) {
				cluster_belry[0][j * 2] = *&data[i][j];
			}
			else if (cluster_belry[0][j * 2 + 1] < data[i][j]) {
				cluster_belry[0][j * 2 + 1] = *&data[i][j];
			}
		}
		cout << endl;
	}

	cout << "------------------------------" << endl;

	cout << "Cluster Boundaries" << endl;
	for (i = 0; i < dim*2; i++) {
		cout << cluster_belry[0][i] << " ";
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

	free(cluster_start);
	free(cluster_size);
	free(cluster_assign);
	free(data);
	free(cluster_belry);
	free(cluster_centroid);

	getchar();

}

int kdtree(int dim, int ndata, double **data, int kk, int *cluster_start, int *cluster_size, double **cluster_belry, double **cluster_centroid, int *cluster_assign) {
	return 0;
}

int bipartition(int dim, int i0, int iN, double **data, int cluster_start[2], int cluster_size[2], double *cluster_belry[2], double *cluster_centroid[2], int *cluster_assign) {
	return 0;
}
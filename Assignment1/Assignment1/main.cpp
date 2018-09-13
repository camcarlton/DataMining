#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>

using namespace std;

int main()
{
	random_device rd;

	mt19937 randomGenerator(rd());

	uniform_int_distribution<int> distributionDim(1, 8);
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
	double **data = (double**)calloc(ndata, sizeof(double));
	double **cluster_belry = (double**)calloc(kk, sizeof(double));
	double **cluster_centroid = (double**)calloc(kk, sizeof(double));

	for (i = 0; i < ndata; i++) {
		data[i] = (double*)calloc(dim, sizeof(double));
	}

	for (i = 0; i < kk; i++) {
		cluster_belry[i] = (double*)calloc(dim * 2, sizeof(double));
		cluster_centroid[i] = (double *)calloc(dim, sizeof(double));
	}

	for (i = 0; i < ndata; i++) {
		cout << i << "- ";
		for (j = 0; j < dim; j++) {
			data[i][j] = distributionData(randomGenerator);
			cout << data[i][j] << " ";
		}
		cout << endl;
	}

	cout << "Dimensions: " << dim << endl;
	cout << "Ndata: " << ndata << endl;

	free(cluster_start);
	free(cluster_size);
	free(cluster_assign);
	free(data);
	free(cluster_belry);
	free(cluster_centroid);

	getchar();

}
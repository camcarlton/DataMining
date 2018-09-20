#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <cmath>

using namespace std;

int bipartition(int dim, int i0, int iN, double **data, int cluster_start1, int cluster_start2, int cluster_size1, int cluster_size2, 
	double *cluster_bdry1_min, double *cluster_bdry1_max, double *cluster_bdry2_min, double *cluster_bdry2_max, double *cluster_centroid1, 
	double *cluster_centroid2, int *cluster_assign, int cluster1, int cluster2) {

	int count = 0;
	int i = i0;
	double *sum = (double*)calloc(2, sizeof(double));
	double *temp;

	cluster_start1 = i0;
	cluster_size1 = iN - i0;
	*cluster_bdry1_min = data[i0][dim];
	*cluster_bdry1_max = data[i0][dim];

	while (count != 2) {
		if (count == 0 && i == iN) {
			//*cluster_centroid1 = sum[count] / cluster_size1;
			cluster_start2 = i;
			cluster_bdry2_min = &data[i][dim];
			cluster_bdry2_max = &data[i][dim];
			iN = (i - i0) + iN;
			count++;
		}
		else if (count == 1 && i == iN) {
			cluster_size2 = i-i0;
			//*cluster_centroid2 = sum[count] / cluster_size2;
			count++;
		}
		else if (count == 0 && i != iN) {
			cluster_assign[i] = cluster1;

			if (cluster_bdry1_min > &data[i][dim]) {
				cluster_bdry1_min = &data[i][dim];
			}
			else if (*cluster_bdry1_max < data[i][dim]) {
				cluster_bdry1_max = &data[i][dim];
			}
		}
		else if (count == 1 && i != iN) {
			cluster_assign[i] = cluster2;

			if (cluster_bdry2_min > &data[i][dim]) {
				cluster_bdry2_min = &data[i][dim];
			}
			else if (cluster_bdry2_max < &data[i][dim]) {
				cluster_bdry2_max = &data[i][dim];
			}
		}
		//sum[count] = sum[count] + data[i][dim];

		i++;
	}

	return 0;
}

int kdtree(int dim, int ndata, double **data, int kk, int *cluster_start, int *cluster_size, double **cluster_bdry, 
	double **cluster_centroid, int *cluster_assign) {
	
	int i = 0;
	int j = 0;
	int skip = 0;
	int temp1, temp2, x = 0;
	int num = 0;
	int depth = pow(2, skip);
	int kd_size = (kk - 1) * 2;
	int temp;
	int *kdtree	= (int*)calloc(kd_size, sizeof(int));

	while (i < kd_size) {
		if (j == kk) {
			j = 0;
			skip++;
			depth = pow(2, skip);
		}
		kdtree[i] = j;
		i++;
		j = j + depth;
	}

	depth = 0;

	for (i = kd_size - 1; i >= 0; i = i - 2) {
		temp1 = kdtree[i];
		temp2 = kdtree[i - 1];

		if (temp2 == 0) {
			depth++;
			num = 0;
			x = ndata / (pow(2, depth));
		}

	
		for (j = 0; j < dim; j++) {
			bipartition(j, num, num + x, data, cluster_start[temp2], cluster_start[temp1], cluster_size[temp2], cluster_size[temp1], 
				&cluster_bdry[temp2][j], &cluster_bdry[temp2][j + 1], &cluster_bdry[temp1][j], &cluster_bdry[temp1][j + 1], &cluster_centroid[temp2][j],
				&cluster_centroid[temp1][j], cluster_assign, temp2, temp1);
		}

		num = num + x;
	}

	delete kdtree;

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
	int ndata = 10;
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

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

void reassignment(int datapoint, int ndata, int dim, int old_cluster, int new_cluster, double *data, vector<int> &cluster_start, vector<int> &cluster_size, vector<int> &cluster_assign, vector<double> &cluster_radius, vector<vector<double> > &cluster_centroid) {
	int x, y, z;
	double radius = 0.0;
	double distance1 = 0.0;
	double distance2 = 0.0;

	//current datapoint is the start for the old_cluster and there is only one datapoint in that cluster
	if (cluster_start[old_cluster] == datapoint and cluster_size[old_cluster] == 1) {
		cluster_start[old_cluster] = -1;
		cluster_radius[old_cluster] = 0;
	}

	//current datapoint is the start for the old_cluster and there is more than one datapoint in that cluster
	else if (cluster_start[old_cluster] == datapoint and cluster_size[old_cluster] > 1) {
		for (x = datapoint + 1; x < ndata; x++) {
			if (cluster_assign[x] == old_cluster) {
				cluster_start[old_cluster] = x;
				for (y = 0; y < dim; y++) {
					distance2 += pow((cluster_centroid[new_cluster][y] - data[y]), 2);
				}
				distance2 = sqrt(distance2);

				if (distance2 > radius) {
					radius = distance2;
				}
			}
		}
		if (cluster_radius[old_cluster] < radius) {
			cluster_radius[old_cluster] = radius;
		}
	}

	//new_cluser's size is zero
	if (cluster_size[new_cluster] == 0) {
		cluster_start[new_cluster] = datapoint;
	}
	
	//new_cluster's start comes after the current datapoint
	 else if (cluster_start[new_cluster] > datapoint) {
		cluster_start[new_cluster] = datapoint;
	}

	for (x = 0; x < dim; x++) {
		distance1 += pow((cluster_centroid[new_cluster][x] - data[x]), 2);
	}
	distance1 = sqrt(distance1);

	if (cluster_radius[new_cluster] < distance1) {
		cluster_radius[new_cluster] = distance1;
	}

	cluster_size[old_cluster] -= 1;
	cluster_size[new_cluster] += 1;
	cluster_assign[datapoint] = new_cluster;
}

bool nearestCluster(int dim, int ndata, int k, int datapoint, double *data, vector<int> &cluster_start, vector<int> &cluster_size, vector<int> &cluster_assign, vector<double> &cluster_radius, vector<vector<double> > &cluster_centroid) {
	int x, y;
	int nearestCluster = -1;
	int old_cluster = cluster_assign[datapoint];
	double distance = 0.0;
	double minDistance = numeric_limits<double>::infinity();
	bool assignmentChange = false;

	for (x = 0; x < k; x++) {
		for (y = 0; y < dim; y++) {
			distance += pow((cluster_centroid[x][y] - data[y]), 2);
		}
		distance = sqrt(distance);

		if (minDistance > distance) {
			minDistance = distance;
			nearestCluster = x;
		}
		distance = 0;
	}

	if (cluster_assign[datapoint] != nearestCluster) {
		reassignment(datapoint, ndata, dim, old_cluster, nearestCluster, data, cluster_start, cluster_size, cluster_assign, cluster_radius, cluster_centroid);
		assignmentChange = true;
	}

	return assignmentChange;
}

void recalculateCentroid(int dim, int ndata, int cluster, double **data, int cluster_size, vector<int> &cluster_assign, vector<double> &cluster_centroid) {
	int x, y;
	double centroid = 0.0;

	for (x = 0; x < dim; x++) {
		for (y = 0; y < ndata; y++) {
			if (cluster_assign[y] == cluster) {
				centroid += data[y][x];
			}
		}
		cluster_centroid[x] = centroid / cluster_size;
	}
}

void intialNearestCluster(int dim, int k, int datapoint, double *data, vector<int> &cluster_start, vector<int> &cluster_size, vector<int> &cluster_assign, vector<double> &cluster_radius, vector<vector<double> > &cluster_centroid) {
	int x, y;
	int nearestCluster = -1;
	double distance = 0.0;
	double radius = 0;
	double minDistance = numeric_limits<double>::infinity();

	for (x = 0; x < k; x++) {
		for (y = 0; y < dim; y++) {
			distance += pow((cluster_centroid[x][y] - data[y]), 2);
		}
		distance = sqrt(distance);

		if (minDistance > distance) {
			radius = distance;
			minDistance = distance;
			nearestCluster = x;
		}
		distance = 0;
	}

	if (cluster_assign[datapoint] == -1 && cluster_size[nearestCluster] == 0) {
		cluster_start[nearestCluster] = datapoint;

		if (cluster_radius[nearestCluster] < radius) {
			cluster_radius[nearestCluster] = radius;
		}
	}

	cluster_assign[datapoint] = nearestCluster;
	cluster_size[nearestCluster] += 1;

}

void intial_centers(int dim, int ndata, double **data, int k, vector<vector<double> > &cluster_centroid) {

	random_device rd1;
	mt19937 randomGenerator1(rd1());
	uniform_int_distribution<int> distribution1(3, ndata - 1);
	uniform_real_distribution<double> distributionData(10, 100);

	int i, j;
	int x = 0;
	int isDuplicate = 0;
	int *dataset = (int*)calloc(ndata, sizeof(int));

	// Pick centroids as random points from the dataset.
	/*for (i = 0; i < k; i++) {
		while (isDuplicate != 1) {
			x = distribution1(randomGenerator1);
			if (dataset[x] == 0) {
				dataset[x] = 1;
				isDuplicate = 1;
			}
		}
		isDuplicate = 0;

		for (j = 0; j < dim; j++) {
			cluster_centroid[i][j] = data[x][j];
			cout << cluster_centroid[i][j] << endl;
		}
	}*/

	// Pick centroids as completey random points
	for (i = 0; i < k; i++) {
		for (j = 0; j < dim; j++) {
			cluster_centroid[i][j] = distributionData(randomGenerator1);
		}
	}

	delete[] dataset;
}

void print(int dim, int ndata, int k, vector<int> &cluster_start, vector<int> &cluster_size, vector<int> &cluster_assign, double **data, vector<double> &cluster_radius, vector<vector<double> > &cluster_centroid) {

	int i, j;

	cout << "KMeans has finished with these following values" << endl;
	cout << "\tdim is " << dim << endl;
	cout << "\tndata is " << ndata << endl;
	cout << "\tk is " << k << endl;

	cout << "\tdata has the following values" << endl;
	for (i = 0; i < ndata; i++) {
		cout << "\n\t\t";
		for (j = 0; j < dim; j++) {
			cout << data[i][j] << "  ";
		}
	}
	cout << "\n\n";

	cout << "\tcluster_size has the following values\n\n\t\t";
	for (i = 0; i < k; i++) {
		cout << cluster_size[i] << "  ";
	}
	cout << "\n\n";

	cout << "\tcluster_start has the following values\n\n\t\t";
	for (i = 0; i < k; i++) {
		cout << cluster_start[i] << "  ";
	}
	cout << "\n\n";

	cout << "\tcluster_assign has the following values\n\n\t\t";
	for (i = 0; i < ndata; i++) {
		cout << cluster_assign[i] << "  ";
	}
	cout << "\n\n";

	cout << "\tcluster_radius has the following values" << endl;
	for (i = 0; i < k; i++) {
		cout << "\n\t\t";
		cout << cluster_radius[i] << "  ";
	}
	cout << "\n\n";

	cout << "\tcluster_centroid has the following values" << endl;
	for (i = 0; i < k; i++) {
		cout << "\n\t\t";
		for (j = 0; j < dim; j++) {
			cout << cluster_centroid[i][j] << "  ";
		}
	}

}

int kmeans(int dim, int ndata, double **data, int k, vector<int> &cluster_start, vector<int> &cluster_size, vector<double> &cluster_radius, vector<vector<double> > &cluster_centroid){
	int i;
	bool exit = false;
	int count = 0;
	vector<int> cluster_assign(ndata, -1);

	intial_centers(dim, ndata, data, k, cluster_centroid);

	for (i = 0; i < ndata; i++) {
		intialNearestCluster(dim, k, i, data[i], cluster_start, cluster_size, cluster_assign, cluster_radius, cluster_centroid);
	}

	while (!exit){
		for (i = 0; i < k; i++) {
			if (cluster_size[i] != 0) {
				recalculateCentroid(dim, ndata, i, data, cluster_size[i], cluster_assign, cluster_centroid[i]);
			}
		}
		for (i = 0; i < ndata; i++) {
			exit = nearestCluster(dim, ndata, k, i, data[i], cluster_start, cluster_size, cluster_assign, cluster_radius, cluster_centroid);
		}
	}

	for (i = 0; i < k; i++) {
		if (cluster_size[i] != 0) {
			count += 1;
		}
	}

	print(dim, ndata, k, cluster_start, cluster_size, cluster_assign, data, cluster_radius, cluster_centroid);
	
	vector<int>().swap(cluster_assign);

	return count;
}

int main()
{
	random_device rd;
	mt19937 randomGenerator(rd());
	uniform_int_distribution<int> distributionDim(3, 8);
	uniform_int_distribution<int> distributionNdata(100, 10000);
	uniform_real_distribution<double> distributionData(10, 100);
	uniform_real_distribution<double> distributionW(0, 1);

	int i, j;
	int dim = distributionDim(randomGenerator);
	int ndata = distributionNdata(randomGenerator);
	int k = 100;
	double **data = (double**)calloc(ndata, sizeof(double));

	vector<int> cluster_start(k, 0);
	vector<int> cluster_size(k, 0);
	vector<double> cluster_radius(k, 0);
	vector<double> centroid(dim, 0);
	vector<vector<double> > cluster_centroid(k, centroid);

	//Allocate data[][]
	for (i = 0; i < ndata; i++) {
		data[i] = (double*)calloc(dim, sizeof(double));
	}

	//Assign random values for data
	for (i = 0; i < ndata; i++) {
		for (j = 0; j < dim; j++) {
			data[i][j] = distributionData(randomGenerator);
		}
	}
	
	cout << "\n" << kmeans(dim, ndata, data, k, cluster_start, cluster_size, cluster_radius, cluster_centroid) << " non-empty clusters found" << endl;

	vector<int>().swap(cluster_start);
	vector<int>().swap(cluster_size);
	vector<double>().swap(cluster_radius);
	vector<vector<double> >().swap(cluster_centroid);
	vector<double>().swap(centroid);
	delete[] data;

	getchar();
}
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

//h is a array of 'unit vectors'
//nclusters_ptr is used for ??
//cluster_hashvalues are going to be growing dependent on how long LSH
//Might be able to use vector within vector assignment so you can grow the size linerally or realloc*/
//hashvalues is found using the dot product

void resizeHashFunctions(vector<vector<double> > &h, const int rows, const int columns, int old_size) {
	
	random_device rd2;

	mt19937 randomGenerator2(rd2());

	uniform_real_distribution<double> distributionData2(10, 100);
	
	int i, j;

	h.resize(rows);
	for (auto &it : h){
		it.resize(columns);
	}

	for (i = old_size; i < rows; i++) {
		for (j = 0; j < columns; j++) {
			h[i][j] = distributionData2(randomGenerator2);
		}
	}
}

void resizeClusterVals(vector<vector<int> > &cluster_hashval, const int rows, const int columns) {
	cluster_hashval.resize(rows);
	for (auto &it : cluster_hashval) {
		it.resize(columns, 0);
	}
}

void print(int dim, int ndata, int m, vector<int> &ncluster_ptr, vector<int> &cluster_start, vector<int> &cluster_size, vector<vector<int> > &cluster_hashval, double W, double **data, vector<vector<double> > &h, size_t old_size) {

	int i, j;

	cout << "LSH has finished with these following values" << endl;
	cout << "\tdim is " << dim << endl;
	cout << "\tndata is " << ndata << endl;
	cout << "\tm is " << m << endl;
	cout << "\tW is " << W << endl;

	cout << "\tdata has the following values" << endl;
	for (i = 0; i < ndata; i++) {
		cout << "\n\t\t";
		for (j = 0; j < dim; j++) {
			cout << data[i][j] << "  ";
		}
	}
	cout << "\n\n";

	cout << "\tncluster_ptr has the following values\n\n\t\t";
	for (i = 0; i < ndata; i++) {
		cout << ncluster_ptr[i] << "  ";
	}
	cout << "\n\n";

	cout << "\tcluster_size has the following values\n\n\t\t";
	for (i = 0; i < m; i++) {
		cout << cluster_size[i] << "  ";
	}
	cout << "\n\n";

	cout << "\tcluster_start has the following values\n\n\t\t";
	for (i = 0; i < m; i++) {
		cout << cluster_start[i] << "  ";
	}
	cout << "\n\n";

	cout << "\tcluster_hashval has the following values" << endl;
	for (i = 0; i < m; i++) {
		cout << "\n\t\t";
		for (j = 0; j < old_size; j++) {
			cout << cluster_hashval[i][j] << "  ";
		}
	}
	cout << "\n\n";

	cout << "\th has the following values" << endl;
	for (i = 0; i < m; i++) {
		cout << "\n\t\t";
		for (j = 0; j < dim; j++) {
			cout << h[i][j] << "  ";
		}
	}

}

void findUnitVector(vector<double> &h, int dim, int x) {
	
	int i;
	double magnitude = 0;

	for (i = 0; i < dim; i++) {
		magnitude += pow(h[i], 2);
	}

	magnitude = sqrt(magnitude);
	for (i = 0; i < dim; i++) {
		h[i] = h[i] / magnitude;
	}
}

void lsh(int dim, int ndata, int m, vector<int> &ncluster_ptr, vector<int> &cluster_start, vector<int> &cluster_size, vector<vector<int> > &cluster_hashval, double W, double **data, vector<vector<double> > &h){

	int x, y, z;
	int i,j;
	int temp_ptr = 0;
	double dtemp_ptr = 0;
	int dotProduct = 0;
	int probe = -1;
	int beginProbe = 0;
	int old_m = m;
	bool realloc = false;

	cout << "Starting LSH" << endl;
	for (x = 0; x < m; x++) {
		findUnitVector(h[x], dim, x);
		for (y = 0; y < ndata; y++) {
			if (ncluster_ptr[y] != -1) {
				continue;
			}
			for (z = 0; z < dim; z++) {
				dotProduct += (data[y][z] * h[x][z]);
			}
			dotProduct = floor(dotProduct/W);
			do {
				if (probe == -1) {
					probe = dotProduct % old_m;
					beginProbe = probe;
				}

				//Set new cluster
				if (cluster_hashval[x][probe] == 0 and cluster_size[x] == 0) {
					cluster_hashval[x][probe] = dotProduct;
					cluster_start[x] = y;
					cluster_size[x] += 1;
					ncluster_ptr[y] = x;
					probe = -1;
					break;
				}

				//New hashval in already set cluster
				else if (cluster_hashval[x][probe] == 0 and cluster_size[x] > 0){
					cluster_hashval[x][probe] = dotProduct;
					cluster_size[x] += 1;
					ncluster_ptr[y] = x;
					probe = -1;
					break;;
				}

				//Same hash value in already set cluster
				else if (cluster_hashval[x][probe] == dotProduct) {
					cluster_size[x] += 1;
					ncluster_ptr[y] = x;
					probe = -1;
					break;
				}

				if (probe +1 == old_m) {
					probe = 0;
				}
				else {
					probe += 1;
				}

			} while (probe != beginProbe);

			if (probe == beginProbe) {
				realloc = true;
			}

			dotProduct = 0;
			probe = -1;
		}
		if (realloc == true and x == m - 1) {
			size_t old_size = m;
			size_t new_size = m * 2;
			cluster_size.resize(new_size);
			cluster_start.resize(new_size);
			cluster_hashval.resize(new_size);
			resizeHashFunctions(h, new_size, dim, old_size);
			resizeClusterVals(cluster_hashval, new_size, old_m);
			realloc = false;
			m = m * 2;
		}
		else if (realloc == false) {
			break;
		}
		else {
			realloc = false;
		}
		//print(dim, ndata, m, ncluster_ptr, cluster_start, cluster_size, cluster_hashval, W, data, h, old_m);
	}
	print(dim, ndata, m, ncluster_ptr, cluster_start, cluster_size, cluster_hashval, W, data, h, old_m);
	cout << "\nLSH Ended" << endl;
}

int main()
{
	random_device rd;

	mt19937 randomGenerator(rd());

	uniform_int_distribution<int> distributionDim(3, 8);
	uniform_int_distribution<int> distributionChoice(0, 1);
	uniform_int_distribution<int> distributionNdata(10, 1000);
	uniform_real_distribution<double> distributionData(10, 100);
	uniform_real_distribution<double> distributionW(0, 1);	

	int i, j;
	int dim = distributionDim(randomGenerator);
	int ndata = distributionNdata(randomGenerator);
	int m = dim - 2;
	double W = distributionW(randomGenerator);
	double **data = (double**)calloc(ndata, sizeof(double));
	
	vector<int> ncluster_ptr(ndata, -1);
	vector<int> cluster_start(m, 0);
	vector<int> cluster_size(m, 0);
	vector<int> cluster_hash(m, 0);
	vector<vector<int> > cluster_hashval(m, cluster_hash);
	vector<double> h1(dim, 0);
	vector<vector<double> > h(m, h1);

	//Allocate data[][]
	for (i = 0; i < ndata; i++) {
		data[i] = (double*)calloc(dim, sizeof(double));
	}
	
	//Assign random values for h
	for (i = 0; i < m; i++) {
		for (j = 0; j < dim; j++) {
			h[i][j] = distributionData(randomGenerator);
		}
	}

	//Assign random values for data
	for (i = 0; i < ndata; i++) {
		for (j = 0; j < dim; j++) {
			data[i][j] = distributionData(randomGenerator);
		}
	}

	lsh(dim, ndata, m, ncluster_ptr, cluster_start, cluster_size, cluster_hashval, W, data, h);

	vector<int>().swap(ncluster_ptr);
	vector<int>().swap(cluster_start);
	vector<int>().swap(cluster_size);
	vector<vector<int> >().swap(cluster_hashval);
	vector<int>().swap(cluster_hash);
	vector<vector<double> >().swap(h);
	vector<double>().swap(h1);
	delete[] data;


	getchar();

}
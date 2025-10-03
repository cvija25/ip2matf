#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>
#include <iomanip>

using namespace std;

using namespace std;

using Matrix = vector<vector<double>>;

// Load CSV into matrix, skipping the header (first row)
// and ignoring the first 3 columns of each row
Matrix loadCSV(const string &filename) {
    ifstream file(filename);
    Matrix data;
    string line;

    bool firstRow = true;
    while (getline(file, line)) {
        if (firstRow) { // skip header row
            firstRow = false;
            continue;
        }
        stringstream ss(line);
        string cell;
        vector<double> row;
        int colIndex = 0;
        while (getline(ss, cell, ',')) {
            if (colIndex >= 3) { // skip first 3 columns
                row.push_back(stod(cell));
            }
            colIndex++;
        }
        if (!row.empty())
            data.push_back(row);
    }
    return data;
}

// Transpose matrix
Matrix transpose(const Matrix &m) {
    Matrix t(m[0].size(), vector<double>(m.size()));
    for (size_t i = 0; i < m.size(); i++)
        for (size_t j = 0; j < m[0].size(); j++)
            t[j][i] = m[i][j];
    return t;
}

// Compute variance of vector
double variance(const vector<double> &v) {
    double mean = 0.0;
    for (double x : v) mean += x;
    mean /= v.size();
    double var = 0.0;
    for (double x : v) var += (x - mean) * (x - mean);
    return var / v.size();
}

// Standardize features (z-score)
void normalize(Matrix &data) {
    size_t rows = data.size(), cols = data[0].size();
    for (size_t j = 0; j < cols; j++) {
        double mean = 0, sd = 0;
        for (size_t i = 0; i < rows; i++) mean += data[i][j];
        mean /= rows;
        for (size_t i = 0; i < rows; i++) sd += pow(data[i][j] - mean, 2);
        sd = sqrt(sd / rows);
        if (sd == 0) sd = 1;
        for (size_t i = 0; i < rows; i++)
            data[i][j] = (data[i][j] - mean) / sd;
    }
}

// K-means clustering
vector<int> kmeans(const Matrix &data, int k, int maxIter = 100) {
    size_t rows = data.size(), cols = data[0].size();
    vector<vector<double>> centroids(k, vector<double>(cols));
    vector<int> labels(rows);

    // Init centroids randomly
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, rows - 1);
    for (int i = 0; i < k; i++) centroids[i] = data[dis(gen)];

    for (int iter = 0; iter < maxIter; iter++) {
        // Assign clusters
        for (size_t i = 0; i < rows; i++) {
            double bestDist = numeric_limits<double>::max();
            int bestCluster = 0;
            for (int c = 0; c < k; c++) {
                double dist = 0;
                for (size_t j = 0; j < cols; j++)
                    dist += pow(data[i][j] - centroids[c][j], 2);
                if (dist < bestDist) {
                    bestDist = dist;
                    bestCluster = c;
                }
            }
            labels[i] = bestCluster;
        }

        // Update centroids
        vector<vector<double>> newCentroids(k, vector<double>(cols, 0.0));
        vector<int> counts(k, 0);
        for (size_t i = 0; i < rows; i++) {
            int c = labels[i];
            counts[c]++;
            for (size_t j = 0; j < cols; j++)
                newCentroids[c][j] += data[i][j];
        }
        for (int c = 0; c < k; c++) {
            if (counts[c] == 0) continue;
            for (size_t j = 0; j < cols; j++)
                newCentroids[c][j] /= counts[c];
        }
        centroids = newCentroids;
    }
    return labels;
}

Matrix select_top_variance(const Matrix &data, int topN) {
    int rows = data.size();
    if (rows == 0) return {};
    int cols = data[0].size();

    // 1. Compute variance for each column
    vector<pair<double, int>> geneVars; // (variance, columnIndex)
    for (int j = 0; j < cols; j++) {
        vector<double> column;
        column.reserve(rows);
        for (int i = 0; i < rows; i++) column.push_back(data[i][j]);
        double var = variance(column);
        geneVars.push_back({var, j});
    }

    // 2. Sort by variance descending
    sort(geneVars.rbegin(), geneVars.rend());

    // 3. Pick topN indices
    int keepN = min(topN, cols);
    vector<int> keepIdx;
    keepIdx.reserve(keepN);
    for (int k = 0; k < keepN; k++) keepIdx.push_back(geneVars[k].second);

    // 4. Build reduced matrix
    Matrix reduced(rows, vector<double>(keepN));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < keepN; j++) {
            reduced[i][j] = data[i][keepIdx[j]];
        }
    }

    return reduced;
}

// Euclidean norm for later if needed
double norm(const vector<double>& v) {
    double sum = 0;
    for (double x : v) sum += x*x;
    return sqrt(sum);
}

// PCA function
vector<vector<double>> PCA(const vector<vector<double>>& data, int k) {
    int n = data.size();
    int m = data[0].size();

    // 1. Compute column means
    vector<double> mean(m, 0.0);
    for(int j=0;j<m;j++)
        for(int i=0;i<n;i++)
            mean[j] += data[i][j];
    for(int j=0;j<m;j++)
        mean[j] /= n;

    // 2. Center data
    vector<vector<double>> X = data;
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            X[i][j] -= mean[j];

    // 3. Compute covariance matrix
    vector<vector<double>> cov(m, vector<double>(m,0.0));
    for(int i=0;i<m;i++)
        for(int j=0;j<=i;j++){
            double sum = 0;
            for(int r=0;r<n;r++) sum += X[r][i]*X[r][j];
            cov[i][j] = cov[j][i] = sum/(n-1);
        }

    // 4. Jacobi eigen decomposition (for small/medium matrices)
    vector<double> eigenvalues(m,0.0);
    vector<vector<double>> eigenvectors(m, vector<double>(m,0.0));
    for(int i=0;i<m;i++) eigenvectors[i][i]=1.0;

    int maxIter = 100;
    for(int iter=0;iter<maxIter;iter++){
        int p=0,q=1; double maxVal=fabs(cov[p][q]);
        for(int i=0;i<m;i++)
            for(int j=i+1;j<m;j++)
                if(fabs(cov[i][j])>maxVal) { maxVal=fabs(cov[i][j]); p=i;q=j; }
        if(maxVal<1e-10) break;

        double theta = 0.5*atan2(2*cov[p][q], cov[q][q]-cov[p][p]);
        double c=cos(theta), s=sin(theta);
        double app = c*c*cov[p][p] - 2*s*c*cov[p][q] + s*s*cov[q][q];
        double aqq = s*s*cov[p][p] + 2*s*c*cov[p][q] + c*c*cov[q][q];
        cov[p][p]=app; cov[q][q]=aqq; cov[p][q]=cov[q][p]=0.0;

        for(int r=0;r<m;r++){
            if(r!=p && r!=q){
                double arp = c*cov[r][p]-s*cov[r][q];
                double arq = s*cov[r][p]+c*cov[r][q];
                cov[r][p]=cov[p][r]=arp;
                cov[r][q]=cov[q][r]=arq;
            }
        }

        for(int r=0;r<m;r++){
            double vrp = c*eigenvectors[r][p]-s*eigenvectors[r][q];
            double vrq = s*eigenvectors[r][p]+c*eigenvectors[r][q];
            eigenvectors[r][p]=vrp;
            eigenvectors[r][q]=vrq;
        }
    }

    for(int i=0;i<m;i++) eigenvalues[i] = cov[i][i];

    // 5. Sort eigenvalues descending
    vector<int> indices(m);
    for(int i=0;i<m;i++) indices[i]=i;
    sort(indices.begin(), indices.end(), [&](int a,int b){return eigenvalues[a]>eigenvalues[b];});

    // 6. Project data onto top-k eigenvectors
    vector<vector<double>> reduced(n, vector<double>(k,0.0));
    for(int i=0;i<n;i++)
        for(int j=0;j<k;j++)
            for(int f=0;f<m;f++)
                reduced[i][j] += X[i][f]*eigenvectors[f][indices[j]];

    return reduced;
}

void saveCSV(const Matrix &data, const string &filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: could not open " << filename << " for writing\n";
        return;
    }

    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data[i].size(); j++) {
            file << fixed << setprecision(6) << data[i][j];
            if (j + 1 < data[i].size()) file << ",";
        }
        file << "\n";
    }
    file.close();
    cout << "Saved matrix to " << filename << " (" 
         << data.size() << " × " << data[0].size() << ")\n";
}

int main() {
    // Load dataset
    Matrix data = loadCSV("all_BS2_1.csv");

    // Transpose to filter low-variance genes
    Matrix dataT = transpose(data);
    cout << "duzina prve " << data.size() << "x" << data[0].size() << endl;
    Matrix filtered;
    for (auto &gene : dataT) {
        if (variance(gene) > 0.01) {
            filtered.push_back(gene);
        }
    }
    
    // Transpose back (samples x genes)
    Matrix cleanData = transpose(filtered);
    cout << "duzina druge " << cleanData.size() << "x" << cleanData[0].size() << endl;


    Matrix newData = select_top_variance(cleanData, 2000);
    cout << "duzina druge " << newData.size() << "x" << newData[0].size() << endl;

    Matrix X_pca = PCA(newData, 50);

    cout << "Reduced (PCA): " << X_pca.size() << " × " << X_pca[0].size() << endl;

    saveCSV(X_pca, "reduced_genes.csv");

    // // Normalize
    // normalize(cleanData);

    // // Cluster into 3 groups
    // vector<int> labels = kmeans(cleanData, 3);

    // // Output results
    // for (size_t i = 0; i < labels.size(); i++)
    //     cout << "Sample " << i << " -> Cluster " << labels[i] << "\n";

    return 0;
}
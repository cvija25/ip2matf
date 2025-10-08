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

using Matrix = vector<vector<double>>;

// ucitavamo podatke iz .csv
Matrix loadCSV(const string &filename) {
    ifstream file(filename);
    Matrix data;
    string line;

    // ucitavamo liniju po liniju
    bool firstRow = true;
    while (getline(file, line)) {
        if (firstRow) { 
            // preskacemo imena atributa
            firstRow = false;
            continue;
        }
        stringstream ss(line);
        string buffer;
        vector<double> row;
        int colIndex = 0;
        
        // citamo element po element sa zarezom kao delimiter
        while (getline(ss, buffer, ',')) {
            // preskacemo prva tri atributa (src,Cell_type,file,row)
            if (colIndex >= 4) {
                row.push_back(stod(buffer));
            }
            colIndex++;
        }
        if (!row.empty())
            data.push_back(row);
    }
    return data;
}

// funkcija koja vraca transponovanu matricu
Matrix transpose(const Matrix &m) {
    Matrix t(m[0].size(), vector<double>(m.size()));
    for (int i = 0; i < m.size(); i++)
        for (int j = 0; j < m[0].size(); j++)
            t[j][i] = m[i][j];
    return t;
}

// funkcija koja racuna varijansu vektora
double variance(const vector<double> &v) {
    double mean = 0.0;
    for (double x : v) mean += x;
    mean /= v.size();
    double var = 0.0;
    for (double x : v) var += (x - mean) * (x - mean);
    return var / v.size();
}

void normalize(Matrix &data) {
    int n = data.size(), m = data[0].size();
    for (int j = 0; j < m; j++) {
        double mean = 0, sd = 0;
        for (int i = 0; i < n; i++) mean += data[i][j];
        mean /= n;
        for (size_t i = 0; i < n; i++) sd += (data[i][j] - mean) * (data[i][j] - mean);
        sd = sqrt(sd / n);
        if (sd == 0) sd = 1;
        for (int i = 0; i < n; i++)
            data[i][j] = (data[i][j] - mean) / sd;
    }
}

Matrix select_top_variance(const Matrix &data, int topN) {
    int n = data.size();
    int m = data[0].size();

    // varijansa svakog atributa
    vector<pair<double, int>> geneVars;
    for (int j = 0; j < m; j++) {
        vector<double> column(n);
        for (int i = 0; i < n; i++) column[i] = data[i][j];
        double var = variance(column);
        geneVars.push_back({var, j});
    }

    sort(geneVars.rbegin(), geneVars.rend());

    vector<int> keepIdx;
    // uzimamo prvih topN
    for (int k = 0; k < topN; k++) keepIdx.push_back(geneVars[k].second);

    Matrix reduced(n, vector<double>(topN));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < topN; j++) {
            reduced[i][j] = data[i][keepIdx[j]];
        }
    }

    return reduced;
}


// PCA function
vector<vector<double>> PCA(const vector<vector<double>>& data, int k) {
    int n = data.size();
    int m = data[0].size();

    // racunamo proseke atributa
    vector<double> mean(m, 0.0);
    for(int j=0;j<m;j++)
        for(int i=0;i<n;i++)
            mean[j] += data[i][j];
    for(int j=0;j<m;j++)
        mean[j] /= n;

    // centriramo ih
    vector<vector<double>> X = data;
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            X[i][j] -= mean[j];

    // matrica kovariansi
    vector<vector<double>> cov(m, vector<double>(m,0.0));
    for(int i=0;i<m;i++)
        for(int j=0;j<=i;j++){
            double sum = 0;
            for(int r=0;r<n;r++) sum += X[r][i]*X[r][j];
            cov[i][j] = cov[j][i] = sum/(n-1);
        }

    // jakovijeva metoda sa sopstvenim vrednostima
    vector<double> eigenvalues(m,0.0);
    vector<vector<double>> eigenvectors(m, vector<double>(m,0.0));
    for(int i=0;i<m;i++) eigenvectors[i][i]=1.0;

    int maxIter = 100;
    for(int iter=0;iter<maxIter;iter++){
        int p=0,q=1; double maxVal=fabs(cov[p][q]);
        for(int i=0;i<m;i++)
            for(int j=i+1;j<m;j++)
                if(fabs(cov[i][j])>maxVal) { maxVal=fabs(cov[i][j]); p=i;q=j; }
        // nadjemo najveci vandijagonalni element
        if(maxVal<1e-10) break;

        // nadjemo ugao rotacije
        double theta = 0.5*atan2(2*cov[p][q], cov[q][q]-cov[p][p]);
        double c=cos(theta), s=sin(theta);

        // primenimo rotaciju
        double app = c*c*cov[p][p] - 2*s*c*cov[p][q] + s*s*cov[q][q];
        double aqq = s*s*cov[p][p] + 2*s*c*cov[p][q] + c*c*cov[q][q];
        cov[p][p]=app; cov[q][q]=aqq; cov[p][q]=cov[q][p]=0.0;


        // azuriramo povezane elemente
        for(int r=0;r<m;r++){
            if(r!=p && r!=q){
                double arp = c*cov[r][p]-s*cov[r][q];
                double arq = s*cov[r][p]+c*cov[r][q];
                cov[r][p]=cov[p][r]=arp;
                cov[r][q]=cov[q][r]=arq;
            }
        }

        // svrstavamo u matricu
        for(int r=0;r<m;r++){
            double vrp = c*eigenvectors[r][p]-s*eigenvectors[r][q];
            double vrq = s*eigenvectors[r][p]+c*eigenvectors[r][q];
            eigenvectors[r][p]=vrp;
            eigenvectors[r][q]=vrq;
        }
    }

    for(int i=0;i<m;i++) eigenvalues[i] = cov[i][i];

    // sortiramo sopstvene vektore
    vector<int> indices(m);
    for(int i=0;i<m;i++) indices[i]=i;
    sort(indices.begin(), indices.end(), [&](int a,int b){return eigenvalues[a]>eigenvalues[b];});

    // pronadjemo k najboljih
    vector<vector<double>> reduced(n, vector<double>(k,0.0));
    for(int i=0;i<n;i++)
        for(int j=0;j<k;j++)
            for(int f=0;f<m;f++)
                reduced[i][j] += X[i][f]*eigenvectors[f][indices[j]];

    return reduced;
}

void saveCSV(const Matrix &data, const string &filename) {
    ofstream file(filename);

    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].size(); j++) {
            file << fixed << setprecision(6) << data[i][j];
            if (j + 1 < data[i].size()) file << ",";
        }
        file << "\n";
    }
    file.close();
}

int main() {
    // ucitavamo skup podataka
    Matrix data = loadCSV("all_BS2_1.csv");

    // transponujemo kako bismo lakse filtrili atribute
    Matrix dataT = transpose(data);
    cout << "dimenzije pocetne " << data.size() << "x" << data[0].size() << endl;
    Matrix filtered;
    for (auto &gene : dataT) {
        if (variance(gene) > 0.01) {
            filtered.push_back(gene);
        }
    }
    
    Matrix cleanData = transpose(filtered);
    cout << "dimenzije posle filtriranja " << cleanData.size() << "x" << cleanData[0].size() << endl;


    Matrix newData = select_top_variance(cleanData, 2000);
    cout << "dimenzije posle biranja najvarijantnijih " << newData.size() << "x" << newData[0].size() << endl;

    normalize(newData);

    Matrix X_pca = PCA(newData, 50);

    cout << "dimenzije posle PCA " << X_pca.size() << " × " << X_pca[0].size() << endl;

    saveCSV(X_pca, "reduced_genes.csv");

    return 0;
}
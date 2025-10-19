#include <iostream>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <cstring>

using namespace std;

struct Point {
    double x, y;
    int cluster;
    
    Point(double x, double y, int cluster) 
        : x(x), y(y), cluster(cluster) {}
};

vector<vector<double>> loadAttributes(const string& filename) {
    vector<vector<double>> data;
    ifstream file(filename.c_str());
    
    if (!file.is_open()) {
        cerr << "Error: cannot open file " << filename << endl;
        return data;
    }
    
    string line;
    int lineNum = 0;
    while (getline(file, line)) {
        lineNum++;
        if (lineNum == 1) continue;
        if (line.empty()) continue;
        
        vector<double> row;
        char* line_copy = new char[line.length() + 1];
        strcpy(line_copy, line.c_str());
        
        char* token = strtok(line_copy, ",");
        while (token != NULL) {
            row.push_back(atof(token));
            token = strtok(NULL, ",");
        }
        delete[] line_copy;
        
        if (!row.empty()) {
            data.push_back(row);
        }
    }
    
    file.close();
    return data;
}

vector<pair<int, int>> loadClusters(const string& filename) {
    vector<pair<int, int>> data;
    ifstream file(filename.c_str());
    
    if (!file.is_open()) {
        cerr << "Error: cannot open file " << filename << endl;
        return data;
    }
    
    string line;
    int lineNum = 0;
    while (getline(file, line)) {
        lineNum++;
        if (lineNum == 1) continue;
        if (line.empty()) continue;
        
        int idx, cluster;
        int parsed = sscanf(line.c_str(), "%d,%d", &idx, &cluster);
        
        if (parsed == 2) {
            data.push_back(make_pair(idx, cluster));
        }
    }
    
    file.close();
    return data;
}

vector<Point> combineData(const vector<vector<double>>& attributes, 
                          const vector<pair<int, int>>& clusters) {
    vector<Point> points;
    
    for (int i = 0; i < clusters.size(); i++) {
        int idx = clusters[i].first;
        int cluster = clusters[i].second;
        
        if (idx >= 0 && idx < attributes.size() && attributes[idx].size() >= 2) {
            double x = attributes[idx][0];
            double y = attributes[idx][1];
            points.push_back(Point(x, y, cluster));
        }
    }
    
    return points;
}

void findBounds(const vector<Point>& data, 
                double& minX, double& maxX, 
                double& minY, double& maxY) {
    minX = data[0].x;
    maxX = data[0].x;
    minY = data[0].y;
    maxY = data[0].y;
    
    for (int i = 0; i < data.size(); i++) {
        if (data[i].x < minX) minX = data[i].x;
        if (data[i].x > maxX) maxX = data[i].x;
        if (data[i].y < minY) minY = data[i].y;
        if (data[i].y > maxY) maxY = data[i].y;
    }
}

void scaleCoords(double x, double y, 
                 double minX, double maxX, 
                 double minY, double maxY,
                 int width, int height, int margin,
                 double& sx, double& sy) {
    sx = margin + (x - minX) / (maxX - minX) * (width - 2 * margin);
    sy = height - margin - (y - minY) / (maxY - minY) * (height - 2 * margin);
}

string getColor(int cluster) {
    switch(cluster) {
        case 0: return "#FF6B6B";
        case 1: return "#4ECDC4";
        case 2: return "#45B7D1";
        case 3: return "#FFA07A";
        case 4: return "#98D8C8";
        case 5: return "#F7DC6F";
        case 6: return "#95E1D3";
        case 7: return "#F38181";
        default: return "#333333";
    }
}

void writeSVGHeader(ofstream& file, int width, int height, const vector<pair<int, int>>& clusters) {
    file << "<svg width='" << width << "' height='" << height 
         << "' xmlns='http://www.w3.org/2000/svg'>\n";
    file << "<rect width='100%' height='100%' fill='white'/>\n";
    
    int clusterCounts[10] = {0};
    int maxCluster = 0;
    for (int i = 0; i < clusters.size(); i++) {
        int c = clusters[i].second;
        clusterCounts[c]++;
        if (c > maxCluster) maxCluster = c;
    }
    
    file << "<text x='10' y='20' font-size='16' font-weight='bold'>Clusters:</text>\n";
    for (int i = 0; i <= maxCluster; i++) {
        file << "<circle cx='15' cy='" << (40 + i * 20) << "' r='6' fill='" << getColor(i) << "'/>\n";
        file << "<text x='30' y='" << (45 + i * 20) << "' font-size='12'>Cluster " << i 
             << " (" << clusterCounts[i] << ")</text>\n";
    }
    
    file << "<text x='10' y='" << (height - 20) << "' font-size='12' fill='#666'>X →</text>\n";
    file << "<text x='" << (width - 30) << "' y='20' font-size='12' fill='#666'>↑ Y</text>\n";
}

void writeSVGFooter(ofstream& file) {
    file << "</svg>";
}

void plotPoints(ofstream& file, const vector<Point>& data,
                double minX, double maxX, double minY, double maxY,
                int width, int height, int margin) {
    for (int i = 0; i < data.size(); i++) {
        double sx, sy;
        scaleCoords(data[i].x, data[i].y, minX, maxX, minY, maxY, 
                    width, height, margin, sx, sy);
        
        file << "<circle cx='" << sx << "' cy='" << sy 
             << "' r='6' fill='" << getColor(data[i].cluster) << "'/>\n";
    }
}

void saveSVG(const vector<Point>& data, const vector<pair<int, int>>& clusters, const string& filename) {
    ofstream file(filename.c_str());
    
    int width = 800;
    int height = 600;
    int margin = 50;
    
    double minX, maxX, minY, maxY;
    findBounds(data, minX, maxX, minY, maxY);
    
    writeSVGHeader(file, width, height, clusters);
    plotPoints(file, data, minX, maxX, minY, maxY, width, height, margin);
    writeSVGFooter(file);
    
    file.close();
}

int main() {
    vector<vector<double>> attributes = loadAttributes("reduced_genes.csv");
    vector<pair<int, int>> clusters = loadClusters("clustering_results.csv");
    vector<Point> data = combineData(attributes, clusters);

    if (data.empty()) {
        cerr << "No data combined." << endl;
        return 1;
    }
    
    saveSVG(data, clusters, "clusters_2d.svg");
    
    return 0;
}
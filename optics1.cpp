#include <iostream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include <limits>
#include <algorithm>
#include <fstream>


using namespace std;

const double UNDEFINED = numeric_limits<double>::max();

struct Point {
    vector<double> coords;
    bool processed = false;
    double reachabilityDist = UNDEFINED;

    Point(vector<double> c) : coords(c) {}
};

double euclideanDistance(const Point& p1, const Point& p2) {
    double dist = 0.0;
    for (size_t i = 0; i < p1.coords.size(); ++i) {
        dist += pow(p1.coords[i] - p2.coords[i], 2);
    }
    return sqrt(dist);
}

vector<Point> loadCsv(string filename) {
    ifstream file(filename);
    vector<Point> data;
    string line, cell;
    // first line
    getline(file, line);

    while (getline(file, line)) {
        stringstream lineStream(line);
        vector<double> row;
        while (getline(lineStream, cell, ',')) {
            row.push_back(stod(cell));
        }
        data.push_back(Point(row));
    }
    return data;
}
struct Optics {
    vector<Point>& points;
    vector<int>& processingOrder;
    int epsilon, minPts;

    Optics(int e, int mP, vector<Point>& p, vector<int>& po) : epsilon(e), minPts(mP), points(p), processingOrder(po) {}

    vector<int> getNeighbors(const vector<Point>& points, int idx) {
        vector<int> neighbors;
        for (size_t i = 0; i < points.size(); ++i) {
            if (i != idx && euclideanDistance(points[idx], points[i]) <= epsilon) {
                neighbors.push_back(i);
            }
        }
        return neighbors;
    }

    void updateReachability(priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> &seeds, const vector<int>& neighbors, int idx) {
        double coreDist = neighbors.size() >= minPts ? euclideanDistance(points[idx], points[neighbors[minPts - 1]]) : UNDEFINED;
        if (coreDist != UNDEFINED) {
            for (int neighbor : neighbors) {
                if (!points[neighbor].processed) {
                    double newReachDist = max(coreDist, euclideanDistance(points[idx], points[neighbor]));
                    if (points[neighbor].reachabilityDist == UNDEFINED || newReachDist < points[neighbor].reachabilityDist) {
                        points[neighbor].reachabilityDist = newReachDist;
                        seeds.push({newReachDist, neighbor});
                    }
                }
            }
        }
    }

    void fit() {
        for (size_t i = 0; i < points.size(); ++i) {
            if (!points[i].processed) {
                points[i].processed = true;
                processingOrder.push_back(i);
                vector<int> neighbors = getNeighbors(points, i);
                if (neighbors.size() >= minPts) {
                    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> seeds;
                    updateReachability(seeds, neighbors, i);
                    while (!seeds.empty()) {
                        int current = seeds.top().second;
                        seeds.pop();
                        if (!points[current].processed) {
                            points[current].processed = true;
                            processingOrder.push_back(current);
                            vector<int> currentNeighbors = getNeighbors(points, current);
                            if (currentNeighbors.size() >= minPts) {
                                updateReachability(seeds, currentNeighbors, current);
                            }
                        }
                    }
                }
            }
        }
    }

    vector<vector<int>> extractClusters(double clusterDistanceThreshold) {
        vector<vector<int>> clusters;
        vector<int> currentCluster;

        for (int idx : processingOrder) {
            if (points[idx].reachabilityDist > clusterDistanceThreshold) {
                if (!currentCluster.empty()) {
                    clusters.push_back(currentCluster);
                    currentCluster.clear();
                }
            } else {
                currentCluster.push_back(idx);
            }
        }

        // Add the last cluster if it exists
        if (!currentCluster.empty()) {
            clusters.push_back(currentCluster);
        }

        return clusters;
    }

    void exportClusters(const vector<vector<int>>& clusters, const string& filename) {
        ofstream file(filename);
        file << "x,y,cluster\n";

        // Assign cluster numbers
        for (size_t clusterId = 0; clusterId < clusters.size(); ++clusterId) {
            for (int idx : clusters[clusterId]) {
                file << points[idx].coords[0] << "," << points[idx].coords[1] << "," << clusterId + 1 << "\n";
            }
        }

        // Add noise points (not in any cluster)
        for (size_t i = 0; i < points.size(); ++i) {
            bool inCluster = false;
            for (const auto& cluster : clusters) {
                if (find(cluster.begin(), cluster.end(), i) != cluster.end()) {
                    inCluster = true;
                    break;
                }
            }
            if (!inCluster) {
                file << points[i].coords[0] << "," << points[i].coords[1] << ",0\n"; // Noise is labeled as cluster 0
            }
        }

        file.close();
    }

};


int main() {
    vector<Point> points = loadCsv("klasterovanje.csv");
    vector<int> processingOrder;
    Optics* model = new Optics(3.0, 7, points, processingOrder);


    double clusterDistanceThreshold = 2.5;

    model->fit();

    cout << "Point\tReachability Distance\n";
    for (size_t i = 0; i < points.size(); ++i) {
        cout << i << "\t";
        if (points[i].reachabilityDist == UNDEFINED)
            cout << "UNDEFINED\n";
        else
            cout << points[i].reachabilityDist << "\n";
    }

    // Output reachability plot
    cout << "\nReachability Plot (Processing Order):\n";
    for (int idx : model->processingOrder) {
        cout << idx << "\t";
        if (points[idx].reachabilityDist == UNDEFINED)
            cout << "UNDEFINED\n";
        else
            cout << points[idx].reachabilityDist << "\n";
    }

    ofstream outFile("reachability_plot.csv");
    outFile << "Index,ReachabilityDist\n";
    for (int idx : model->processingOrder) {
        outFile << idx << ",";
        if (points[idx].reachabilityDist == UNDEFINED)
            outFile << "UNDEFINED\n";
        else
            outFile << points[idx].reachabilityDist << "\n";
    }
    outFile.close();

    vector<vector<int>> clusters = model->extractClusters(clusterDistanceThreshold);

    cout << "\nClusters:\n";
    for (size_t i = 0; i < clusters.size(); ++i) {
        cout << "Cluster " << i + 1 << ": ";
        for (int idx : clusters[i]) {
            cout << idx << " ";
        }
        cout << "\n";
    }
    model->exportClusters(clusters, "clusters.csv");

    return 0;
}

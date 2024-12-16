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

vector<int> getNeighbors(const vector<Point>& points, int idx, double epsilon) {
    vector<int> neighbors;
    for (size_t i = 0; i < points.size(); ++i) {
        if (i != idx && euclideanDistance(points[idx], points[i]) <= epsilon) {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

void updateReachability(vector<Point>& points, priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> &seeds, const vector<int>& neighbors, int idx, double epsilon, int minPts) {
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

void optics(vector<Point>& points, double epsilon, int minPts, vector<int>& processingOrder) {
    for (size_t i = 0; i < points.size(); ++i) {
        if (!points[i].processed) {
            points[i].processed = true;
            processingOrder.push_back(i);
            vector<int> neighbors = getNeighbors(points, i, epsilon);
            if (neighbors.size() >= minPts) {
                priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> seeds;
                updateReachability(points, seeds, neighbors, i, epsilon, minPts);
                while (!seeds.empty()) {
                    int current = seeds.top().second;
                    seeds.pop();
                    if (!points[current].processed) {
                        points[current].processed = true;
                        processingOrder.push_back(current);
                        vector<int> currentNeighbors = getNeighbors(points, current, epsilon);
                        if (currentNeighbors.size() >= minPts) {
                            updateReachability(points, seeds, currentNeighbors, current, epsilon, minPts);
                        }
                    }
                }
            }
        }
    }
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

vector<vector<int>> extractClusters(const vector<Point>& points, const vector<int>& processingOrder, double clusterDistanceThreshold) {
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


#include <fstream>

void exportClusters(const vector<Point>& points, const vector<vector<int>>& clusters, const string& filename) {
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


int main() {
    vector<Point> points = loadCsv("klasterovanje.csv");

    double epsilon = 3.0;
    int minPts = 7;
    double clusterDistanceThreshold = 2.5;

    vector<int> processingOrder;
    optics(points, epsilon, minPts, processingOrder);

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
    for (int idx : processingOrder) {
        cout << idx << "\t";
        if (points[idx].reachabilityDist == UNDEFINED)
            cout << "UNDEFINED\n";
        else
            cout << points[idx].reachabilityDist << "\n";
    }

    ofstream outFile("reachability_plot.csv");
    outFile << "Index,ReachabilityDist\n";
    for (int idx : processingOrder) {
        outFile << idx << ",";
        if (points[idx].reachabilityDist == UNDEFINED)
            outFile << "UNDEFINED\n";
        else
            outFile << points[idx].reachabilityDist << "\n";
    }
    outFile.close();

    vector<vector<int>> clusters = extractClusters(points, processingOrder, clusterDistanceThreshold);

    cout << "\nClusters:\n";
    for (size_t i = 0; i < clusters.size(); ++i) {
        cout << "Cluster " << i + 1 << ": ";
        for (int idx : clusters[i]) {
            cout << idx << " ";
        }
        cout << "\n";
    }
    exportClusters(points, clusters, "clusters.csv");

    return 0;
}

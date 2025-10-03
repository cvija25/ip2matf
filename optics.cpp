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
    for (int i = 0; i < p1.coords.size(); ++i) {
        dist += pow(p1.coords[i] - p2.coords[i], 2);
    }
    return sqrt(dist);
}

vector<Point>* loadCsv(string filename) {
    ifstream file(filename);
    vector<Point>* data = new vector<Point>;
    string line, cell;
    // first line
    getline(file, line);

    while (getline(file, line)) {
        stringstream lineStream(line);
        vector<double> row;
        while (getline(lineStream, cell, ',')) {
            row.push_back(stod(cell));
        }
        data->push_back(Point(row));
    }
    return data;
}
struct Optics {
    vector<Point>* m_points;
    vector<int>* m_processingOrder;
    int epsilon, minPts;

    Optics(int e, int mP) : epsilon(e), minPts(mP), m_points(nullptr), m_processingOrder(nullptr) {}

    ~Optics() {
        delete m_points;
        delete m_processingOrder;
    }

    vector<int> getNeighbors(int idx) {
        vector<int> neighbors;
        vector<Point>& points = *m_points;
        for (size_t i = 0; i < points.size(); ++i) {
            if (i != idx && euclideanDistance(points[idx], points[i]) <= epsilon) {
                neighbors.push_back(i);
            }
        }
        return neighbors;
    }

    // goes through each unprocessed neighbor and updates it's reachability and puts them in priority queue
    void updateReachability(priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> &seeds, const vector<int>& neighbors, int idx) {
        vector<Point>& points = *m_points;
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

    void fit(vector<Point>& points) {
        m_points = &points;
        m_processingOrder = new vector<int>;

        // goes through each point
        for (int i = 0; i < points.size(); ++i) {
            if (!points[i].processed) {
                points[i].processed = true;
                m_processingOrder->push_back(i);
                vector<int> neighbors = getNeighbors(i);

                // if it's a core point
                if (neighbors.size() >= minPts) {
                    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> seeds;
                    updateReachability(seeds, neighbors, i);

                    // porcesses all "near" points
                    while (!seeds.empty()) {
                        int current = seeds.top().second;
                        seeds.pop();
                        if (!points[current].processed) {
                            points[current].processed = true;
                            m_processingOrder->push_back(current);
                            vector<int> currentNeighbors = getNeighbors(current);
                            if (currentNeighbors.size() >= minPts) {
                                updateReachability(seeds, currentNeighbors, current);
                            }
                        }
                    }
                }
            }
        }
    }

    // we define cluster by setting a threshold
    vector<vector<int>> extractClusters(double clusterDistanceThreshold) {
        vector<vector<int>> clusters;
        vector<int> currentCluster;
        vector<Point>& points = *m_points;

        for (int idx : *m_processingOrder) {
            if (points[idx].reachabilityDist < clusterDistanceThreshold) {
                currentCluster.push_back(idx);
            } else if (!currentCluster.empty()) {
                clusters.push_back(currentCluster);
                currentCluster.clear();
            }
        }

        // add the last cluster if it exists
        if (!currentCluster.empty()) {
            clusters.push_back(currentCluster);
        }

        return clusters;
    }

    // we export clusers to file so we can plot them in python, TODO in c++
    void exportClusters(const vector<vector<int>>& clusters, const string& filename) {
        ofstream file(filename);
        // we do 2 dimensions for siplicity
        file << "x,y,cluster\n";
        vector<Point>& points = *m_points;

        for (int clusterId = 0; clusterId < clusters.size(); ++clusterId) {
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


int main(int argc, char* argv[]) {
    vector<Point>* points = loadCsv(string(argv[1]));
    vector<int> processingOrder;
    Optics* model = new Optics(atof(argv[2]), atoi(argv[3]));


    double clusterDistanceThreshold = atof(argv[4]);

    model->fit(*points);

    ofstream outFile("reachability_plot.csv");
    outFile << "Index,ReachabilityDist\n";
    for (int idx : *model->m_processingOrder) {
        outFile << idx << ",";
        if ((*points)[idx].reachabilityDist == UNDEFINED)
            outFile << "UNDEFINED\n";
        else
            outFile << (*points)[idx].reachabilityDist << "\n";
    }
    outFile.close();

    vector<vector<int>> clusters = model->extractClusters(clusterDistanceThreshold);

    model->exportClusters(clusters, "clusters.csv");

    delete model;

    return 0;
}

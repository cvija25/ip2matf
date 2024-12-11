#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <queue>
#include <algorithm>

using DataFrame = std::vector<std::vector<double>>;

DataFrame loadCsv(std::string filename) {
    std::ifstream file(filename);
    DataFrame data;
    std::string line, cell;
    // first line
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::vector<double> row;
        while (std::getline(lineStream, cell, ',')) {
            row.push_back(stod(cell));
        }
        data.push_back(row);
    }
    return data;
}

double dist(std::vector<double> p1, std::vector<double> p2) {
    double suma = 0.0;
    for (int i = 0; i<p1.size(); i++) {
        suma += (p1[i]-p2[i])*(p1[i]-p2[i]);
    }
    return sqrt(suma);
}

class OpticsModel {
    double m_eps;
    int m_minPts;
    // {pointId, reachaiblity}
    std::vector<std::pair<int, double>> m_reachabilityPlot;
public:
    OpticsModel(double eps, int minPts) : m_eps(eps), m_minPts(minPts) {}

    std::vector<std::pair<int, double>> reachabilityPlot() { return m_reachabilityPlot; }

    void fit(DataFrame& data) {
        // preprocessing..
        std::vector<double> coreDist(data.size()), reachabilites(data.size(),std::numeric_limits<int>::max());
        for (int i = 0; i<data.size(); i++) {
            std::priority_queue<double> pq;
            for (int j = 0; j<data.size(); j++) {
                if (i == j) continue;
                double distance = dist(data[i], data[j]);
                if (distance > m_eps) continue;
                pq.push(distance);
            }

            int minPts = m_minPts-1;
            if (pq.size() > minPts) {
                while (minPts--) pq.pop();
                coreDist[i] = pq.top();
            }
        }

        // start algorithm..
        int startingPoint = 7;
        std::priority_queue<std::pair<double, int>> pq;
        std::vector<int> order;
        pq.push({0.0, startingPoint});
        while (!pq.empty()) {
            auto current = pq.top();
            double pointReach = current.first;
            int point = current.second;
            order.push_back(point);
            pq.pop();
            for (int j = 0; j<data.size(); j++) {
                if (startingPoint == j) continue;
                reachabilites[j] = std::min(reachabilites[j], std::max(coreDist[point], dist(data[point], data[j])));
            }
            pq.push({reachabilites[point], point});
        }

        for (int point : order) {
            m_reachabilityPlot.push_back({point, reachabilites[point]});
        }
    }

};

int main() {
    DataFrame data = loadCsv("C:\\Users\\l_cvijic\\Desktop\\personal\\ip2matf\\klasterovanje.csv");
    OpticsModel model = OpticsModel(1.0, 4);
    model.fit(data);
    for (auto row : model.reachabilityPlot()) {
        std::cout << row.second << " " << std::endl;
    }
    return 0;
}

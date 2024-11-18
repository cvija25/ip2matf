#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using Attribute = std::vector<std::string>;
using DataFrame = std::vector<Attribute>;

DataFrame loadCsv(std::string filename) {
    std::ifstream file(filename);
    DataFrame data;
    std::string line, cell;
    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::vector<std::string> row;
        while (std::getline(lineStream, cell, ',')) {
            row.push_back(cell);
        }
        // for (auto dd:row)
        //     std::cout<<dd;
        // std::cout<<std::endl;
        data.push_back(row);
    }
    return data;
}

class kmeansModel {
    std::vector<int> m_labels;
    int m_k;
public:
    kmeansModel(int k) : m_k(k) {}

    void dist(Attribute a, Attribute b) {
        double suma = 0;
        for (int i = 0; i<a.size(); i++) {
            suma += (b[i] - a[i])*(b[i] - a[i])
        }
    }

    void fit(DataFrame data, int generations) {
        std::vector<Attribute> centroids(m_k);
        while (generations--) {
            std::vector<double> sums, counts;
            for (int i = 0; i<data.size();i++) {
                int mini = 0;
                for (int j = 0; j<centroids.size();j++) {
                    if (dist(centroids[j],data[i]) < dist(centroids[mini], data[i])) {
                        mini = j;
                    }
                }
                m_labels[i] = mini
                sums[mini] += data[i];
                counts[mini]++;
            }
        }
    }
};

int main() {
    DataFrame data = loadCsv("klasterovanje.csv");
    for (auto row : data) {
        for (auto cell : row) {
            std::cout << cell;
        }
        std::cout << std::endl;
    }
    return 0;
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <queue>
#include <limits>
#include <algorithm>
#include <unordered_set>

using namespace std;

// Konstanta za nedefinisanu vrednost (beskonacnost)
const double UNDEFINED = numeric_limits<double>::max();

// Struktura za predstavljanje tačke u prostoru
struct Point {
    vector<double> coords;
    bool processed = false;
    double reachabilityDist = UNDEFINED;
};

// Računanje Pearson distance (1 - Pearson korelacija)
// Koristi se za merenje udaljenosti između dve tačke na osnovu njihove korelacije
double pearsonDistance(const Point& p1, const Point& p2) {
    const auto& x = p1.coords;
    const auto& y = p2.coords;
    int n = x.size();

    double meanX = accumulate(x.begin(), x.end(), 0.0) / n;
    double meanY = accumulate(y.begin(), y.end(), 0.0) / n;

    // Racunanje koeficijenta korelacije
    double num = 0.0, denX = 0.0, denY = 0.0;
    for (int i = 0; i < n; ++i) {
        double dx = x[i] - meanX;
        double dy = y[i] - meanY;
        num += dx * dy;      // Brojilac
        denX += dx * dx;     // Imenilac za X
        denY += dy * dy;     // Imenilac za Y
    }

    // Pearson korelacija
    double corr = (denX == 0 || denY == 0) ? 0.0 : num / sqrt(denX * denY);
    
    return 1.0 - corr;
}

vector<Point> loadCsv(const string& filename) {
    ifstream file(filename);
    vector<Point> data;
    string line, cell;
    
    // Preskacemo header red
    getline(file, line);

    // Citamo red po red
    while (getline(file, line)) {
        stringstream lineStream(line);
        vector<double> row;
        
        // Parsiramo svaki element odvojen zarezom
        while (getline(lineStream, cell, ',')) {
            row.push_back(stod(cell));
        }
        data.push_back({row});
    }
    return data;
}

class Optics {
    vector<Point>& points;
    vector<int> processingOrder;
    double epsilon;
    int minPts;

public:
    Optics(vector<Point>& pts, double e, int mP) 
        : points(pts), epsilon(e), minPts(mP) {}

    // Pronalaženje svih suseda tačke koji su u epsilon rastojanju
    vector<pair<int, double>> getNeighbors(int idx) {
        vector<pair<int, double>> neighbors;
        
        for (int i = 0; i < points.size(); ++i) {
            if (i != idx) {
                double dist = pearsonDistance(points[idx], points[i]);
                
                if (dist <= epsilon) {
                    neighbors.push_back({i, dist});
                }
            }
        }
        
        // Sortiramo po rastojanju (najbliži prvo)
        sort(neighbors.begin(), neighbors.end(), 
             [](const auto& a, const auto& b) { return a.second < b.second; });
        return neighbors;
    }

    // Racunanje core distance - rastojanje do minPts-og najbližeg komšije
    // Ako tačka nema dovoljno suseda, nije "core point"
    double getCoreDistance(const vector<pair<int, double>>& neighbors) {
        if (neighbors.size() < minPts) return UNDEFINED;
    
        return neighbors[minPts - 1].second;
    }

    // Azuriranje reachability distance za susede
    // seeds je priority queue tačaka koje treba obraditi
    void updateReachability(priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>>& seeds, 
                           unordered_set<int>& inSeeds,
                           const vector<pair<int, double>>& neighbors, 
                           double coreDist) {
        
        for (const auto& [neighbor, dist] : neighbors) {
            if (!points[neighbor].processed) {
                // Reachability distance je maksimum od core distance i stvarnog rastojanja
                double newReachDist = max(coreDist, dist);
                
                if (newReachDist < points[neighbor].reachabilityDist) {
                    points[neighbor].reachabilityDist = newReachDist;
                    
                    if (inSeeds.find(neighbor) == inSeeds.end()) {
                        seeds.push({newReachDist, neighbor});
                        inSeeds.insert(neighbor);
                    }
                }
            }
        }
    }

    void fit() {
        for (int i = 0; i < points.size(); ++i) {
            if (points[i].processed) continue;

            points[i].processed = true;
            processingOrder.push_back(i);
            
            auto neighbors = getNeighbors(i);
            double coreDist = getCoreDistance(neighbors);

            // Ako je tacka core point (ima dovoljno suseda)
            if (coreDist != UNDEFINED) {
                // Kreiramo priority queue za ekspanziju klastera
                priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> seeds;
                unordered_set<int> inSeeds;
                updateReachability(seeds, inSeeds, neighbors, coreDist);

                // Obradjujemo te tacke
                while (!seeds.empty()) {
                    auto [reachDist, current] = seeds.top();
                    seeds.pop();
                    inSeeds.erase(current);
                    
                    if (points[current].processed) continue;

                    points[current].processed = true;
                    processingOrder.push_back(current);
                    
                    auto currentNeighbors = getNeighbors(current);
                    double currentCoreDist = getCoreDistance(currentNeighbors);
                    
                    // Ako je i ova tačka core point, nastavljamo ekspanziju
                    if (currentCoreDist != UNDEFINED) {
                        updateReachability(seeds, inSeeds, currentNeighbors, currentCoreDist);
                    }
                }
            }
        }
    }

    // Ekstrakcija klastera na osnovu threshold vrednosti
    // Grupišemo uzastopne tačke čija je reachability distance <= threshold
    vector<vector<int>> extractClusters(double threshold) {
        vector<vector<int>> clusters;
        vector<int> currentCluster;

        // Prolazimo kroz tačke u redosledu kako su obrađene
        for (int idx : processingOrder) {
            // Ako je reachability distance ispod threshold, dodajemo u trenutni klaster
            if (points[idx].reachabilityDist <= threshold) {
                currentCluster.push_back(idx);
            } else {
                // Inače, završavamo trenutni klaster i počinjemo novi
                if (!currentCluster.empty()) {
                    clusters.push_back(move(currentCluster));
                    currentCluster.clear();
                }
            }
        }

        // Dodajemo poslednji klaster ako postoji
        if (!currentCluster.empty()) {
            clusters.push_back(move(currentCluster));
        }

        return clusters;
    }

    // Eksportovanje reachability plot-a u CSV fajl
    // Ovaj grafik se koristi za vizuelizaciju gustine podataka
    void exportReachability(const string& filename) {
        ofstream file(filename);
        file << "Index,ReachabilityDist\n";
        
        // Zapisujemo svaku tačku sa njenom reachability distance
        for (int idx : processingOrder) {
            file << idx << "," 
                 << (points[idx].reachabilityDist == UNDEFINED ? "UNDEFINED" : to_string(points[idx].reachabilityDist)) 
                 << "\n";
        }
    }

    // Eksportovanje rezultata klasterovanja u CSV fajl
    // Format: element (indeks tačke), cluster (ID klastera)
    void exportClusters(const vector<vector<int>>& clusters, const string& filename) {
        ofstream file(filename);
        file << "element,cluster\n";

        // Kreiramo mapiranje: indeks elementa -> ID klastera
        // 0 znači da element nije u nijednom klasteru (šum/noise)
        vector<int> elementToCluster(points.size(), 0);

        // Popunjavamo mapiranje - klasteri su numerisani od 1
        for (int clusterId = 0; clusterId < clusters.size(); ++clusterId) {
            for (int idx : clusters[clusterId]) {
                elementToCluster[idx] = clusterId + 1;
            }
        }

        // Zapisujemo sve elemente sa njihovim klasterima
        for (int i = 0; i < points.size(); ++i) {
            file << i << "," << elementToCluster[i] << "\n";
        }
    }
};

int main(int argc, char* argv[]) {
    // Ucitavanje podataka iz CSV fajla
    vector<Point> points = loadCsv(string(argv[1]));
    
    double epsilon = atof(argv[2]);
    int minPts = atoi(argv[3]);
    double clusterThreshold = atof(argv[4]);

    Optics model(points, epsilon, minPts);
    model.fit();
    
    model.exportReachability("reachability_plot.csv");
    
    vector<vector<int>> clusters = model.extractClusters(clusterThreshold);
    model.exportClusters(clusters, "clusters.csv");

    return 0;
}
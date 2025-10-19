#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <sstream>
#include <cstdlib>
#include <ctime>

struct Point {
    std::vector<double> coords;
};

class CLARA {
private:
    std::vector<Point> data;
    std::vector<int> medoids;
    std::vector<int> assignments;
    int k;
    int sample_size;
    int num_samples;
    
    double pearsonCorrelationDistance(const Point& a, const Point& b) const {
        int n = a.coords.size();
        if (n < 2) return 1.0;
        
        double mean_a = 0.0, mean_b = 0.0;
        for (int i = 0; i < n; ++i) {
            mean_a += a.coords[i];
            mean_b += b.coords[i];
        }
        mean_a /= n;
        mean_b /= n;
        
        double cov = 0.0, var_a = 0.0, var_b = 0.0;
        for (int i = 0; i < n; ++i) {
            double dev_a = a.coords[i] - mean_a;
            double dev_b = b.coords[i] - mean_b;
            cov += dev_a * dev_b;
            var_a += dev_a * dev_a;
            var_b += dev_b * dev_b;
        }
        
        double denom = std::sqrt(var_a * var_b);
        if (denom < 1e-10) return 1.0;
        
        double correlation = cov / denom;
        correlation = std::max(-1.0, std::min(1.0, correlation));
        return 1.0 - correlation;
    }
    
    // Dodeli svaki gen njegovoj najbližoj medoidu (centru klastera)
    // assignments[i] = redni broj klastera kojem pripada gen i
    void buildClusters() {
        assignments.assign(data.size(), 0);
        for (int i = 0; i < (int)data.size(); ++i) {
            double min_dist = std::numeric_limits<double>::max();
            // Pronađi medoidu sa najmanjom distancom
            for (int j = 0; j < k; ++j) {
                double dist = pearsonCorrelationDistance(data[i], data[medoids[j]]);
                if (dist < min_dist) {
                    min_dist = dist;
                    assignments[i] = j; // Dodeli genu j-ti klaster
                }
            }
        }
    }
    
    // Izracunaj ukupnu cenu klastera - zbir rastojanja svakog gena do njegove medoide
    // Manja cena = bolji klasteri
    double computeCost(const std::vector<int>& sample_indices) {
        double cost = 0.0;
        for (int idx : sample_indices) {
            double min_dist = std::numeric_limits<double>::max();
            // Pronađi distancu do najbliže medoide
            for (int m : medoids) {
                double dist = pearsonCorrelationDistance(data[idx], data[m]);
                if (dist < min_dist) {
                    min_dist = dist;
                }
            }
            cost += min_dist;
        }
        return cost;
    }
    
    // Optimizuj medoide zamenom sa drugim genima iz uzorka
    // Ako zamena smanji cenu, cuva se nova medoida, inace se vraca stara
    // Nastavlja dok se više ne moze poboljsati
    void optimizeMedoids(const std::vector<int>& sample_indices) {
        bool improved = true;
        while (improved) {
            improved = false;
            // Za svaku trenutnu medoidu
            for (int i = 0; i < (int)medoids.size(); ++i) {
                int old_medoid = medoids[i];
                // Pokusaj zamenu sa svakim drugim genom iz uzorka
                for (int idx : sample_indices) {
                    if (std::find(medoids.begin(), medoids.end(), idx) != medoids.end()) {
                        continue; // Ako je vec medoida, preskoci
                    }
                    // Zameni privremeno
                    medoids[i] = idx;
                    double new_cost = computeCost(sample_indices);
                    medoids[i] = old_medoid;
                    double old_cost = computeCost(sample_indices);
                    
                    // Ako je nova cena bolja, cuva zamenu
                    if (new_cost < old_cost) {
                        medoids[i] = idx;
                        improved = true;
                        break;
                    }
                }
                if (improved) break;
            }
        }
    }
    
    // Nasumicno odaberi k medoida iz uzorka
    // Osiguraj da nema duplikata
    void initializeMedoids(const std::vector<int>& sample_indices) {
        medoids.clear();
        while ((int)medoids.size() < k) {
            int random_idx = rand() % sample_indices.size();
            int idx = sample_indices[random_idx];
            // Proverite da li je vec u listi
            if (std::find(medoids.begin(), medoids.end(), idx) == medoids.end()) {
                medoids.push_back(idx);
            }
        }
    }
    
public:
    CLARA(int k, int sample_size, int num_samples)
        : k(k), sample_size(sample_size), num_samples(num_samples) {}
    
    bool loadCSV(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            return false;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            Point p;
            std::stringstream ss(line);
            std::string val;
            while (std::getline(ss, val, ',')) {
                try {
                    p.coords.push_back(std::stod(val));
                } catch (...) {
                    continue;
                }
            }
            if (!p.coords.empty()) {
                data.push_back(p);
            }
        }
        
        file.close();
        return !data.empty();
    }
    
    // CLARA glavni algoritam - Clustering Large Applications
    // 1. Nasumično uzorkuj sample_size gena
    // 2. Inicijalizuj k medoida na uzorku
    // 3. Optimizuj medoide na uzorku (brzo, jer je mali uzorak)
    // 4. Ponovi num_samples puta i čuva najbolji rezultat
    // 5. Primeni najbolje medoide na sve gene
    void cluster() {
        if (data.empty()) {
            return;
        }
        
        double best_cost = std::numeric_limits<double>::max();
        std::vector<int> best_medoids;
        
        sample_size = std::min(sample_size, (int)data.size());
        
        // Ponovi num_samples puta sa razlicitim uzorcima
        for (int sample = 0; sample < num_samples; ++sample) {
            // Kreiraj random permutaciju svih indeksa
            std::vector<int> indices(data.size());
            for (int i = 0; i < (int)indices.size(); ++i) {
                indices[i] = i;
            }
            
            for (int i = (int)indices.size() - 1; i > 0; --i) {
                int j = rand() % (i + 1);
                std::swap(indices[i], indices[j]);
            }
            
            // Uzmi prvih sample_size indeksa
            std::vector<int> sample_indices(indices.begin(), indices.begin() + sample_size);
            
            // Inicijalizuj i optimizuj medoide na ovom uzorku
            initializeMedoids(sample_indices);
            optimizeMedoids(sample_indices);
            
            // Ako je ova kombinacija bolja, cuva je
            double cost = computeCost(sample_indices);
            if (cost < best_cost) {
                best_cost = cost;
                best_medoids = medoids;
            }
        }
        
        // Primeni najbolje medoide na sve gene
        medoids = best_medoids;
        buildClusters();
    }
    
    bool saveResultsToCSV(const std::string& output_file) const {
        std::ofstream file(output_file);
        if (!file.is_open()) {
            return false;
        }
        
        file << "element_id,cluster\n";
        for (int i = 0; i < (int)assignments.size(); ++i) {
            file << i << "," << assignments[i] << "\n";
        }
        
        file.close();
        return true;
    }
    
    const std::vector<int>& getAssignments() const {
        return assignments;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 5) {
        return 1;
    }
    
    srand(std::time(0));
    
    std::string filename = argv[1];
    int k = std::stoi(argv[2]);
    int sample_size = std::stoi(argv[3]);
    int num_samples = std::stoi(argv[4]);
    
    CLARA clara(k, sample_size, num_samples);
    
    if (!clara.loadCSV(filename)) {
        return 1;
    }
    
    clara.cluster();
    clara.saveResultsToCSV("clustering_results.csv");
    
    return 0;
}
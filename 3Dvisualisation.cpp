#include <iostream>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;

struct Point3D {
    double x, y, z;
    int cluster;
    
    Point3D(double x, double y, double z, int cluster) 
        : x(x), y(y), z(z), cluster(cluster) {}
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
        int col = 0;
        double val;
        int pos = 0;
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

vector<Point3D> combineData(const vector<vector<double>>& attributes, 
                            const vector<pair<int, int>>& clusters) {
    vector<Point3D> points;
    
    for (int i = 0; i < clusters.size(); i++) {
        int idx = clusters[i].first;
        int cluster = clusters[i].second;
        
        if (idx >= 0 && idx < attributes.size() && attributes[idx].size() >= 3) {
            double x = attributes[idx][0];
            double y = attributes[idx][1];
            double z = attributes[idx][2];
            points.push_back(Point3D(x, y, z, cluster));
        }
    }
    
    return points;
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

void countClusters(const vector<Point3D>& data, 
                   vector<int>& counts, int& maxCluster) {
    for (int i = 0; i < data.size(); i++) {
        if (data[i].cluster >= counts.size()) {
            counts.resize(data[i].cluster + 1, 0);
        }
        counts[data[i].cluster]++;
        if (data[i].cluster > maxCluster) {
            maxCluster = data[i].cluster;
        }
    }
}

void writeHTMLHeader(ofstream& file, int totalPoints, 
                     const vector<int>& counts, int maxCluster) {
    file << "<!DOCTYPE html>\n<html>\n<head>\n";
    file << "    <title>3D Cluster Visualization</title>\n";
    file << "    <style>\n";
    file << "        body { margin: 0; overflow: hidden; background: #1a1a1a; }\n";
    file << "        #info { position: absolute; top: 10px; left: 10px; ";
    file << "color: white; background: rgba(0,0,0,0.7); ";
    file << "padding: 15px; border-radius: 8px; }\n";
    file << "        #info h2 { margin: 0 0 10px 0; }\n";
    file << "        .legend-item { display: flex; align-items: center; margin: 5px 0; }\n";
    file << "        .legend-color { width: 15px; height: 15px; border-radius: 50%; ";
    file << "margin-right: 8px; }\n";
    file << "        canvas { display: block; }\n";
    file << "    </style>\n";
    file << "</head>\n<body>\n";
    file << "    <div id=\"info\">\n";
    file << "        <h2>3D Cluster Visualization</h2>\n";
    file << "        <p>Total: " << totalPoints << " points</p>\n";
    file << "        <hr style=\"border-color: #444; margin: 10px 0;\">\n";
    file << "        <div style=\"font-size: 13px; font-weight: bold;\">Clusters:</div>\n";
    
    for (int i = 0; i <= maxCluster; i++) {
        file << "        <div class='legend-item'>\n";
        file << "            <div class='legend-color' style='background-color: " 
             << getColor(i) << ";'></div>\n";
        file << "            <span>Cluster " << i << " (" << counts[i] << ")</span>\n";
        file << "        </div>\n";
    }
    
    file << "    </div>\n";
}

void writeScriptStart(ofstream& file) {
    file << "    <script src=\"https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js\"></script>\n";
    file << "    <script>\n";
    file << "        const scene = new THREE.Scene();\n";
    file << "        scene.background = new THREE.Color(0x1a1a1a);\n";
    file << "        const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 10000);\n";
    file << "        camera.position.set(50, 50, 50);\n";
    file << "        camera.lookAt(0, 0, 0);\n";
    file << "        const renderer = new THREE.WebGLRenderer({ antialias: true });\n";
    file << "        renderer.setSize(window.innerWidth, window.innerHeight);\n";
    file << "        document.body.appendChild(renderer.domElement);\n";
    file << "        const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);\n";
    file << "        scene.add(ambientLight);\n";
    file << "        const pointLight = new THREE.PointLight(0xffffff, 0.8);\n";
    file << "        pointLight.position.set(50, 50, 50);\n";
    file << "        scene.add(pointLight);\n";
    file << "        const axesHelper = new THREE.AxesHelper(100);\n";
    file << "        scene.add(axesHelper);\n";
}

void writeScriptEnd(ofstream& file) {
    file << "        let isDragging = false;\n";
    file << "        let prevX = 0, prevY = 0;\n";
    file << "        renderer.domElement.addEventListener('mousedown', (e) => {\n";
    file << "            isDragging = true; prevX = e.clientX; prevY = e.clientY;\n";
    file << "        });\n";
    file << "        renderer.domElement.addEventListener('mousemove', (e) => {\n";
    file << "            if (isDragging) {\n";
    file << "                scene.rotation.y += (e.clientX - prevX) * 0.005;\n";
    file << "                scene.rotation.x += (e.clientY - prevY) * 0.005;\n";
    file << "                prevX = e.clientX; prevY = e.clientY;\n";
    file << "            }\n";
    file << "        });\n";
    file << "        document.addEventListener('mouseup', () => { isDragging = false; });\n";
    file << "        renderer.domElement.addEventListener('wheel', (e) => {\n";
    file << "            e.preventDefault();\n";
    file << "            const delta = e.deltaY > 0 ? 1.1 : 0.9;\n";
    file << "            camera.position.multiplyScalar(delta);\n";
    file << "        });\n";
    file << "        function animate() {\n";
    file << "            requestAnimationFrame(animate);\n";
    file << "            if (!isDragging) scene.rotation.y += 0.0005;\n";
    file << "            renderer.render(scene, camera);\n";
    file << "        }\n";
    file << "        window.addEventListener('resize', () => {\n";
    file << "            camera.aspect = window.innerWidth / window.innerHeight;\n";
    file << "            camera.updateProjectionMatrix();\n";
    file << "            renderer.setSize(window.innerWidth, window.innerHeight);\n";
    file << "        });\n";
    file << "        animate();\n";
    file << "    </script>\n";
    file << "</body>\n</html>";
}

void plotPoints(ofstream& file, const vector<Point3D>& data) {
    for (int i = 0; i < data.size(); i++) {
        const Point3D& p = data[i];
        file << "        {\n";
        file << "            const geom = new THREE.SphereGeometry(0.5, 8, 8);\n";
        file << "            const mat = new THREE.MeshPhongMaterial({ color: '" 
             << getColor(p.cluster) << "', emissive: '" << getColor(p.cluster) 
             << "', emissiveIntensity: 0.2 });\n";
        file << "            const sphere = new THREE.Mesh(geom, mat);\n";
        file << "            sphere.position.set(" << p.x << ", " << p.y << ", " << p.z << ");\n";
        file << "            scene.add(sphere);\n";
        file << "        }\n";
    }
}

void saveHTML(const vector<Point3D>& data, const string& filename) {
    ofstream file(filename.c_str());
    
    vector<int> counts(10, 0);
    int maxCluster = 0;
    countClusters(data, counts, maxCluster);
    
    writeHTMLHeader(file, data.size(), counts, maxCluster);
    writeScriptStart(file);
    plotPoints(file, data);
    writeScriptEnd(file);
    
    file.close();
}

int main() {
    vector<vector<double>> attributes = loadAttributes("reduced_genes.csv");
    vector<pair<int, int>> clusters = loadClusters("clusters.csv");
    vector<Point3D> data = combineData(attributes, clusters);
    
    if (data.empty()) {
        cerr << "No data combined." << endl;
        return 1;
    }
    
    saveHTML(data, "clusters_3d.html");
    
    return 0;
}
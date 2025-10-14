#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

struct Point3D
{
    double x, y, z;
    int cluster;

    Point3D(double x, double y, double z, int cluster) : x(x), y(y), z(z), cluster(cluster) {}
};

// Dummy dataset - već klasterovani 3D podaci
vector<Point3D> loadDummyData3D()
{
    vector<Point3D> data;

    // Klaster 0 (crveni) - dole levo napred
    data.push_back(Point3D(2.0, 2.0, 2.0, 0));
    data.push_back(Point3D(2.3, 1.8, 2.2, 0));
    data.push_back(Point3D(1.8, 2.2, 1.9, 0));
    data.push_back(Point3D(2.1, 2.1, 2.3, 0));
    data.push_back(Point3D(1.9, 1.9, 2.1, 0));
    data.push_back(Point3D(2.2, 2.3, 1.8, 0));
    data.push_back(Point3D(2.4, 1.7, 2.0, 0));
    data.push_back(Point3D(1.7, 2.4, 2.2, 0));

    // Klaster 1 (plavi) - gore desno nazad
    data.push_back(Point3D(8.0, 8.0, 8.0, 1));
    data.push_back(Point3D(8.3, 7.8, 8.2, 1));
    data.push_back(Point3D(7.8, 8.2, 7.9, 1));
    data.push_back(Point3D(8.1, 8.1, 8.3, 1));
    data.push_back(Point3D(7.9, 7.9, 8.1, 1));
    data.push_back(Point3D(8.2, 8.3, 7.8, 1));
    data.push_back(Point3D(8.4, 7.7, 8.0, 1));
    data.push_back(Point3D(7.7, 8.4, 8.2, 1));

    // Klaster 2 (zeleni) - dole desno napred
    data.push_back(Point3D(8.0, 2.0, 8.0, 2));
    data.push_back(Point3D(8.3, 1.8, 8.2, 2));
    data.push_back(Point3D(7.8, 2.2, 7.9, 2));
    data.push_back(Point3D(8.1, 2.1, 8.3, 2));
    data.push_back(Point3D(7.9, 1.9, 8.1, 2));
    data.push_back(Point3D(8.2, 2.3, 7.8, 2));
    data.push_back(Point3D(8.4, 1.7, 8.0, 2));
    data.push_back(Point3D(7.7, 2.4, 8.2, 2));

    // Klaster 3 (narandžasti) - gore levo nazad
    data.push_back(Point3D(2.0, 8.0, 2.0, 3));
    data.push_back(Point3D(2.3, 7.8, 2.2, 3));
    data.push_back(Point3D(1.8, 8.2, 1.9, 3));
    data.push_back(Point3D(2.1, 8.1, 2.3, 3));
    data.push_back(Point3D(1.9, 7.9, 2.1, 3));
    data.push_back(Point3D(2.2, 8.3, 1.8, 3));
    data.push_back(Point3D(2.4, 7.7, 2.0, 3));
    data.push_back(Point3D(1.7, 8.4, 2.2, 3));

    return data;
}

// Funkcija za učitavanje klastera (preskače header, uzima samo 2. kolonu)
vector<int> loadClusters(const string &filename)
{
    vector<int> clusters;
    ifstream file(filename);
    string line;

    // Preskoči header
    getline(file, line);

    while (getline(file, line))
    {
        stringstream ss(line);
        string index, clusterStr;
        getline(ss, index, ',');      // preskoči prvi stupac (element)
        getline(ss, clusterStr, ','); // drugi stupac = cluster
        clusters.push_back(stoi(clusterStr));
    }
    return clusters;
}

// Funkcija za učitavanje PCA podataka (50 kolona po redu, bez headera)
vector<vector<double>> loadPCA(const string &filename)
{
    vector<vector<double>> data;
    ifstream file(filename);
    string line;
    while (getline(file, line))
    {
        vector<double> row;
        stringstream ss(line);
        string cell;
        while (getline(ss, cell, ','))
        {
            row.push_back(stod(cell));
        }
        data.push_back(row);
    }
    return data;
}

vector<Point3D> combineData(const std::vector<std::vector<double>> &pca, const std::vector<int> &clusters)
{
    vector<Point3D> points;
    for (size_t i = 0; i < pca.size(); i++)
    {
        points.push_back(Point3D(pca[i][0], pca[i][1], pca[i][2], clusters[i]));
    }
    return points;
}

#include <cmath>

void saveToHTML(const vector<Point3D> &data, const string &filename)
{
    ofstream file(filename);

    string colors[] = {"#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8", "#F7DC6F"};

    // Prebroj po klasterima
    vector<int> clusterCounts(10, 0);
    int maxCluster = 0;
    for (const auto &p : data)
    {
        clusterCounts[p.cluster]++;
        maxCluster = max(maxCluster, p.cluster);
    }

    file << R"(<!DOCTYPE html>
<html>
<head>
    <title>3D Vizualizacija klastera</title>
    <style>
        body { 
            margin: 0; 
            overflow: hidden; 
            font-family: Arial, sans-serif;
            background: #1a1a1a;
        }
        #info { 
            position: absolute; 
            top: 10px; 
            left: 10px; 
            color: white; 
            background: rgba(0,0,0,0.7); 
            padding: 15px; 
            border-radius: 8px;
            max-width: 250px;
        }
        #info h2 { margin: 0 0 10px 0; font-size: 20px; }
        #info p { margin: 5px 0; font-size: 14px; }
        .legend-item { 
            display: flex; 
            align-items: center; 
            margin: 5px 0; 
        }
        .legend-color { 
            width: 15px; 
            height: 15px; 
            border-radius: 50%; 
            margin-right: 8px;
            border: 2px solid white;
        }
        canvas { display: block; }
    </style>
</head>
<body>
    <div id="info">
        <h2>3D Vizualizacija klastera</h2>
        <p>Ukupno: )"
         << data.size() << R"( tačaka</p>
        <hr style="border-color: #444; margin: 10px 0;">
        <div style="font-size: 13px; font-weight: bold; margin-bottom: 5px;">Klasteri:</div>
)";

    // Dodaj legendu
    for (int i = 0; i <= maxCluster; i++)
    {
        file << "        <div class='legend-item'>\n";
        file << "            <div class='legend-color' style='background-color: " << colors[i % 6] << ";'></div>\n";
        file << "            <span>Klaster " << i << " (" << clusterCounts[i] << ")</span>\n";
        file << "        </div>\n";
    }

    file << R"(        <hr style="border-color: #444; margin: 10px 0;">
        <p style="font-size: 12px; color: #aaa;">
            <strong>Kontrole:</strong><br>
            🖱️ Prevuci mišem - rotacija<br>
            🔍 Scroll - zum
        </p>
    </div>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script>
        // Setup
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0x1a1a1a);
        
        const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
        camera.position.set(15, 15, 15);
        camera.lookAt(5, 5, 5);
        
        const renderer = new THREE.WebGLRenderer({ antialias: true });
        renderer.setSize(window.innerWidth, window.innerHeight);
        document.body.appendChild(renderer.domElement);
        
        // Osvetljenje
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
        scene.add(ambientLight);
        const pointLight = new THREE.PointLight(0xffffff, 0.8);
        pointLight.position.set(15, 15, 15);
        scene.add(pointLight);
        
        // Ose
        const axesHelper = new THREE.AxesHelper(12);
        scene.add(axesHelper);
        
        // Grid
        const gridHelper = new THREE.GridHelper(20, 20, 0x444444, 0x222222);
        scene.add(gridHelper);
        
        // Dodaj labele za ose
        const canvas = document.createElement('canvas');
        const context = canvas.getContext('2d');
        canvas.width = 64;
        canvas.height = 64;
        
        function createTextTexture(text, color) {
            const canvas = document.createElement('canvas');
            const context = canvas.getContext('2d');
            canvas.width = 64;
            canvas.height = 64;
            context.fillStyle = color;
            context.font = 'Bold 32px Arial';
            context.textAlign = 'center';
            context.fillText(text, 32, 40);
            return new THREE.CanvasTexture(canvas);
        }
        
        // X osa (crvena)
        const xTexture = createTextTexture('X', '#ff0000');
        const xSprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: xTexture }));
        xSprite.position.set(13, 0, 0);
        xSprite.scale.set(1.5, 1.5, 1);
        scene.add(xSprite);
        
        // Y osa (zelena)
        const yTexture = createTextTexture('Y', '#00ff00');
        const ySprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: yTexture }));
        ySprite.position.set(0, 13, 0);
        ySprite.scale.set(1.5, 1.5, 1);
        scene.add(ySprite);
        
        // Z osa (plava)
        const zTexture = createTextTexture('Z', '#0000ff');
        const zSprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: zTexture }));
        zSprite.position.set(0, 0, 13);
        zSprite.scale.set(1.5, 1.5, 1);
        scene.add(zSprite);
)";

    // Dodaj podatke
    file << "\n        // Tačke podataka\n";
    for (size_t i = 0; i < data.size(); i++)
    {
        const auto &p = data[i];
        file << "        {\n";
        file << "            const geom = new THREE.SphereGeometry(0.2, 16, 16);\n";
        file << "            const mat = new THREE.MeshPhongMaterial({ \n";
        file << "                color: '" << colors[p.cluster % 6] << "',\n";
        file << "                emissive: '" << colors[p.cluster % 6] << "',\n";
        file << "                emissiveIntensity: 0.2\n";
        file << "            });\n";
        file << "            const sphere = new THREE.Mesh(geom, mat);\n";
        file << "            sphere.position.set(" << p.x << ", " << p.z << ", " << p.y << ");\n";
        file << "            scene.add(sphere);\n";
        file << "        }\n";
    }

    file << R"(
        // Kontrole misa
        let isDragging = false;
        let previousMousePosition = { x: 0, y: 0 };
        
        renderer.domElement.addEventListener('mousedown', (e) => {
            isDragging = true;
            previousMousePosition = { x: e.clientX, y: e.clientY };
        });
        
        renderer.domElement.addEventListener('mousemove', (e) => {
            if (isDragging) {
                const deltaX = e.clientX - previousMousePosition.x;
                const deltaY = e.clientY - previousMousePosition.y;
                
                scene.rotation.y += deltaX * 0.005;
                scene.rotation.x += deltaY * 0.005;
                
                previousMousePosition = { x: e.clientX, y: e.clientY };
            }
        });
        
        document.addEventListener('mouseup', () => {
            isDragging = false;
        });
        
        // Zoom sa scroll-om
        renderer.domElement.addEventListener('wheel', (e) => {
            e.preventDefault();
            const zoomSpeed = 0.1;
            const delta = e.deltaY > 0 ? 1 + zoomSpeed : 1 - zoomSpeed;
            camera.position.multiplyScalar(delta);
        });
        
        // Animacija
        function animate() {
            requestAnimationFrame(animate);
            // Lagana automatska rotacija kad ne vuces
            if (!isDragging) {
                scene.rotation.y += 0.001;
            }
            renderer.render(scene, camera);
        }
        
        // Resize
        window.addEventListener('resize', () => {
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        });
        
        animate();
    </script>
</body>
</html>)";

    file.close();
    cout << "HTML sacuvan: " << filename << endl;
}

void scaleX(vector<Point3D> &points, double factor = 0.5)
{
    if (points.empty())
        return;

    // Nađi min i max po x
    double minX = points[0].x, maxX = points[0].x;
    for (const auto &p : points)
    {
        if (p.x < minX)
            minX = p.x;
        if (p.x > maxX)
            maxX = p.x;
    }

    double midX = (minX + maxX) / 2.0;
    double halfRange = (maxX - minX) / 2.0 * factor;

    // Skaliraj x oko sredine
    for (auto &p : points)
    {
        p.x = midX + (p.x - midX) * factor;
    }
}

int main()
{
    cout << "=== 3D Vizualizacija klastera ===" << endl;

    vector<vector<double>> pcaData = loadPCA("reduced_genes.csv");
    vector<int> clusterData = loadClusters("clusters.csv");

    cout << pcaData.size() << " PCA tacaka ucitano." << endl;
    cout << clusterData.size() << " klaster ID-jeva ucitano." << endl;

    /*if (pcaData.size() != clusterData.size())
    {
        cerr << "Greska: broj PCA tacaka i klaster ID-jeva se ne poklapa!" << endl;
        return 1;
    }*/

    vector<Point3D> data = combineData(pcaData, clusterData);
    scaleX(data);

    cout << "Ucitano " << data.size() << " 3D tacaka" << endl;

    saveToHTML(data, "clusters_3d_n.html");

    cout << "Gotovo! Otvori 'clusters_3d.html' u browseru." << endl;
    return 0;
}
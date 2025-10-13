#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

struct Point
{
    double x, y;
    int cluster;

    Point(double x, double y, int cluster) : x(x), y(y), cluster(cluster) {}
};

// Dummy dataset - već klasterovani podaci
vector<Point> loadDummyData()
{
    vector<Point> data;

    // Klaster 0 (crveni) - dole levo
    data.push_back(Point(1.5, 2.0, 0));
    data.push_back(Point(2.0, 1.8, 0));
    data.push_back(Point(1.8, 2.3, 0));
    data.push_back(Point(2.2, 2.1, 0));
    data.push_back(Point(1.7, 1.9, 0));
    data.push_back(Point(2.1, 2.4, 0));
    data.push_back(Point(1.9, 2.2, 0));
    data.push_back(Point(2.3, 1.7, 0));

    // Klaster 1 (plavi) - gore desno
    data.push_back(Point(7.5, 8.0, 1));
    data.push_back(Point(8.0, 7.8, 1));
    data.push_back(Point(7.8, 8.3, 1));
    data.push_back(Point(8.2, 8.1, 1));
    data.push_back(Point(7.7, 7.9, 1));
    data.push_back(Point(8.1, 8.4, 1));
    data.push_back(Point(7.9, 8.2, 1));
    data.push_back(Point(8.3, 7.7, 1));

    // Klaster 2 (zeleni) - dole desno
    data.push_back(Point(7.5, 2.0, 2));
    data.push_back(Point(8.0, 1.8, 2));
    data.push_back(Point(7.8, 2.3, 2));
    data.push_back(Point(8.2, 2.1, 2));
    data.push_back(Point(7.7, 1.9, 2));
    data.push_back(Point(8.1, 2.4, 2));
    data.push_back(Point(7.9, 2.2, 2));
    data.push_back(Point(8.3, 1.7, 2));

    return data;
}

void saveToSVG(const vector<Point> &data, const string &filename)
{
    ofstream file(filename);

    int width = 800, height = 600;
    int margin = 50;

    // Pronađi min/max za skaliranje
    double minX = data[0].x, maxX = data[0].x;
    double minY = data[0].y, maxY = data[0].y;

    for (const auto &p : data)
    {
        minX = min(minX, p.x);
        maxX = max(maxX, p.x);
        minY = min(minY, p.y);
        maxY = max(maxY, p.y);
    }

    auto scale = [&](double x, double y, double &sx, double &sy)
    {
        sx = margin + (x - minX) / (maxX - minX) * (width - 2 * margin);
        sy = height - margin - (y - minY) / (maxY - minY) * (height - 2 * margin);
    };

    string colors[] = {"#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8", "#F7DC6F"};
    string clusterNames[] = {"Crveni", "Plavi", "Zeleni", "Narandzasti", "Tirkiz", "Zuti"};

    file << "<svg width='" << width << "' height='" << height << "' xmlns='http://www.w3.org/2000/svg'>\n";
    file << "<rect width='100%' height='100%' fill='#f8f9fa'/>\n";

    // Naslov
    file << "<text x='" << width / 2 << "' y='30' text-anchor='middle' font-size='24' font-weight='bold'>2D Vizualizacija klastera</text>\n";

    // Nacrtaj podatke
    for (const auto &p : data)
    {
        double sx, sy;
        scale(p.x, p.y, sx, sy);
        file << "<circle cx='" << sx << "' cy='" << sy
             << "' r='8' fill='" << colors[p.cluster % 6]
             << "' stroke='#333' stroke-width='1.5' opacity='0.85'/>\n";
    }

    // Dodaj legendu
    int legendY = height - 100;
    file << "<text x='20' y='" << legendY - 10 << "' font-size='16' font-weight='bold'>Legenda:</text>\n";

    // Prebroj koliko ima tačaka po klasteru
    vector<int> clusterCounts(6, 0);
    int maxCluster = 0;
    for (const auto &p : data)
    {
        clusterCounts[p.cluster]++;
        maxCluster = max(maxCluster, p.cluster);
    }

    for (int i = 0; i <= maxCluster; i++)
    {
        int yPos = legendY + i * 25;
        file << "<circle cx='30' cy='" << yPos << "' r='6' fill='" << colors[i % 6] << "' stroke='#333' stroke-width='1'/>\n";
        file << "<text x='45' y='" << yPos + 5 << "' font-size='14'>Klaster " << i
             << " (" << clusterCounts[i] << " tacaka)</text>\n";
    }

    file << "</svg>";

    file.close();
    cout << "SVG sacuvan: " << filename << endl;
}

int main()
{
    cout << "=== 2D Vizualizacija klastera ===" << endl;

    // Ucitaj dummy podatke
    vector<Point> data = loadDummyData();
    cout << "Ucitano " << data.size() << " tacaka" << endl;

    // Prebroj po klasterima
    vector<int> clusterCounts(10, 0);
    for (const auto &p : data)
    {
        clusterCounts[p.cluster]++;
    }

    cout << "\nRaspodela po klasterima:" << endl;
    for (int i = 0; i < 10; i++)
    {
        if (clusterCounts[i] > 0)
        {
            cout << "  Klaster " << i << ": " << clusterCounts[i] << " tacaka" << endl;
        }
    }

    // Sacuvaj vizualizaciju
    saveToSVG(data, "clusters_2d.svg");

    cout << "\nGotovo! Otvori 'clusters_2d.svg' u browseru." << endl;

    return 0;
}
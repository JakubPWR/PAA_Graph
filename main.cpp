#include <iostream>
#include <vector>
#include <queue>
#include <ctime>
#include <cstdlib>
#include <chrono>

using namespace std;
using namespace std::chrono;

// Klasa reprezentująca graf za pomocą macierzy sąsiedztwa
class MatrixGraph {
private:
    int V;
    vector<vector<int>> graph;

public:
    MatrixGraph(int vertices) : V(vertices) {
        graph.resize(V, vector<int>(V, 0)); // Inicjalizacja macierzy zerami
    }

    // Dodaj krawędź z wagą
    void addEdge(int u, int v, int weight) {
        graph[u][v] = weight;
        graph[v][u] = weight; // Graf nieskierowany, więc dodajemy krawędź w drugą stronę
    }

    // Algorytm Dijkstry
    vector<int> dijkstra(int src) {
        vector<int> dist(V, INT_MAX);
        vector<bool> visited(V, false);
        dist[src] = 0;

        for (int count = 0; count < V - 1; ++count) {
            int u = -1;
            // Znajdź wierzchołek o najmniejszej odległości, który nie został jeszcze odwiedzony
            for (int i = 0; i < V; ++i) {
                if (!visited[i] && (u == -1 || dist[i] < dist[u]))
                    u = i;
            }

            visited[u] = true;

            // Zaktualizuj odległości sąsiadów wierzchołka u
            for (int v = 0; v < V; ++v) {
                if (!visited[v] && graph[u][v] && dist[u] != INT_MAX && dist[u] + graph[u][v] < dist[v]) {
                    dist[v] = dist[u] + graph[u][v];
                }
            }
        }

        return dist;
    }
};

// Klasa reprezentująca graf za pomocą listy sąsiedztwa
class ListGraph {
private:
    int V;
    vector<vector<pair<int, int>>> adjList;

public:
    ListGraph(int vertices) : V(vertices) {
        adjList.resize(V);
    }

    // Dodaj krawędź z wagą
    void addEdge(int u, int v, int weight) {
        adjList[u].push_back({v, weight});
        adjList[v].push_back({u, weight}); // Graf nieskierowany, więc dodajemy krawędź w drugą stronę
    }

    // Algorytm Dijkstry
    vector<int> dijkstra(int src) {
        vector<int> dist(V, INT_MAX);
        vector<bool> visited(V, false);
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
        dist[src] = 0;
        pq.push({0, src});

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

            if (visited[u]) continue;
            visited[u] = true;

            // Zaktualizuj odległości sąsiadów wierzchołka u
            for (auto neighbor : adjList[u]) {
                int v = neighbor.first;
                int weight = neighbor.second;
                if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    pq.push({dist[v], v});
                }
            }
        }

        return dist;
    }
};

// Funkcja generująca losowy graf o zadanej liczbie wierzchołków i gęstości
void generateRandomGraph(MatrixGraph &matrixGraph, ListGraph &listGraph, int numVertices, double density) {
    for (int u = 0; u < numVertices; ++u) {
        for (int v = u + 1; v < numVertices; ++v) {
            double prob = (double)rand() / RAND_MAX; // Losowa liczba z przedziału [0, 1]
            if (prob < density) { // Jeśli losowa liczba jest mniejsza niż gęstość, dodaj krawędź
                int weight = rand() % 10 + 1; // Losowa waga z przedziału [1, 10]
                matrixGraph.addEdge(u, v, weight);
                listGraph.addEdge(u, v, weight);
            }
        }
    }
}

// Funkcja mierząca czas wykonania algorytmu Dijkstry dla macierzy sąsiedztwa
double measureTimeMatrix(MatrixGraph &graph, int src) {
    auto start = high_resolution_clock::now();
    graph.dijkstra(src);
    auto end = high_resolution_clock::now();
    duration<double> duration = end - start;
    return duration.count(); // Zwróć czas wykonania w sekundach
}

// Funkcja mierząca czas wykonania algorytmu Dijkstry dla listy sąsiedztwa
double measureTimeList(ListGraph &graph, int src) {
    auto start = high_resolution_clock::now();
    graph.dijkstra(src);
    auto end = high_resolution_clock::now();
    duration<double> duration = end - start;
    return duration.count(); // Zwróć czas wykonania w sekundach
}

int main() {
    srand(time(nullptr)); // Inicjalizacja generatora liczb pseudolosowych

    // Parametry eksperymentu
    vector<int> numVertices = {10, 50, 100, 500, 1000};
    vector<double> densities = {0.25, 0.5, 0.75, 1.0}; // 25%, 50%, 75%, graf pełny

    // Liczba powtórzeń dla każdej kombinacji parametrów
    const int repetitions = 10;

    // Wyniki pomiarów czasu dla macierzy sąsiedztwa
    vector<vector<double>> matrixTimes(numVertices.size(), vector<double>(densities.size(), 0.0));

    // Wyniki pomiarów czasu dla listy sąsiedztwa
    vector<vector<double>> listTimes(numVertices.size(), vector<double>(densities.size(), 0.0));

    // Przeprowadź eksperymenty
    for (int i = 0; i < numVertices.size(); ++i) {
        for (int j = 0; j < densities.size(); ++j) {
            for (int k = 0; k < repetitions; ++k) {
                int numV = numVertices[i];
                double dens = densities[j];

                MatrixGraph matrixGraph(numV);
                ListGraph listGraph(numV);

                generateRandomGraph(matrixGraph, listGraph, numV, dens);

                int src = rand() % numV; // Wylosuj wierzchołek źródłowy

                matrixTimes[i][j] += measureTimeMatrix(matrixGraph, src);
                listTimes[i][j] += measureTimeList(listGraph, src);
            }

            matrixTimes[i][j] /= repetitions; // Uśrednij czasy wykonania
            listTimes[i][j] /= repetitions;
        }
    }

    // Wyświetl wyniki
    cout << "Uśrednione czasy wykonania algorytmu Dijkstry (w sekundach):" << endl;
    cout << "Liczba wierzchołków \\ Gęstość\t25%\t\t50%\t\t75%\t\tPełny graf" << endl;
    for (int i = 0; i < numVertices.size(); ++i) {
        cout << numVertices[i] << "\t\t";
        for (int j = 0; j < densities.size(); ++j) {
            cout << matrixTimes[i][j] << "\t\t";
        }
        cout << endl;
    }

    return 0;
}

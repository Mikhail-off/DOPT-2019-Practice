#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>

typedef long long ull;

struct Point {
    ull x;
    ull y;
};

struct Edge {
    double weight;
    size_t from;
    size_t to;
};

bool operator < (Edge e1, Edge e2) { return e1.weight < e2.weight; }
bool operator == (Edge e1, Edge e2) {
    return e1.from == e2.from && e1.to == e2.to && e1.weight == e2.weight;
}
bool operator != (Edge e1, Edge e2) { return !(e1 == e2); }

double dist(Point p1, Point p2) {
    return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

std::vector<Edge> MST(std::vector<Edge>& edges, size_t vertSize) {
    std::sort(edges.begin(), edges.end());
    std::vector<size_t> subtrees(vertSize);
    std::vector<Edge> MST_edges;
    for (size_t i = 0; i != vertSize; ++i) subtrees[i] = i;
    size_t labels = vertSize;
    for (auto edge : edges) {
        if (labels == 1) break;
        if (subtrees[edge.from] != subtrees[edge.to]) {
            labels--;
            MST_edges.push_back(edge);
            size_t oldLabel = subtrees[edge.to], newLabel = subtrees[edge.from];
            for (size_t j = 0; j != subtrees.size(); ++j) {
                if (subtrees[j] == oldLabel) subtrees[j] = newLabel;
            }
        }
    }
    std::sort(MST_edges.begin(), MST_edges.end());
    return MST_edges;
}

double komi(std::vector<Point> countries, std::vector<size_t>& vertexes) {
    double pathweight = 0;
    for (size_t i = 0; i != countries.size(); ++i)
        pathweight += dist(countries[vertexes[i]], countries[vertexes[(i + 1) % vertexes.size()]]);
    return pathweight;
}


int main() {
    size_t n;
    std::cin >> n;
    /*if (n == 25) {
    std::cout << "ЧТО ЭТО ЗА ТЕСТ ТАКОЙ???????????????????????" << std::endl;
    return 0;
    }*/
    std::vector<Edge> edges;
    std::vector<Point> countries;
    for (size_t i = 0; i != n; ++i) {
        Point country;
        std::cin >> country.x >> country.y;
        for (size_t j = 0; j != countries.size(); ++j) {
            edges.push_back({ dist(countries[j], country), j, i });
        }
        countries.push_back(country);
    }
    auto MST_edges = MST(edges, countries.size());
    std::vector<size_t> ham_path;
    std::vector<bool> visited(countries.size(), false);
    ham_path.push_back(0);
    for (size_t i = 0; i != MST_edges.size(); ++i) {
        size_t last_vertex = ham_path.back();
        visited[ham_path.back()] = true;
        size_t best_candidate = last_vertex;
        double min_weight;
        bool was = false;
        for (auto edge : MST_edges) {
            if (edge.to == last_vertex) std::swap(edge.to, edge.from);
            if (edge.from == last_vertex && visited[edge.to] == false
                && (was == false || edge.weight < min_weight)) {
                min_weight = edge.weight;
                best_candidate = edge.to;
                was = true;
            }
        }
        if (last_vertex != best_candidate) {
            ham_path.push_back(best_candidate);
        } else {
            Edge min_edge;
            bool was_not_initialized = true;
            for (auto edge : edges) {
                if (edge.to == last_vertex) std::swap(edge.to, edge.from);
                if (edge.from == last_vertex && visited[edge.to] == false
                    && (was_not_initialized || edge.weight < min_edge.weight)) {
                    min_edge = edge;
                    was_not_initialized = false;
                }
            }
            if (was_not_initialized == false) ham_path.push_back(min_edge.to);
        }
    }
    ull constant = 50000000;
    if (n < 30) {
        double best_weight = komi(countries, ham_path);
        size_t constant = 750 * static_cast<size_t>(std::sqrt(ham_path.size() + edges.size()));
        for (size_t i = 0
            ; i != constant; ++i) {
            size_t left = std::rand() % ham_path.size(), right = std::rand() % ham_path.size();
            if (left == 0) left++;
            if (right == 0) right++;
            if (left > right) std::swap(left, right);
            if (right - left < 2) continue;
            std::reverse(ham_path.begin() + left, ham_path.begin() + right + 1);
            double current_weight = komi(countries, ham_path);
            if (current_weight < best_weight)
                best_weight = current_weight;
            else
                std::reverse(ham_path.begin() + left, ham_path.begin() + right + 1);
        }
    } else {
        if (n < 2000) constant = 135000000;
        for (size_t i = 0; i != constant; ++i) {
            size_t left = std::rand() % n, right = std::rand() % n;
            if (left == 0) left++;
            if (right == 0) right++;
            if (left > right) std::swap(left, right);
            if (right == left) continue;
            double difference = 0;
            difference += dist(countries[ham_path[left]], countries[ham_path[left - 1]]);
            difference += dist(countries[ham_path[right]], countries[ham_path[(right + 1) % n]]);
            difference -= dist(countries[ham_path[left]], countries[ham_path[(right + 1) % n]]);
            difference -= dist(countries[ham_path[right]], countries[ham_path[left - 1]]);
            if (difference > 0.0001) {
                std::reverse(ham_path.begin() + left, ham_path.begin() + right + 1);
            }
        }
    }
    if (ham_path.size() != 0) ham_path.push_back(0);
    for (auto vertex : ham_path) std::cout << vertex + 1 << ' ';
    return 0;
}


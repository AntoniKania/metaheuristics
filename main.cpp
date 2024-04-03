#include <iostream>
#include <vector>
#include <random>
#include <functional>

namespace mhe {
    std::random_device rd;
    std::mt19937 rdgen(rd());

    using adjacency_matrix_t = std::vector<char>;

    int count_nodes_in_graph(const adjacency_matrix_t &graph) {
        return static_cast<int>(sqrt(static_cast<int>(graph.size()) * 8 + 1) + 1) / 2;
    }

    int nth_triangle(int n){
        return (n * n + n) / 2;
    }

    void generate_graphviz_output(const adjacency_matrix_t &adjacency_matrix) {
        int node_count = count_nodes_in_graph(adjacency_matrix);
        std::cout << "graph {" << std::endl;
        for (int i = 0; i < node_count - 1; i++) {
            for (int j = i; j < node_count - 1; j++) {
                int array_index = (node_count - 1) * i - nth_triangle(i - 1) + j - i;
                if (adjacency_matrix[array_index] == 1) {
                    std::cout << "    " << i << " -- " << j + 1 << ";" << std::endl;
                }
            }
        }
        std::cout << "}" << std::endl;
    }

    adjacency_matrix_t generate_random_graph(const int &num_nodes) {
        int adjacency_matrix_size = num_nodes * (num_nodes - 1) / 2;
        adjacency_matrix_t adjacency_matrix(adjacency_matrix_size);
        std::uniform_int_distribution<char> dist(0, 1);
        for (auto &e : adjacency_matrix) e = dist(rdgen);
        return adjacency_matrix;
    }
}

int main() {
    using namespace mhe;
    adjacency_matrix_t graph = {
            1,1,1,1,1,
            1,1,1,1,
            1,1,1,
            1,1,
            1
    };

    generate_graphviz_output(generate_random_graph(8));

    return 0;
}

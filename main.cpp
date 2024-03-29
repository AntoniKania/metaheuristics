#include <iostream>
#include <vector>
#include <random>

namespace mhe {
    using adjacency_matrix = std::vector<char>;

    int count_nodes_in_graph(const adjacency_matrix &graph) {
        return static_cast<int>(sqrt(static_cast<int>(graph.size()) * 8 + 1) + 1) / 2;
    }

    int nth_triangle(int n){
        return (n * n + n) / 2;
    }

    void generate_graphviz_output(const adjacency_matrix &adjacency_matrix) {
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
}

int main() {
    using namespace mhe;
    adjacency_matrix graph = {
            1,1,1,1,1,
            1,1,1,1,
            1,1,1,
            1,1,
            1
    };

    generate_graphviz_output(graph);

    return 0;
}

#include <iostream>
#include <vector>
#include <random>
#include <functional>

namespace mhe {
    std::random_device rd;
    std::mt19937 rdgen(rd());

    using adjacency_matrix_t = std::vector<char>;
    using subgraph_t = std::vector<char>;

    // LLM: I have a formula: E = n(n-1)/2. How can I transform it to calculate n?
    int count_nodes_in_graph(const adjacency_matrix_t &graph) {
        return static_cast<int>(sqrt(static_cast<int>(graph.size()) * 8 + 1) + 1) / 2;
    }

    int nth_triangle(int n) {
        return (n * n + n) / 2;
    }

    // node_y should be always < than node_x
    int get_index_in_adjacency_matrix(int node_y, int node_x, const adjacency_matrix_t &problem) {
        int node_y_without_shift = node_y - 1;
        int node_x_without_shift = node_x - 2;
        int adjacencies_missing_from_previous_rows = nth_triangle(node_y_without_shift);
        int base_row_index = node_y_without_shift * (count_nodes_in_graph(problem) - 1);
        return base_row_index - adjacencies_missing_from_previous_rows + node_x_without_shift;
    }

    void generate_graphviz_output(const adjacency_matrix_t &adjacency_matrix) {
        int node_count = count_nodes_in_graph(adjacency_matrix);
        std::cout << "graph {" << std::endl;
        for (int nodeY = 1; nodeY < node_count; nodeY++) {
            for (int nodeX = nodeY + 1; nodeX <= node_count; nodeX++) {
                int array_index = get_index_in_adjacency_matrix(nodeY, nodeX, adjacency_matrix);
                if (adjacency_matrix[array_index] == 1) {
                    std::cout << "    " << nodeY << " -- " << nodeX << ";" << std::endl;
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

    subgraph_t generate_random_subgraph(const adjacency_matrix_t &problem) {
        int node_count = count_nodes_in_graph(problem);
        subgraph_t subgraph = subgraph_t(node_count);
        std::uniform_int_distribution<char> dist(0, 1);
        for (auto &node : subgraph) node = dist(rdgen);
        return subgraph;
    }

    subgraph_t generate_first_subgraph(const adjacency_matrix_t &problem) {
        subgraph_t subgraph(count_nodes_in_graph(problem));
        return subgraph;
    }

    subgraph_t generate_next_solution(subgraph_t subgraph) {
        for (int i = 0; i < subgraph.size(); i++) {
            if (subgraph[i] == 0) {
                subgraph[i] = 1;
                return subgraph;
            }
            subgraph[i] = 0;
        }
        return subgraph;
    }

    adjacency_matrix_t create_subgraph_adjacency_matrix(const subgraph_t &subgraph, const adjacency_matrix_t &problem) {
        adjacency_matrix_t subgraph_adjacency_matrix;
        for (int node_y = 1; node_y < subgraph.size(); node_y++) {
            if (subgraph.at(node_y - 1) == 0) {
                continue; // skipping excluded nodes
            }
            for (int node_x = node_y + 1; node_x <= subgraph.size(); node_x++) {
                if (subgraph.at(node_x - 1) == 0) { // +1 as self edges were removed from matrix, and it's shifted by +1 on x-axis
                    continue; // skipping excluded nodes
                }
                int array_index = get_index_in_adjacency_matrix(node_y, node_x, problem);
                subgraph_adjacency_matrix.push_back(problem.at(array_index));
            }
        }
        return subgraph_adjacency_matrix;
    }

    std::function<int(const subgraph_t &)> goal_factory(const adjacency_matrix_t &problem) {
        return [=](const subgraph_t &subgraph) {
            if (subgraph.size() != count_nodes_in_graph(problem)) throw std::invalid_argument("Bad subgraph size!");
            int score = 0;
            adjacency_matrix_t subgraph_adjacency_matrix = create_subgraph_adjacency_matrix(subgraph, problem);
            for (char i : subgraph_adjacency_matrix) {
                if (i == 1) {
                    score++;
                } else {
                    score -= 100;
                }
            }
            return score;
        };
    }

    subgraph_t solve(const adjacency_matrix_t &problem) {
        auto packing = generate_first_subgraph(problem);
        auto goal = goal_factory(problem);
        auto best_goal_value = goal(packing);
        auto best_solution = packing;
        while (true) {
            packing = generate_next_solution(packing);
            double next_goal_value = goal(packing);
            if (next_goal_value > best_goal_value) {
                best_goal_value = next_goal_value;
                best_solution = packing;
            }
            int s = 0;
            for (auto e : packing) s += e;
            if (s == 0) break;
        }
        return best_solution;
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
    generate_graphviz_output(graph);

//    generate_graphviz_output(generate_random_graph(8));
    auto solution = solve(graph);

    return 0;
}

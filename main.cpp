#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <algorithm>
#include <list>
#include <set>
#include <fstream>
#include <map>
#include <chrono>

std::fstream logger("mhe.log");

namespace mhe {
    std::random_device rd;
    std::mt19937 rdgen(rd());

    using adjacency_matrix_t = std::vector<char>;
    using subgraph_t = std::vector<char>;

    std::ostream & operator<<(std::ostream &o, const subgraph_t &p) {
        for (auto &e: p) {
            o << ((e == 0) ? "0" : "1");
        }
        return o;
    }

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
        logger << "graph {" << std::endl;
        for (int nodeY = 1; nodeY < node_count; nodeY++) {
            for (int nodeX = nodeY + 1; nodeX <= node_count; nodeX++) {
                int array_index = get_index_in_adjacency_matrix(nodeY, nodeX, adjacency_matrix);
                if (adjacency_matrix[array_index] == 1) {
                    logger << "    " << nodeY << " -- " << nodeX << ";" << std::endl;
                }
            }
        }
        logger << "}" << std::endl;
    }

    void generate_graphviz_output(const adjacency_matrix_t &adjacency_matrix, const subgraph_t subgraph) {
        int node_count = count_nodes_in_graph(adjacency_matrix);
        std::vector<int> nodes_in_subgraph;
        for (int i = 1; i <= subgraph.size(); i++) {
            if (subgraph.at(i - 1) == 1) {
                nodes_in_subgraph.push_back(i);
            }
        }
        logger << "graph {" << std::endl;
        for (int nodeY = 1; nodeY < node_count; nodeY++) {
            for (int nodeX = nodeY + 1; nodeX <= node_count; nodeX++) {
                int array_index = get_index_in_adjacency_matrix(nodeY, nodeX, adjacency_matrix);
                if (adjacency_matrix[array_index] == 1) {
                    logger << "    " << nodes_in_subgraph.at(nodeY - 1) << " -- " << nodes_in_subgraph.at(nodeX - 1) << ";" << std::endl;
                }
            }
        }
        logger << "}" << std::endl;
    }

    adjacency_matrix_t generate_random_graph(const int &num_nodes) {
        int adjacency_matrix_size = num_nodes * (num_nodes - 1) / 2; // upper triangle, self connections skipped
        adjacency_matrix_t adjacency_matrix(adjacency_matrix_size);
        std::uniform_int_distribution<char> dist(0, 1);
        for (auto &e : adjacency_matrix) e = dist(rdgen);
        return adjacency_matrix;
    }

    subgraph_t generate_random_subgraph(const adjacency_matrix_t &problem, double p_1 = 0.5) {
        int node_count = count_nodes_in_graph(problem);
        subgraph_t subgraph = subgraph_t(node_count);
        std::uniform_real_distribution<double> dist(0.0,1.0);
        for (auto &e : subgraph) e = (dist(rdgen) < p_1) ? 1 : 0;
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

    std::vector<subgraph_t> generate_neighbor_solutions(const subgraph_t & subgraph) {
        std::vector<subgraph_t> result;
        for (int i = 0; i < subgraph.size(); i++) {
            auto modified = subgraph;
            modified[i] = 1 - modified[i];
            result.push_back(modified);
        }
        return result;
    }

    subgraph_t generate_random_neighbor_norm(subgraph_t p) {
        std::normal_distribution<double> distr;

        int count = distr(rdgen) + 1;
        if (count >= (p.size() * 2)) count = p.size();
        for (int i = 0; i < count; i++) {
            std::uniform_int_distribution<int> selbit(0,p.size() - 1);
            int sel = selbit(rdgen);
            p[sel] = 1 - p[sel];
        }

        return p;
    }

    adjacency_matrix_t create_subgraph_adjacency_matrix(const subgraph_t &subgraph, const adjacency_matrix_t &problem) {
        adjacency_matrix_t subgraph_adjacency_matrix;
        for (int node_y = 1; node_y < subgraph.size(); node_y++) {
            if (subgraph.at(node_y - 1) == 0) {
                continue; // skipping excluded nodes
            }
            for (int node_x = node_y + 1; node_x <= subgraph.size(); node_x++) {
                if (subgraph.at(node_x - 1) == 0) { // -1 to get the index in array
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
        auto subgraph = generate_first_subgraph(problem);
        auto goal = goal_factory(problem);
        auto best_goal_value = goal(subgraph);
        auto best_solution = subgraph;
        while (true) {
            subgraph = generate_next_solution(subgraph);
            double next_goal_value = goal(subgraph);
            if (next_goal_value > best_goal_value) {
                best_goal_value = next_goal_value;
                best_solution = subgraph;
            }
            int s = 0;
            for (auto e : subgraph) s += e;
            if (s == 0) break;
        }
        return best_solution;
    }

    subgraph_t solve_random(const adjacency_matrix_t &problem, int iterations = 20, double p_1 = 0.5) {
        using namespace std;
        auto packing = generate_random_subgraph(problem, p_1);
        auto goal = goal_factory(problem);

        for ( int i = 0; i <  iterations; i++) {
            auto t = generate_random_subgraph(problem, p_1);
            if (goal(t) >= goal(packing) ) {
                packing = t;
                logger << " --> " << packing << " -> " << goal(packing);
                logger << "  BEST! ";
                logger << std::endl;
            }
        }
        return packing;
    }

    subgraph_t solve_hill_climbing(const adjacency_matrix_t &problem, int iterations = 20) {
        auto subgraph = generate_random_subgraph(problem);
        auto goal = goal_factory(problem);
        auto best_goal_value = goal(subgraph);
        auto best_solution = subgraph;
        for (int i = 0; i < iterations; i++) {
            auto subgraphs = generate_neighbor_solutions(subgraph);
            subgraph = *std::max_element(subgraphs.begin(), subgraphs.end(),
                                         [=](auto a, auto b) {
                                            return goal(a) < goal(b);
                                        });
            double next_goal_value = goal(subgraph);
            if (next_goal_value > best_goal_value) {
                best_goal_value = next_goal_value;
                best_solution = subgraph;
                logger << " --> " << subgraph << " -> " <<  best_goal_value;
                logger << "  BEST! ";
                logger << std::endl;
            }
        }
        return best_solution;
    }

    subgraph_t solve_tabu_set(const adjacency_matrix_t &problem, int iterations = 20, int tabu_size = 1000) {
        auto subgraph = generate_random_subgraph(problem);
        auto best_subgraph = subgraph;
        auto goal = goal_factory(problem);

        std::list<subgraph_t> tabu_list = {subgraph};
        std::set<subgraph_t> tabu = {subgraph};

        auto is_in_tabu = [&](subgraph_t &elem) -> bool {
            return tabu.find(elem) != tabu.end() ;
        };

        for (int i = 0; i <  iterations; i++) {
            std::vector<subgraph_t> subgraphs;
            for (auto p: generate_neighbor_solutions(subgraph)) {
                if (!is_in_tabu(p)) subgraphs.push_back(p);
            }

            if (subgraphs.empty()) {
                logger << "SNAKE CONDITION - TABU IS EMPTY!!!!!!!!! ";
                logger << std::endl;
                return best_subgraph;
            }

            auto subgraph_next = *std::max_element(subgraphs.begin(), subgraphs.end(),
                                                   [=](auto a, auto b)
                                                  {
                                                      return goal(a) < goal(b);
                                                  });
            subgraph = subgraph_next;

            if (goal(subgraph) > goal(best_subgraph)) {
                best_subgraph = subgraph;
                logger << " --> " << subgraph_next << " -> " << goal(subgraph_next);
                logger << "  BEST! ";
                logger << std::endl;
            }
            tabu_list.push_back(subgraph);
            tabu.insert(subgraph);
            if (tabu.size() > tabu_size) {
                tabu.erase(tabu_list.front());
                tabu_list.pop_front();
            }

        }
        return best_subgraph;
    }

    subgraph_t solve_tabu_list(const adjacency_matrix_t &problem, int iterations = 20, int tabu_size = 1000) {
        auto subgraph = generate_random_subgraph(problem);
        auto best_subgraph = subgraph;
        auto goal = goal_factory(problem);

        std::list<subgraph_t> tabu_list = {subgraph};

        auto is_in_tabu = [&](subgraph_t &elem) -> bool {
            return std::find(tabu_list.begin(), tabu_list.end(),elem) != tabu_list.end();
        };

        for (int i = 0; i <  iterations; i++) {
            std::vector<subgraph_t> subgraphs;
            for (auto p: generate_neighbor_solutions(subgraph)) {
                if (!is_in_tabu(p)) subgraphs.push_back(p);
            }

            if (subgraphs.empty()) {
                logger << "SNAKE CONDITION - TABU IS EMPTY!!!!!!!!! ";
                logger << std::endl;
                return best_subgraph;
            }

            auto subgraph_next = *std::max_element(subgraphs.begin(), subgraphs.end(),
                                                   [=](auto a, auto b)
                                                   {
                                                       return goal(a) < goal(b);
                                                   });
            subgraph = subgraph_next;

            if (goal(subgraph) > goal(best_subgraph)) {
                best_subgraph = subgraph;
                logger << " --> " << subgraph_next << " -> " << goal(subgraph_next);
                logger << "  BEST! ";
                logger << std::endl;
            }
            tabu_list.push_back(subgraph);
            if (tabu_list.size() > tabu_size) {
                tabu_list.pop_front();
            }

        }
        return best_subgraph;
    }

    subgraph_t solve_tabu_avoid_snake(const adjacency_matrix_t &problem, int iterations = 20, int tabu_size = 1000) {
        auto subgraph = generate_random_subgraph(problem);
        auto best_subgraph = subgraph;
        auto goal = goal_factory(problem);

        std::list<subgraph_t> last_visited_list = {subgraph};
        std::set<subgraph_t> tabu = {subgraph};

        auto is_in_tabu = [&](subgraph_t &elem) -> bool {
            return tabu.find(elem) != tabu.end() ;
        };

        for (int i = 0; i <  iterations; i++) {
            std::vector<subgraph_t> subgraphs;
            subgraph_t subgraph_next;
            for (auto p: generate_neighbor_solutions(subgraph)) {
                if (!is_in_tabu(p)) subgraphs.push_back(p);
            }

            if (subgraphs.empty()) {
                subgraph_next = last_visited_list.back();
                last_visited_list.pop_back();
            } else {
                subgraph_next = *std::max_element(subgraphs.begin(), subgraphs.end(),
                                                  [=](auto a, auto b)
                                                  {
                                                      return goal(a) < goal(b);
                                                  });
            }
            subgraph = subgraph_next;

            if (goal(subgraph) > goal(best_subgraph)) {
                best_subgraph = subgraph;
                logger << " --> " << subgraph_next << " -> " << goal(subgraph_next);
                logger << "  BEST! ";
                logger << std::endl;
            }
            last_visited_list.push_back(subgraph);
            tabu.insert(subgraph);
            if (tabu.size() > tabu_size) {
                tabu.erase(last_visited_list.front());
                last_visited_list.pop_front();
            }

        }
        return best_subgraph;
    }

    subgraph_t solve_sim_annealing(const adjacency_matrix_t &problem, int iterations = 20,
                                   std::function<double(int)> T = [](auto i){return 1.0 / (i + 1);}
                                   ) {
        using namespace  std;
        auto subgraph = generate_random_subgraph(problem);
        auto best_subgraph = subgraph;
        auto goal = goal_factory(problem);

        for (int i = 0; i <  iterations; i++) {
            auto t = generate_random_neighbor_norm(subgraph);
            if (goal(t) >= goal(subgraph) ) {
                subgraph = t;
                if (goal(subgraph) > goal(best_subgraph)) {
                    best_subgraph = subgraph;
                    logger << " --> " << best_subgraph << " -> " << goal(best_subgraph);
                    logger << "  BEST! ";
                    logger << std::endl;
                }
            } else {
                uniform_real_distribution<double> u(0.0,1.0);
                if (u(rdgen) < exp((abs(goal(t) - goal(subgraph)) / T(i))) ) {
                    subgraph = t;
                }
            }
        }
        return best_subgraph;
    }

    struct subgraph_with_score {
        subgraph_t subgraph;
        int score;
    };

    void evaluate_population(std::vector<subgraph_with_score>& population,
                             const std::function<int(const subgraph_t)>& goal_function
                             ) {
        for (auto& individual : population) {
            individual.score = goal_function(individual.subgraph);
        }
    }

    std::vector<subgraph_with_score> select_using_truncation_selection(const std::vector<subgraph_with_score>& population) {
        auto sorted_population = population;
        std::sort(sorted_population.begin(), sorted_population.end(),
                  [](const subgraph_with_score& a, const subgraph_with_score& b){ return a.score > b.score; });

        auto first = sorted_population.begin();
        auto last = sorted_population.begin() + population.size() / 2;
        std::vector<subgraph_with_score> selected(first, last);
        return selected;
    }

    subgraph_with_score uniform_crossover(const subgraph_with_score& parent1, const subgraph_with_score& parent2) {
        subgraph_with_score offspring = subgraph_with_score{subgraph_t(parent1.subgraph.size()), 0};
        std::uniform_int_distribution<> dist(0, 1);
        for (int i = 0; i < offspring.subgraph.size(); i++) {
            if (dist(rdgen) == 1) {
                offspring.subgraph[i] = parent1.subgraph[i];
            } else {
                offspring.subgraph[i] = parent2.subgraph[i];
            }
        }

        return offspring;
    }

    subgraph_with_score one_point_crossover(const subgraph_with_score& parent1,
                                            const subgraph_with_score& parent2,
                                            int crossover_point) {
        subgraph_with_score offspring = subgraph_with_score{subgraph_t(parent1.subgraph.size()), 0};

        for (int i = 0; i < crossover_point; i++) {
            offspring.subgraph[i] = parent1.subgraph[i];
        }

        for (int i = crossover_point; i < parent2.subgraph.size(); i++) {
            offspring.subgraph[i] = parent2.subgraph[i];
        }

        return offspring;
    }

    std::vector<subgraph_with_score> perform_one_point_crossover(const std::vector<subgraph_with_score>& selected, int population_size) {
        std::vector<subgraph_with_score> new_population;
        int size_individual = selected.at(0).subgraph.size();
        std::uniform_int_distribution<> dist_parent_index(0, selected.size() - 1);
        std::uniform_int_distribution<char> dist_crossover_point(0, size_individual - 1);
        while (new_population.size() < population_size) {
            int crossover_point = dist_crossover_point(rdgen);
            int parent1_id = dist_parent_index(rdgen);
            int parent2_id = dist_parent_index(rdgen);
            subgraph_with_score offspring1 = one_point_crossover(selected[parent1_id], selected[parent2_id], crossover_point);
            subgraph_with_score offspring2 = one_point_crossover(selected[parent2_id], selected[parent1_id], crossover_point);
            new_population.push_back(offspring1);
            if (new_population.size() < population_size) {
                new_population.push_back(offspring2);
            }
        }
        return new_population;
    }

    std::vector<subgraph_with_score> perform_uniform_crossover(const std::vector<subgraph_with_score>& selected, int population_size) {
        std::vector<subgraph_with_score> new_population;
        std::uniform_int_distribution<> dist_parent_index(0, selected.size() - 1);
        while (new_population.size() < population_size) {
            int parent1_id = dist_parent_index(rdgen);
            int parent2_id = dist_parent_index(rdgen);
            subgraph_with_score offspring = uniform_crossover(selected[parent1_id], selected[parent2_id]);
            new_population.push_back(offspring);
        }
        return new_population;
    }

    void perform_bit_flip_mutation(std::vector<subgraph_with_score>& population) {
        std::uniform_int_distribution<> dist_mutate(0, 1);
        std::uniform_int_distribution<> dist_gene_to_mutate(0, population.at(0).subgraph.size());
        for (auto &e : population) {
            if (dist_mutate(rdgen) == 1) {
                int bit_to_flip = dist_gene_to_mutate(rdgen);
                e.subgraph[bit_to_flip] = 1 - e.subgraph[bit_to_flip];
            }
        }
    }

    void perform_bit_swap_mutation(std::vector<subgraph_with_score>& population) {
        std::uniform_int_distribution<> dist_mutate(0, 1);
        std::uniform_int_distribution<> dist_gene_to_mutate(0, population.at(0).subgraph.size());
        for (auto &e : population) {
            if (dist_mutate(rdgen) == 1) {
                int first_bit = dist_gene_to_mutate(rdgen);
                int second_bit = dist_gene_to_mutate(rdgen);
                char tmp = e.subgraph[first_bit];
                e.subgraph[first_bit] = e.subgraph[second_bit];
                e.subgraph[second_bit] = tmp;
            }
        }
    }

    enum crossover_type {
        one_point,
        uniform
    };

    enum mutation_type {
        bit_flip,
        bit_swap
    };

    subgraph_t solve_genetic_algorithm(const adjacency_matrix_t &problem,
                                       crossover_type crossover_type = one_point,
                                       mutation_type mutation_type = bit_flip,
                                       bool terminate_on_generations_number = false,
                                       int iterations = 1000
                                       ) {
        std::vector<subgraph_with_score> population;
        auto goal = goal_factory(problem);
        for (int i = 0; i < count_nodes_in_graph(problem); i++) {
            const subgraph_t &subgraph = generate_random_subgraph(problem);
            population.push_back(subgraph_with_score{subgraph, 0});
        }

        int max_generations_without_finding_fitter = 0;
        if (!terminate_on_generations_number) {
            max_generations_without_finding_fitter = count_nodes_in_graph(problem) * 20;
        }

        int i = 0;
        int generations_without_fitter = 0;
        int best_score = 0;
        while (true) {
            evaluate_population(population, goal);
            auto selected = select_using_truncation_selection(population);
            std::vector<subgraph_with_score> new_population;
            crossover_type == one_point ?
                    new_population = perform_one_point_crossover(selected, population.size()) :
                    new_population = perform_uniform_crossover(selected, population.size());
            mutation_type == bit_flip ?
                    perform_bit_flip_mutation(new_population) :
                    perform_bit_swap_mutation(new_population);
            population = new_population;

            if (!terminate_on_generations_number) {
                evaluate_population(new_population, goal);
                auto sorted_population = new_population;
                std::sort(sorted_population.begin(), sorted_population.end(),
                          [](const subgraph_with_score& a, const subgraph_with_score& b){ return a.score > b.score; });
                if (sorted_population[0].score > best_score) {
                    best_score = sorted_population[0].score;
                    generations_without_fitter = 0;
                } else {
                    generations_without_fitter++;
                    if (generations_without_fitter >= max_generations_without_finding_fitter) {
                        break;
                    }
                }
            } else {
                i++;
                if (i >= iterations) {
                    break;
                }
            }
        }

        evaluate_population(population, goal);
        auto sorted_population = population;
        std::sort(sorted_population.begin(), sorted_population.end(),
                  [](const subgraph_with_score& a, const subgraph_with_score& b){ return a.score > b.score; });

        return sorted_population.at(0).subgraph;
    }
}

int main(int argc, char **argv) {
    using namespace mhe;
    std::vector<std::string> args(argv, argv+argc);
    std::string  selected_solver = "solve_random";
    if (args.size() >= 2) {
        selected_solver = args[1];
    }
    adjacency_matrix_t example_graph = {
            1,1,1,1,1,
            1,0,1,1,
            1,1,1,
            0,0,
            0
    };

    auto enormous_graph = generate_random_graph(10);
    generate_graphviz_output(enormous_graph);

    auto solution = solve_genetic_algorithm(enormous_graph, uniform, bit_swap);
//    auto solution_hill_climbing = solve_hill_climbing(enormous_graph);
//    auto solution_tabu = solve_tabu_set(enormous_graph, 5000);
//    auto solution_random = solve_random(enormous_graph, 10000, 0.1);
//    auto solution_sim_annealing = solve_sim_annealing(enormous_graph, 10000, [](int i){return 1000*std::pow(0.99,(double)i);});

//    const mhe::adjacency_matrix_t &matrix = create_subgraph_adjacency_matrix(solution, enormous_graph);
//    generate_graphviz_output(matrix);
//    const mhe::adjacency_matrix_t &matrix2 = create_subgraph_adjacency_matrix(solution_random, enormous_graph);
//    generate_graphviz_output(matrix2);
//    const mhe::adjacency_matrix_t &matrix3 = create_subgraph_adjacency_matrix(solution_sim_annealing, enormous_graph);
//    generate_graphviz_output(matrix3, solution_sim_annealing);

    std::map<std::string, std::function<subgraph_t(adjacency_matrix_t)>> solvers;
    solvers["solve_hill_climbing"] = [&](auto knapsack){return solve_hill_climbing(knapsack, 10000);};
    solvers["solve_tabu"] = [&](auto knapsack){return solve_tabu_set(knapsack, 5000);};
    solvers["solve_tabu_list"] = [&](auto knapsack){return solve_tabu_list(knapsack, 5000);};
    solvers["solve_tabu_avoid_snake"] = [&](auto knapsack){return solve_tabu_avoid_snake(knapsack, 5000);};
    solvers["solve_sim_annealing"] = [&](auto knapsack){return solve_sim_annealing(knapsack, 1000000, [](int i){return 1000*std::pow(0.99,(double)i);});};
    solvers["solve_random"] = [&](auto knapsack){return solve_random(knapsack, 10000000);};
    solvers["solve_random_n"] = [&](auto knapsack){return solve_random(knapsack, 100000, 0.1);};

    auto goal = goal_factory(enormous_graph);
    auto start_time = std::chrono::system_clock::now();
    auto result = solvers[selected_solver](enormous_graph);
    auto end_time = std::chrono::system_clock::now();
    auto computation_time =  std::chrono::nanoseconds(end_time - start_time);

    std::cout << selected_solver << " " << goal(result) << " " << computation_time.count() << " " << result << std::endl;

    return 0;
}

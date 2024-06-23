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
std::ofstream conv_log("conv.log");

namespace mhe {
    std::random_device rd;
    std::mt19937 rdgen(rd());

    using adjacency_matrix_t = std::vector<char>;
    using subgraph_t = std::vector<char>;

    adjacency_matrix_t string_to_adjacency_matrix(const std::string &str) {
        adjacency_matrix_t matrix;
        for (char c : str) {
            if (c == '1') {
                matrix.push_back(1);
            } else {
                matrix.push_back(0);
            }
        }
        return matrix;
    }

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

    void generate_graphviz_output(const adjacency_matrix_t &problem, const subgraph_t subgraph) {
        std::vector<int> nodes_in_subgraph;
        for (int i = 1; i <= subgraph.size(); i++) {
            if (subgraph.at(i - 1) == 1) {
                nodes_in_subgraph.push_back(i);
            }
        }
        logger << "graph {" << std::endl;
        for (auto nodeY : nodes_in_subgraph) {
            if (nodeY == nodes_in_subgraph.at(nodes_in_subgraph.size() - 1)) {
                continue;
            }
            for (auto nodeX : nodes_in_subgraph) {
                if (nodeX == nodes_in_subgraph.at(0) || nodeY == nodeX) {
                    continue;
                }
                int array_index = get_index_in_adjacency_matrix(nodeY, nodeX, problem);
                if (problem[array_index] == 1) {
                    logger << "    " << nodeY << " -- " << nodeX << ";" << std::endl;
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
        auto subgraph = generate_random_subgraph(problem, p_1);
        auto goal = goal_factory(problem);
        conv_log << goal(subgraph) << std::endl;

        for (int i = 0; i <  iterations; i++) {
            auto t = generate_random_subgraph(problem, p_1);
            if (goal(t) >= goal(subgraph) ) {
                subgraph = t;
                conv_log << goal(subgraph) << std::endl;
                logger << " --> " << subgraph << " -> " << goal(subgraph);
                logger << "  BEST! ";
                logger << std::endl;
            }
        }
        return subgraph;
    }

    subgraph_t solve_hill_climbing(const adjacency_matrix_t &problem, int iterations = 20) {
        auto subgraph = generate_random_subgraph(problem);
        auto goal = goal_factory(problem);
        auto best_goal_value = goal(subgraph);
        auto best_solution = subgraph;
        conv_log << best_goal_value << std::endl;
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
            conv_log << goal(best_solution) << std::endl;
        }
        return best_solution;
    }

    subgraph_t solve_tabu_set(const adjacency_matrix_t &problem, int iterations = 20, int tabu_size = 1000) {
        auto subgraph = generate_random_subgraph(problem);
        auto best_subgraph = subgraph;
        auto goal = goal_factory(problem);
        std::list<subgraph_t> tabu_list = {subgraph};
        std::set<subgraph_t> tabu = {subgraph};
        conv_log << goal(subgraph) << std::endl;

        auto is_in_tabu = [&](subgraph_t &elem) -> bool {
            return tabu.find(elem) != tabu.end() ;
        };
        for (int i = 0; i <  iterations; i++) {
            std::vector<subgraph_t> subgraphs;
            for (auto p: generate_neighbor_solutions(subgraph)) {
                if (!is_in_tabu(p)) subgraphs.push_back(p);
            }
            if (subgraphs.empty()) {
                logger << "SNAKE CONDITION - ALL NEIGHBOURS ARE IN TABU!!!!!!!!!";
                logger << std::endl;
                return best_subgraph;
            }
            auto subgraph_next = *std::max_element(subgraphs.begin(), subgraphs.end(),
                                                   [=](auto a, auto b)
                                                   {
                                                       return goal(a) < goal(b);
                                                   });
            subgraph = subgraph_next;

            conv_log << goal(subgraph) << std::endl;
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
        conv_log << goal(subgraph) << std::endl;

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
            conv_log << goal(subgraph) << std::endl;

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
        conv_log << goal(subgraph) << std::endl;

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
            conv_log << goal(subgraph) << std::endl;
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
        conv_log << goal(subgraph) << std::endl;

        for (int i = 0; i <  iterations; i++) {
            auto t = generate_random_neighbor_norm(subgraph);
            if (goal(t) >= goal(subgraph) ) {
                subgraph = t;
                conv_log << goal(subgraph) << std::endl;
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
                    conv_log << goal(subgraph) << std::endl;
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
        std::uniform_int_distribution<> dist_gene_to_mutate(0, population.at(0).subgraph.size() - 1);
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

    std::vector<subgraph_with_score> tournament_selection(const std::vector<subgraph_with_score> &population) {
        std::uniform_int_distribution<int> u(0, population.size() - 1);
        std::vector<subgraph_with_score> result;
        for (int i = 0; i < population.size(); i++) {
            int idx1 = u(rdgen);
            int idx2 = u(rdgen);
            if (population.at(idx1).score > population.at(idx2).score) population.at(idx1);
            else result.push_back(population.at(idx2));
        }
        return result;
    }

    std::vector<subgraph_with_score> tournament_selection(const std::vector<subgraph_with_score> &population, int result_size) {
        std::uniform_int_distribution<int> u(0, population.size() - 1);
        std::vector<subgraph_with_score> result;
        for (int i = 0; i < result_size; i++) {
            int idx1 = u(rdgen);
            int idx2 = u(rdgen);
            if (population.at(idx1).score > population.at(idx2).score) population.at(idx1);
            else result.push_back(population.at(idx2));
        }
        return result;
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
        std::function<double(const subgraph_t &)> fitness = [goal](const subgraph_t & specimen){
            double g = goal(specimen);
            if (g < 0) return 0.0;
            return g;
        };

        for (int i = 0; i < count_nodes_in_graph(problem) * count_nodes_in_graph(problem); i++) {
            const subgraph_t &subgraph = generate_random_subgraph(problem);
            population.push_back(subgraph_with_score{subgraph, 0});
        }

        int max_generations_without_finding_fitter = 0;
        if (!terminate_on_generations_number) {
            max_generations_without_finding_fitter = count_nodes_in_graph(problem) * 20;
        }

        auto get_best = [&](){
            std::sort(population.begin(), population.end(),
                      [](const subgraph_with_score& a, const subgraph_with_score& b){ return a.score > b.score; });
            return population.at(0);
        };

        int i = 0;
        int generations_without_fitter = 0;
        int best_score = 0;
        while (true) {
            evaluate_population(population, fitness);
            conv_log << get_best().score << std::endl;
            auto selected = tournament_selection(population);
            std::vector<subgraph_with_score> new_population;
            crossover_type == one_point ?
                    new_population = perform_one_point_crossover(selected, population.size()) :
                    new_population = perform_uniform_crossover(selected, population.size());
            mutation_type == bit_flip ?
            perform_bit_flip_mutation(new_population) :
            perform_bit_swap_mutation(new_population);
            population = new_population;

            if (!terminate_on_generations_number) {
                evaluate_population(population, fitness);
                if (get_best().score > best_score) {
                    best_score = get_best().score;
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
        evaluate_population(population, fitness);
        conv_log << get_best().score << std::endl;

        return get_best().subgraph;
    }

    subgraph_t solve_genetic_algorithm_elit(const adjacency_matrix_t &problem,
                                            crossover_type crossover_type = one_point,
                                            mutation_type mutation_type = bit_flip,
                                            bool terminate_on_generations_number = false,
                                            int iterations = 1000,
                                            int elit_size = 1
    ) {
        std::vector<subgraph_with_score> population;
        auto goal = goal_factory(problem);
        std::function<double(const subgraph_t &)> fitness = [goal](const subgraph_t & specimen){
            double g = goal(specimen);
            if (g < 0) return 0.0;
            return g;
        };

        for (int i = 0; i < count_nodes_in_graph(problem) * count_nodes_in_graph(problem); i++) {
            const subgraph_t &subgraph = generate_random_subgraph(problem);
            population.push_back(subgraph_with_score{subgraph, 0});
        }

        int max_generations_without_finding_fitter = 0;
        if (!terminate_on_generations_number) {
            max_generations_without_finding_fitter = count_nodes_in_graph(problem) * 20;
        }

        auto get_best = [&](){
            std::sort(population.begin(), population.end(),
                      [](const subgraph_with_score& a, const subgraph_with_score& b){ return a.score > b.score; });
            return population.at(0);
        };

        int i = 0;
        int generations_without_fitter = 0;
        int best_score = 0;
        while (true) {
            evaluate_population(population, fitness);
            conv_log << get_best().score << std::endl;
            std::vector<subgraph_with_score> elite;
            for (int i = 0; i < elit_size; i++) {
                elite.push_back(population.at(i));
            }
            auto selected = tournament_selection(population, population.size() - elit_size);
            std::vector<subgraph_with_score> new_population;
            crossover_type == one_point ?
                    new_population = perform_one_point_crossover(selected, selected.size()) :
                    new_population = perform_uniform_crossover(selected, selected.size());
            mutation_type == bit_flip ?
                    perform_bit_flip_mutation(new_population) :
                    perform_bit_swap_mutation(new_population);
            std::vector<subgraph_with_score> new_population_with_elite;
            new_population_with_elite.reserve(population.size());
            new_population_with_elite.insert(new_population_with_elite.end(), new_population.begin(), new_population.end() );
            new_population_with_elite.insert(new_population_with_elite.end(), elite.begin(), elite.end());
            population = new_population_with_elite;

            if (!terminate_on_generations_number) {
                evaluate_population(population, fitness);
                if (get_best().score > best_score) {
                    best_score = get_best().score;
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
        evaluate_population(population, fitness);
        conv_log << get_best().score << std::endl;

        return get_best().subgraph;
    }
}

int main(int argc, char **argv) {
    using namespace mhe;
    adjacency_matrix_t example_graph = {
            1,1,1,1,
            1,0,1,
            1,1,
            0,
    };

    generate_graphviz_output(example_graph);

    std::vector<std::string> args(argv, argv+argc);
    std::string  selected_solver = "solve_random";
    int iterations = 10000;
    adjacency_matrix_t problem = example_graph;

    for (size_t i = 1; i < args.size(); ++i) {
        if (args[i] == "-s" || args[i] == "--solver") {
            if (i + 1 < args.size()) {
                selected_solver = args[++i];
            } else {
                return 1;
            }
        } else if (args[i] == "-i" || args[i] == "--iterations") {
            if (i + 1 < args.size()) {
                iterations = stoi(args[++i]);
            } else {
                return 1;
            }
        } else if (args[i] == "-p" || args[i] == "--problem") {
            if (i + 1 < args.size()) {
                problem = string_to_adjacency_matrix(args[++i]);
            } else {
                return 1;
            }
        } else {
            return 1;
        }
    }

    std::map<std::string, std::function<subgraph_t(adjacency_matrix_t, int)>> solvers;
    solvers["solve_hill_climbing"] = [&](auto problem, int iterations){return solve_hill_climbing(problem, iterations);};
    solvers["solve_tabu"] = [&](auto problem, int iterations){return solve_tabu_set(problem, iterations);};
    solvers["solve_tabu_list"] = [&](auto problem, int iterations){return solve_tabu_list(problem, iterations);};
    solvers["solve_tabu_avoid_snake"] = [&](auto problem, int iterations){return solve_tabu_avoid_snake(problem, iterations);};
    solvers["solve_sim_annealing"] = [&](auto problem, int iterations){return solve_sim_annealing(problem, iterations, [](int i){return 1000 * std::pow(0.99, (double)i);});};
    solvers["solve_random"] = [&](auto problem, int iterations){return solve_random(problem, iterations);};
    solvers["solve_random_n"] = [&](auto problem, int iterations){return solve_random(problem, iterations, 0.1);};
    solvers["solve_genetic_algorithm_iterations"] = [&](auto problem, int iterations){return solve_genetic_algorithm(problem, crossover_type::one_point, mutation_type::bit_flip, true, iterations);};
    solvers["solve_genetic_algorithm"] = [&](auto problem, int iterations){return solve_genetic_algorithm(problem, crossover_type::one_point, mutation_type::bit_flip);};

    generate_graphviz_output(problem);
    auto goal = goal_factory(problem);
    auto start_time = std::chrono::system_clock::now();
    auto result = solvers[selected_solver](problem, iterations);
    auto end_time = std::chrono::system_clock::now();
    auto computation_time =  std::chrono::nanoseconds(end_time - start_time);

    logger << "solution:" << std::endl;
    generate_graphviz_output(problem, result);
    int score = goal(result);
    std::cout << selected_solver << " " << score << " " << computation_time.count() << " " << result << std::endl;

    return 0;
}

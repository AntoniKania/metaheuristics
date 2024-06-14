for METHOD in solve_hill_climbing solve_tabu solve_tabu_list solve_tabu_avoid_snake solve_sim_annealing solve_random solve_random_n solve_genetic_algorithm_iterations solve_genetic_algorithm; do
  cmake-build-debug/metaheuristics --solver $METHOD --problem 011101111110101011000110101011100111011110010010100101000101110101 -i 100000 >> results.txt
done
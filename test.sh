for METHOD in solve_hill_climbing solve_tabu solve_tabu_list solve_tabu_avoid_snake solve_sim_annealing solve_random solve_random_n; do
  cmake-build-debug/metaheuristics $METHOD >> results.txt
done
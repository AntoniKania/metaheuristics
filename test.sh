for METHOD in solve_hill_climbing solve_tabu solve_tabu_list solve_tabu_avoid_snake solve_sim_annealing solve_random solve_random_n solve_genetic_algorithm_iterations solve_genetic_algorithm; do
  cmake-build-debug/metaheuristics --solver $METHOD --problem 101101011010110110011111011010100010110101111001000010011011101000 -i 10000 >> results.txt
done
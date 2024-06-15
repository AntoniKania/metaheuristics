# Metaheuristics
This repository contains code from the Metaheuristics course at my university. The codebase includes various metaheuristic algorithms designed to solve the clique problem. Implemented algorithms:
* Comprehensive review
* Random solution
* Hill Climbing
* Tabu Search
* Simulated Annealing
* Genetic Algorithm

### Clique Problem
The clique problem involves identifying the largest clique within a given graph. A clique is a subset of vertices where every pair is connected by an edge, forming a complete subgraph. The challenge is to determine the size of the largest clique, known as the graph's maximum clique. This problem is significant in various fields such as social network analysis, bioinformatics, and computer science, particularly due to its NP-complete complexity, making it computationally challenging to solve for large graphs.

### Graph Representation

        [ 0 1 1 1 1 ]
        [ 1 0 1 0 1 ]
    A = [ 1 1 0 1 1 ]
        [ 1 0 1 0 0 ]
        [ 1 1 1 0 0 ]

graph {
1 -- 2;
1 -- 3;
1 -- 4;
1 -- 5;
2 -- 3;
2 -- 5;
3 -- 4;
3 -- 5;
}

In this example, it can be seen that the matrix is mirrored across the diagonal, as the value at [4,2] is the same as at [2,4], representing the same connection, assuming that graphs are undirected.

        [ 0 1 1 1 1 ]     [ 0 1 1 1 1 ]
        [ 1 0 1 * 1 ]     [ 1 0 1 0 1 ]
    A = [ 1 1 0 1 1 ]  =  [ 1 1 0 1 1 ]
        [ 1 0 1 0 0 ]     [ 1 * 1 0 0 ]
        [ 1 1 1 0 0 ]     [ 1 1 1 0 0 ]

By dropping the lower triangle, we do not lose any information.

        [ 0 1 1 1 1 ]
        [   0 1 0 1 ]
    A = [     0 1 1 ]
        [       0 0 ]
        [         0 ]

This can be further simplified to the following matrix, as we are not interested in self-connections.

             (2)(3)(4)(5)
        (1) [ 1  1  1  1 ]
    A = (2) [    1  0  1 ]
        (3) [       1  1 ]
        (4) [          0 ]

### Goal Function
*This section explains why the modifications were made to the matrix in the first place.*

The adjacency matrix in its current state represents every possible edge in the graph, so if the graph is complete, it will only contain ones in the matrix. This simplifies the score functions as the

Let's view it in the previous example. The following subgraph will be checked and scored `{1, 1, 1, 0, 1}`, so nodes 1, 2, 3, and 5 are considered. The following adjacency matrix represents this subgraph:

             (2)(3)(5)
        (1) [ 1  1  1 ]
    A = (2) [    1  1 ]
        (3) [       1 ]

graph {
1 -- 2;
1 -- 3;
1 -- 5;
2 -- 3;
2 -- 5;
3 -- 5;
}

Using the following code from the `goal_factory` function, the adjacency matrix for this subgraph will be scored as `6`:

```c++
for (char i : subgraph_adjacency_matrix) {
    if (i == 1) {
        score++;
    } else {
        score -= 100;
    }
}
```c++
for (char i : subgraph_adjacency_matrix) {
    if (i == 1) {
        score++;
    } else {
        score -= 100;
    }
}
```
Here, it can also be seen that each missing connection in the subgraph decreases the score by 100 points.
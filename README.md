# PageRank
PageRank serial and parallel implementation

This code is written in C. It implements the pageRank algorithm and solves the linear problem of Ax=B using 
serial Power method, serial Gauss-Seidel method and parallel Gauss-Seidel method.

The tool used for the parallelization is openMP.

For the part of pagerank, a pseudo-network is constructed that resembles to the real internet problem.

Only the linear system is solved in parallel.

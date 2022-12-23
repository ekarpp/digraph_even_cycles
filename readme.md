# Introduction
This experimental software is done as part of my master's thesis [[1]](http://urn.fi/URN:NBN:fi:aalto-202212187112). It is an implementation of the randomized polynomial time algorithm presented in Section 4.2 of the publication by Björklund, Husfeldt and Kaski [[2]](https://doi.org/10.1145/3519935.3520030). The algorithm in question can determine the length of a shortest even cycle in a directed graph of $n$ vertices in time $O(n^{3+m})$ with high probability. Here $m$ is the square matrix multiplication exponent.

The algorithm uses a randomized algorithm design technique called *algebraic fingerprinting*. For a given graph we compute a *fingerprint* from which the length of a shortest even cycle can be determined with high probability. A fingerprint is a partial evaluation of the multivariate polynomial $per(A) - det(A)$, where $y$ is an indeterminate and $A$ is the adjacency matrix of the given digraph with loops appended to all vertices and each arc weight treated as an independent indeterminate. The correspondence between the fingerprints and shortest even cycles is due to Vazirani and Yannakakis [[3]](https://doi.org/10.1016/0166-218X(89)90053-X).

To compute the permanent of a matrix in polynomial time, the algorithm utilizes a specifically engineered characteristic 4 extension of a characteristic 2 finite field. The characteristic 4 extension, specifically a Galois ring, admits an $O(n^{2+m})$ time algorithm for the matrix permanent by an extension of techniques described by Valiant [[4]](https://doi.org/10.1016/0304-3975(79)90044-6).

# Running
Clone repository with `nauty` submodule included
```
git clone --recurse-submodules https://github.com/ekarpp/digraph_shortest_even_cycle && cd digraph_shortest_even_cycle
```
Run `make` to build everything or alternatively `make help` for other options. Requires `g++` and x86-64 microarchitecture with support for `PCLMULQDQ`, `BMI2`, and `AVX2` instruction set extensions.

Only mandatory argument to the algorithm binary, namely `digraph`, is `-f <path to graph file>`. It requires a file where line $i$ in the file (starting from zero) lists zero or more numbers separated with a space. Each number corresponds to an endpoint of an arc starting from vertex $i$. Some example files and graph generators can be found in `graphs` folder.

In addition to the implementation of the algorithm, this repository provides:
- Unit testing software, `digraph-tests` binary
- Nauty as a submodule to generate all digraphs with $n$ vertices (up to isomorphism)
- Various graph generator scripts
- Graphs and SLURM scripts used to perform experiments
- Benchmarking software for our implementations finite fields and the characteristic 4 extension
- Software to benchmark memory bandwidth

###### Bibliography
<sup><sub>
[1] Karppinen, E. Engineering an algorithm for the shortest even cycle problem. Master's thesis, School of Science at Aalto University, 2022.
</sub></sup>

<sup><sub>
[2] Björklund, A., Husfeldt, T., and Kaski, P. The shortest even cycle problem is tractable. In *STOC ’22: 54th Annual ACM SIGACT Symposium on Theory of Computing, Rome, Italy, June 20 - 24, 2022* (2022), ACM, pages 117–130.
</sub></sup>

<sup><sub>
[3] Vazirani, V. V., and Yannakakis, M. Pfaffian orientations, 0/1 permanents, and even cycles in directed graphs. In *Automata, Languages and Programming (Tampere, 1988)*, volume 317 of *Lecture Notes in Comput. Sci.* Springer, Berlin, 1988, pages 667–681.
</sub></sup>

<sup><sub>
[4] Valiant, L. G. The complexity of computing the permanent. *Theoret. Comput. Sci.* 8, 2 (1979), 189–201.
</sub></sup>

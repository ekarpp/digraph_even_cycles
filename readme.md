# Introduction
This experimental software is done as part of my master's thesis. It is an implementation of the randomized polynomial time algorithm presented in Section 4.2 of the publication by Björklund, Husfeldt and Kaski [[1]](https://doi.org/10.1145/3519935.3520030). The algorithm in question can determine the length of a shortest even cycle in a directed graph of $n$ vertices in time $O(n^{3+m})$ with high probability. Here $m$ is the matrix multiplication exponent.

The algorithm uses a randomized algorithm design technique called *algebraic fingerprinting*. For a given graph we compute a *fingerprint* from which the length of a shortest even cycle can be determined with high probability. A fingerprint is a partial evaluation of the multivariate polynomial $per(A + Iy) - det(A + Iy)$, where $y$ is an indeterminate and $A$ is the adjacency matrix of the given loopless digraph with each arc weight treated as an independent indeterminate. The correspondence between the fingerprints and shortest even cycles is due to Vazirani and Yannakakis [[2]](https://doi.org/10.1016/0166-218X(89)90053-X).

To compute the permanent of a matrix in polynomial time, the algorithm utilizes a specifically engineered characteristic 4 extension of a characteristic 2 finite field. The characteristic 4 extension admits a $O(n^{2+m})$ time algorithm for the matrix permanent by an extension of techniques described by Valiant [[3]](https://doi.org/10.1016/0304-3975(79)90044-6).

# Compilation
Executing `make` will build all binaries, while `make ${binary_name}` builds a specific binary. Requires `g++` and x86-64 microarchitecture with support for `PCLMULQDQ` and `BMI2` instruction set extensions. Additionally, the binaries utilizing $GF(2^{16})$ require support for `AVX` and `AVX2` instruction set extensions.

The number after the name of the binary corresponds to the exponent of the underlying finite field. Parallel binaries have the `PAR` suffix and non-vectorized $GF(2^{16})$ binaries have the `NOVEC` suffix. The binaries can be ran without arguments for further instructions.

In addition to the implementation of the algorithm, this repository provides:
- Unit testing software
- Various graph generator scripts
- Graphs and SLURM scripts used to perform experiments
- Benchmarking software for our implementations finite fields and the characteristic 4 extension
- Software to benchmark memory bandwidth

###### Bibliography
<sup><sub>
[1] Bjorklund, A., Husfeldt, T., and Kaski, P. The shortest even
cycle problem is tractable. In *STOC ’22: 54th Annual ACM SIGACT
Symposium on Theory of Computing, Rome, Italy, June 20 - 24, 2022*
(2022), ACM, pages 117–130.
</sub></sup>

<sup><sub>
[2] Vazirani, V. V., and Yannakakis, M. Pfaffian orientations, 0/1
permanents, and even cycles in directed graphs. In *Automata, Languages
and Programming (Tampere, 1988)*, volume 317 of *Lecture Notes in
Comput. Sci.* Springer, Berlin, 1988, pages 667–681.
</sub></sup>

<sup><sub>
[3] Valiant, L. G. The complexity of computing the permanent. *Theoret.
Comput. Sci.* 8, 2 (1979), 189–201.
</sub></sup>

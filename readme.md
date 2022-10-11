# Introduction
This experimental software is done as part of my Master's thesis. It implements the randomized polynomial-time algorithm presented in chapter 4.2 of a publication by Björklund, Husfeldt and Kaski [1]. The algorithm in question can determine the length of a shortest even cycle in a directed graph of $n$ vertices in time $O(n^{3+m})$ where $m$ is the matrix multiplication exponent.

The algorithm uses a technique called *algebraic fingerprinting*. For a given graph we compute a *fingerprint* from which the length of a shortest even cycle can be determined with high probability. A fingerprint is a partial evaluation of the multivariate polynomial per$(A + Iy) - $ det$(A + Iy)$, where $A$ is the adjacency matrix of the given graph with each arc weight treated as an independent indeterminate and $y$ an indeterminate. The correspondence between the fingerprint and shortest even cycles is due to Vazirani and Yannakakis [2].

To compute the permanent of a matrix in polynomial time, the algorithm utilizes a specifically engineered characteristic 4 extension of a characteristic 2 finite field. The characteristic 4 extension admits a $O(n^{2+m})$ time algorithm for the matrix permanent by an extension to techniques described by Valiant [3].

# Compilation
Executing `make` will build all binaries. Requires `g++` and x86-64 microarchitecture with support for `PCLMULQDQ` and `BMI2` instruction set extensions. Additionally, the binaries utilizing GF($2^{16}$) require support for `AVX` and `AVX2` instruction set extensions. The binaries can be ran for further instructions.

# Bibliography

[1] Bjorklund, A., Husfeldt, T., and Kaski, P. The shortest even
cycle problem is tractable. In *STOC ’22: 54th Annual ACM SIGACT
Symposium on Theory of Computing, Rome, Italy, June 20 - 24, 2022*
(2022), ACM, pages 117–130.

[2] Vazirani, V. V., and Yannakakis, M. Pfaffian orientations, 0{1
permanents, and even cycles in directed graphs. In *Automata, Languages
and Programming (Tampere, 1988)*, volume 317 of *Lecture Notes in
Comput. Sci.* Springer, Berlin, 1988, pages 667–681.

[3] Valiant, L. G. The complexity of computing the permanent. *Theoret.
Comput. Sci.* 8, 2 (1979), 189–201.

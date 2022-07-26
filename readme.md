# Introduction
This experimental software is done as implements the algorithm presented in chapter 4.2 of [The Shortest Even Cycle Problem is Tractable](https://arxiv.org/abs/2111.02992). The algorithm in question can determine the length of a shortest even cycle in a directed graph of $n$ vertices in time $O(n^{3+m})$ where $m$ is the matrix multiplication exponent.

The algorithm is based on algebraic fingerprinting. For a directed graph we compute its "fingerprint", namely the *parity cycle cover enumerator*. The fingerprint is engineered such that it can be used to determine the length of a shortest even cycle. Efficient computation of the fingerprint is the core of the algorithm. The computation is done by extending a finite field of characteristic two, to a finite ring of characteristic four and computing the determinant and permanent of the directed graph's adjacency matrix in the extension ring.

# Compilation
Executing `make` will build all binaries. Requires `g++` and x86-64 microarchitecture with support for `PCLMULQDQ` and `BMI2` instruction set extensions. Additionally, the binaries utilizing GF($2^{16}$) require support for `AVX` and `AVX2` instruction set extensions. The binaries can be ran for further instructions.

Dynamic Influence Analysis
========================

This is a real-time fully-dynamic index data structure designed for influence analysis on evolving networks.

## Usage
To use our index algorithm, please include "dim.hpp".

### Index construction
Our index can be build as follows:

    #include "dim.hpp"
    DIM dim;
    dim.init();
    dim.set_beta(0);
    for all vertex v
    	dim.insert(v); // Insert vertex v
    for all edge (u, v)
    	dim.insert(u, v); // Insert edge (u,v)
    dim.set_beta(beta); // Set \beta=32 (See paper)
	dim._adjust();

### Dynamic updates
* insert(int v): Insert vertex v
* erase(int v): Delete vertex v (and its incident edges)
* insert(int u, int v, double p): Insert edge (u,v) with probability p
* erase(int u, int v): Delete edge (u,v)
* change(int u, int v, double p): Change probability of edge (u,v) to p

### Influence queries
* infmax(int k): Extract an influential set of k vertices
* infest_naive(vector&lt;int&gt; &S): Estimate the influence of S (slow implementation)
* infest(vector&lt;int&gt; &S): Estimate the influence of S (fast implementation)
* infest(int v): Estimate the influence of v

### Example
    #include "dim.hpp"
    DIM dim;
    dim.init();
    dim.set_beta(0);

    dim.insert(0);
    dim.insert(1);
    dim.insert(0, 1, 0.4);
	dim.set_beta(32);
	dim._adjust();

    dim.insert(2);
    dim.insert(3);

    dim.insert(1, 2, 0.5);
    dim.insert(2, 0, 0.6);
    dim.insert(2, 3, 0.7);

    dim.infest(0); 
    dim.infmax(1);

## References

* Naoto Ohsaka, Takuya Akiba, Yuichi Yoshida, and Ken-ichi Kawarabayashi. **[Dynamic Influence Analysis in Evolving Networks](http://www.vldb.org/pvldb/vol9/p1077-ohsaka.pdf)**.
Proceedings of the VLDB Endowment, 9(12), pages 1077–-1088, 2016.

# scaled-minhash-large-metagenome-expt

Suppose we have a small genome G, and a large metagenome M. In this experiment, we will take C% of G, put it together with M, and call it M'.

Naturally, containment of G with in M', C(G, M') should ideally be C.

In this experiment, we will then determine C(G, M') using Scaled MinHash, and classic MinHash. For classic MinHash, we will be using Mash.
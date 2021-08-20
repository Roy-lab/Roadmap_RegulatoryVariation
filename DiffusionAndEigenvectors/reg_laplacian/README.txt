chasman@cs.wisc.edu, September 2014
-----------------------------------

Basic code related to constructing a regularized Laplacian kernel from a 
sif-format network (assumed UNDIRECTED). 

Requires numpy and SciPy's linear algebra goodness.

The regularized Laplacian kernel normalizes diffusion by node degree
and provides a parameter, lambda, which controls the smoothness.
Higher lambda results in greater smoothness. Lambda=1 is a fine place to start.

I've observed some weird effects when I've applied this to human-scale networks
with extremely exponential degree distributions, but that's another story for another time.

README contents:
[A] Main script: get_kernel_scores.py
[B] Helpers: kernel.py and graphy.py
[C] Input file format: graph in sif format
[D] Input file format: Node scores
[E] Examples (run with run_examples.sh; choose from star, chain, bg, degree)
[E1] Simple chain network
[E2] Simple star network 
[E3] General test network
[E4] Dramatic degree distribution network


[A] Main script: get_kernel_scores.py
---------------------------------
Two usage options:

(1) Constructs new kernel from graph.sif and lambda value (> 0), 
saves the kernel as graph_lambda.npy, and writes smoothed scores to stdout:

python get_kernel_scores.py input_scores.tab graph.sif lambda > output_scores.tab

(2) Reads kernel from graph_lambda.npy, writes smoothed scores to stdout. 
Need to re-supply lambda and graph as a basic compatibility check and in order 
to match node names properly.

python get_kernel_scores.py input_scores.tab graph.sif lambda graph_lambda.npy > output_scores.tab


[B] Helpers: kernel.py and graphy.py
---------------------------------
Contain various methods for reading files and doing kernel-related algebra.

[C] Input file format: graph in sif format
------------------------------------------
We read the graph in space-delimited sif format, which is easily read by Cytoscape. 
It looks like this:
NODE1 interaction NODE2

This code does nothing with the interaction type.
In graphy's readInteractions method, you can provide a different set of columns 
in order to read from different file formats.

By default, all node names are converted to upper case.

[D] Input file format: Node scores
-------------------------------
We read the scores from a simple tab-delimited file format.

NODE1	1.0
NODE2	2.0
# Comment line

By default, all node names are converted to upper case.

Output file format: Tab-delimited with header and comments
----------------------------------------------------------
The output file produced by get_kernel_scores.py looks like this:

# Lots of comments about progress and input data filenames
# more info about progress
name	input_score	output_score
A	0.5	0.381387
B	1.0	0.635834

[E] Examples (run_examples.sh; choose from chain, star, bg, degree)
-----------------------------------------
I've provided a couple of examples that may help you see
how the diffusion operates. 

[E1] Simple chain network
----------------------------------------
Input graph: tiny_chain.sif; input scores: tiny_chain_scores.tab
The graph is a chain of nodes A-B-C-D-E. A has a score of 1.0. 

The script run_chain.sh varies lambda; observe the change in smoothness in the
resulting scores.

[E2] Simple star network
---------------------------------------
Input graph: tiny_star.sif; input scores: tiny_star_scores.tab
The graph consists of two stars: A-{B,C,D,E,F} and G-{H,I}.
A and G have the same score (1.0). Observe that A's neighbors receive
lower scores than G, because the kernel normalizes diffusion by node degree.

[E3] General test network
---------------------------------------
Input graph: tiny_bg.sif; input scores: tiny_bg_scores.tab
Network comes from my thesis in the related work chapter.
Arbitrary scores chosen for testing this script.

[E4] Dramatic degree distribution network
-----------------------------------------
Input graph: tiny_degree_test.sif; input scores: tiny_degree_scores.tab
I just made this one to see what happens when I hide high-degree nodes.
Possibly a more detailed example to come; possibly not.








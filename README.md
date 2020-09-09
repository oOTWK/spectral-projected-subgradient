# Spectral Projected Subgradient
Implementation of subgradient methods based on research papers and computational results tested on benchmark set-covering problems.

## Table of contents
* [How to run](#how-to-run)
* [References](#references)
* [Introduction](#introduction)
* [Pseudo-code](#pseudo-code)
  * [Subgradient vector](#subgradient-vector)
  * [Initial Lagrangian multiplier vector](#initial-lagrangian-multiplier-vector)
  * [Spectral Projected Subgradient](#spectral-projected-subgradient-pseudo-code)
  * [Basic Subgradient Method](#basic-subgradient-method-pseudo-code)
* [Computational results](#computational-results)

## How to run
1. `cd spectral-projected-subgradient`
1. `make`
1. `./build/bin/subgradient file_path` to run spectral projected subgradient
1. `./build/bin/subgradient file_path -b upperbound` to run basic subgradient
1. `make clean`

## References
- Spectral Projected Subgradient (SPS)  
[1] Crema, A., Loreto, M., & Raydan, M. (2007) Spectral projected subgradient with a momentum term for the Lagrangian dual approach. *Computers and Operations Research*, 34(10), 3174–3186.
- Basic Subgradient Method (BSM)  
[2] Beasley, J.E. (1990) A Lagrangian heuristic for set-covering problems. *Naval Research Logistics*, 37(1), 151–164.
- Test problem sets  
scpnre, scpnrf, scpnrg, scpnrh from [Beasley’s OR library](http://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html)

## Introduction
The [set-covering problem](https://en.wikipedia.org/wiki/Set_cover_problem) (SCP) can be described as below.
<p align="center"><img src="../image/latex/original-scp.png" /></p>

One popular approach to solve SCP is branch-and-bound technique.  
**Branch**: The original problem can be reduced to two subproblems by fixing one variable at 0 or 1. This can be applied recursively and forms a binary tree.  
**Bound**: Since there are 2^n possible nodes in the binary tree, it is almost impossible to evaluate all nodes. We can prune or backtrack at a node if it is guaranteed that no better feasible primal solution can be found in its subtree. We can check this property by comparing lowerbound of the local optimal value of the subproblem with the current best primal solution which is a valid upperbound of the global optimal value.  

**Lowerbound**: The optimal value of linear programming (LP) relaxation of the original problem is a valid lowerbound of the optimal value of the original problem. LP relaxation problem is the same as the original problem except variable can be a real number (i.e., x \in [0, 1])  
Often, it is time-consuming to compute the optimal value of LP relaxation problem. Instead, we use dual of LP relaxation problem.  
<p align="center"><img src="../image/latex/dual-lp-relaxation.png" /></p>

Because the dual problem is maximization problem, any feasible solution of the dual problem is less than or equal to the optimal value of the LP relaxation problem. Hence, it is a valid lowerbound.

**Lagragian dual problem**: In this project, we will use a subgradient method to compute the lowerbound. A subgradient method is an iterative algorithm to optimize a non-differential convex function. A better formulation for a subgradient method is the Lagragian dual problem.  

It is known that the optimal value of the dual of LP relaxation is equal to the optimal value of the Lagragian dual problem. Hence, any feasible solution of the Lagragian dual problem is also a valid lowerbound.  
<p align="center"><img src="../image/latex/Lagrangian-dual.png" /></p>

In this project, I implemented two subgradient methods proposed in the literature, [1] and [2], and compared their performances. The experiments were conducted on the benchmark problems which are available on [Beasley’s OR library](http://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html).

The full SCP solver is under development. SPS is one componenet of the full SCP solver.

## Pseudo-code
### Subgradient vector
<p align="center"><img src="../image/latex/subg-vector.png" /></p>

### Initial Lagrangian multiplier vector
<p align="center"><img src="../image/latex/init-vector.png" /></p>

### Spectral Projected Subgradient pseudo-code 
<p align="center"><img src="../image/latex/sps-pseudo.png" /></p>

### Basic Subgradient Method pseudo-code 
<p align="center"><img src="../image/latex/bsm-pseudo.png" /></p>

## Computational results
##### scpnre1
![scpnre1](../image/graph/scpnre1.png)
##### scpnre2
![scpnre1](../image/graph/scpnre2.png)
##### scpnre3
![scpnre1](../image/graph/scpnre3.png)
##### scpnre4
![scpnre1](../image/graph/scpnre4.png)
##### scpnre5
![scpnre1](../image/graph/scpnre5.png)
##### scpnrf1
![scpnre1](../image/graph/scpnrf1.png)
##### scpnrf2
![scpnre1](../image/graph/scpnrf2.png)
##### scpnrf3
![scpnre1](../image/graph/scpnrf3.png)
##### scpnrf4
![scpnre1](../image/graph/scpnrf4.png)
##### scpnrf5
![scpnre1](../image/graph/scpnrf5.png)
##### scpnrg1
![scpnre1](../image/graph/scpnrg1.png)
##### scpnrg2
![scpnre1](../image/graph/scpnrg2.png)
##### scpnrg3
![scpnre1](../image/graph/scpnrg3.png)
##### scpnrg4
![scpnre1](../image/graph/scpnrg4.png)
##### scpnrg5
![scpnre1](../image/graph/scpnrg5.png)
##### scpnrh1
![scpnre1](../image/graph/scpnrh1.png)
##### scpnrh2
![scpnre1](../image/graph/scpnrh2.png)
##### scpnrh3
![scpnre1](../image/graph/scpnrh3.png)
##### scpnrh4
![scpnre1](../image/graph/scpnrh4.png)
##### scpnrh5
![scpnre1](../image/graph/scpnrh5.png)


MAXITER=300 | SPS Zmax | SPS time | BSM Zmax | BSM time
:-: | :-: | :-: | :-: | :-:
scpnre1 | 21.175433 | 0.157 | 20.955923 | 0.051
scpnre2 | 21.970871 | 0.227 | 21.626564 | 0.056  
scpnre3 | 20.265204 | 0.221 | 20.128360 | 0.053  
scpnre4 | 21.212269 | 0.172 | 20.708464 | 0.052  
scpnre5 | 21.010122 | 0.232 | 20.818565 | 0.042  
scpnrf1 | 8.303401 | 0.700 | 8.095496 | 0.138  
scpnrf2 | 9.415759 | 0.713 | 9.236443 | 0.130  
scpnrf3 | 8.647725 | 0.547 | 8.681868 | 0.131  
scpnrf4 | 7.769377 | 0.706 | 7.887587 | 0.141  
scpnrf5 | 6.878050 | 0.800 | 7.189729 | 0.173  
scpnrg1 | 158.884052 | 0.044 | 158.734726 | 0.041  
scpnrg2 | 141.383817 | 0.042 | 141.135389 | 0.040  
scpnrg3 | 147.512255 | 0.061 | 147.295610 | 0.041  
scpnrg4 | 147.716041 | 0.040 | 147.513275 | 0.041  
scpnrg5 | 147.329375 | 0.042 | 147.069810 | 0.041  
scpnrh1 | 47.411162 | 0.171 | 46.908672 | 0.152  
scpnrh2 | 48.297180 | 0.326 | 47.515954 | 0.152  
scpnrh3 | 44.647453 | 0.409 | 44.068299 | 0.157 
scpnrh4 | 43.397665 | 0.202 | 42.454367 | 0.155  
scpnrh5 | 41.305773 | 0.188 | 41.188819 | 0.146 

# anytime pmnk-landscapes
## *Anytime analysis of algorithms for the pMNK-landscapes problem*

### Index

* [Problem Description](#problem-description)
* [References](#references)
* [Project Supervisors](#project-supervisors)

### Problem Description

Multi-objective combinatorial optimization problems are combinatorial problems
that have more than one objective function to be optimized, e.g. minimizing
both the time and cost required to travel between two cities. Due to the 
conflicting nature of these problems there is often no single solution
that is optimal on all objective functions. Instead, there is a set of optimal
solutions, namely the Pareto set, that offer a trade-off between the objective
functions [[1]](#1). Unfortunately, finding the Pareto set in an acceptable amount
of time is often infeasible. So instead, many algorithms focus instead on 
finding a good approximation in a short amount of time, and indeed many algorithms 
have been proposed for multi-objective combinatorial optimization
problems over the past few years. In particular, there have been a large 
number of evolutionary algorithm and local search heuristics proposed, which
can often quickly achieve good solutions even if they can not often achieve
the complete Pareto set.

However, the study of these algorithms often focuses on analyzing the
quality of the approximation obtained within a given time budget, or, in
reverse, the time needed to achieve a given target quality. Unfortunately,
the real time budgets or target qualities that are relevant for the use of such
algorithms in real-world applications is often unknown and can vary widely.
So instead, it can be more useful to analyze the algorithms from an anytime
perspective, i.e. to study the evolution of the approximation quality with
respect to the elapsed time to allow for a broader analysis. Moreover, it
is relevant to study how certain properties of the problem, e.g. it's fitness
landscape [[2]](#2), affect the anytime behavior of the algorithms, which can be
invaluable for the tasks of automatic algorithm configuration and selection.
The aim of this project is to implement heuristic algorithms for the
ρMNK-landscapes problem [[3]](#3), and to study their performance from an any-
time perspective. The ρMNK-landscapes problem is a multi-objective 
combinatorial problem that allows us to tune its fitness landscape to 
some extent by setting the parameters: ρ, M, N, and K [[2]](#4). As such,
it is an interesting problem to study the anytime behavior of the algorithms 
and how this behavior is affected by the landscape of the problem.

### References

<a name="1">[1]</a> 
M. Ehrgott, Multicriteria Optimization. Springer-Verlag, second ed., 2005.

<a name="2">[2]</a>
F. Daolio, A. Liefooghe, S. Verel, H. Aguirre, and K. Tanaka, “Problem
features versus algorithm performance on rugged multiobjective combi-
natorial fitness landscapes,” Evolutionary Computation, vol. 25, no. 4,
pp. 555–585, 2017.

<a name="3">[3]</a>
S. Verel, A. Liefooghe, L. Jourdan, and C. Dhaenens, “On the structure
of multiobjective combinatorial search space: MNK-landscapes with cor-
related objectives,”European Journal of Operational Research, vol. 227,
no. 2, pp. 331–342, 2013.

<a name="4">[4]</a>
M. Laumanns, L. Thiele, and E. Zitzler, “Running time analysis of evo-
lutionary algorithms on a simplified multiobjective knapsack problem,”
Natural Computing, vol. 3, no. 1, pp. 37–51, 2004.

<a name="5">[5]</a>
E. Zitzler and L. Thiele, “Multiobjective optimization using evolutionary
algorithms — A comparative case study,” inProceedings of the 5th In-
ternational Conference on Parallel Problem Solving from Nature (PPSN
1998), pp. 292–301

<a name="6">[6]</a>
L. Paquete, M. Chiarandini, and T. Stützle,Pareto Local Optimum Sets
in the Biobjective Traveling Salesman Problem: An Experimental Study,
pp. 177–199. Lecture Notes in Economics and Mathematical Systems,
Springer Berlin Heidelberg, 2004.

<a name="7">[7]</a>
Jérémie Dubois-Lacoste, Manuel López-Ibáñez, Thomas Stützle, Anytime Pareto local search,
European Journal of Operational Research, Volume 243, Issue 2, 2015 Pages 369-385, ISSN 0377-2217
DOI: https://doi.org/10.1016/j.ejor.2014.10.062.

<a name="8">[8]</a>
Eckart Zitler, Simon Künzli, "Indicator-Based Selection in Multiobjective Search"

<a name="9">[9]</a>
Arnaud Liefooghe, Bilel Derbel. A Correlation Analysis of Set Quality Indicator Values in
Multiobjective Optimization. Genetic and Evolutionary Computation Conference (GECCO 2016), Jul 2016,
Denver, United States. hal-01159961v2A

<a name="10">[10]</a>
A.E. Eiben, J.E. Smith, Introduction to Evolutionary Computing, ISSN 1619-7127, 
DOI: https://doi.org/10.1007/978-3-662-44874-8
Alexandre D. Jesus, July 23, 2020

### Project Supervisors
- Alexandre D. Jesus
- Luis Paquete

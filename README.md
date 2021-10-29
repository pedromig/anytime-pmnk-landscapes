# anytime pmnk-landscapes

## Dependencies
CLI11

## Compilation

cmake -B build-S . -DCMAKE_BUILD_TYPE=Release  
cmake --build build

## Usage 

### General 

```
Usage: anytime-pmnk-landscapes [OPTIONS] instance ALGORITHM

Positionals:
  instance .dat:FILE REQUIRED           
      = pmnk-landscapes instance file path
      (instances can be generated using the rmnkGenerator.R script).

Options:
  -m,--maxeval UINT:NONNEGATIVE REQUIRED Needs: instance
            = maximum number of evaluations to be performed (stopping criterion).
  -s,--seed UINT:NONNEGATIVE Needs: instance
            = pseudo random generator seed used by the search heuristics.
  -o,--output TEXT Needs: instance      
            = specifies the file to which the output stream should be redirected.
  -r,--hvref FLOAT ... Needs: instance  
            = reference point considered in the hypervolume calculation.
  -h,--help                             
            = print this help message and exit.
  -H,--help-all                         
            = expand all help.

Algorithms:
  GSEMO    Run the global simple evolutionary multiobjective optimizer algorithm
           on the instance.
  PLS      Run the pareto local search algorithm on the instance.
  IBEA     Run the indicator-based evolutionary algorithm on the instance.
```

<details>
<summary>Examples</summary>

* Common Usage 
```
./anytime-pmnk-landscapes
```

</details>

### GSEMO

```
Run the global simple evolutionary multiobjective optimizer algorithm
on the instance.
Usage: anytime-pmnk-landscapes [OPTIONS] instance GSEMO
```

<details>
<summary>Examples</summary>
  
</details>

### PLS
```
Run the pareto local search algorithm on the instance.
Usage: anytime-pmnk-landscapes [OPTIONS] instance PLS [OPTIONS]

Options:
  -h,--help       
      = print this help message and exit.
  -H,--help-all                         
      = expand all help.
  -a,--pls-acceptance-criterion 
    ENUM:value in {BOTH->2,DOMINATING->1,NON_DOMINATING->0} OR {2,1,0}
        = acceptance criterion considered whilst running pls.
            => (NON_DOMINATING): accept every non-dominated neighbor.
            => (DOMINATING): accept only neighbors that dominate current solution.
            => (BOTH): first try to accept only neighbors that dominate
                the current solution, if none exist accept non-dominated solutions.
  -e,--pls-neighborhood-exploration 
    ENUM:value in {BEST_IMPROVEMENT->0,BOTH->2,FIRST_IMPROVEMENT->1} OR {0,2,1}
        = neighborhood exploration criterion considered whilst running pls.
            => (BEST_IMPROVEMENT): explore every acceptable neighboor.
            => (FIRST_IMPROVEMENT): stop once on neighbor is accepted.
            => (BOTH): use FIRST_IMPROVEMENT until PLS stops, afterwards use BEST_IMPROVEMENT
```

### IBEA

```
Run the indicator-based evolutionary algorithm on the instance.
Usage: anytime-pmnk-landscapes IBEA [OPTIONS] INDICATOR MUTATION CROSSOVER SELECTION

Options:
  -h,--help   
    = print this help message and exit.
  -H,--help-all                         
    = expand all help.
  -p,--pop-size UINT:NONNEGATIVE REQUIRED
    = max population size.
  -g,--generations UINT:NONNEGATIVE REQUIRED
    = number of generations. (stopping criterion)
  -k,--scaling-factor FLOAT:NONNEGATIVE REQUIRED
    = scaling factor.
  -a,--adaptive                         
    = use the adaptive version of the algorithm

Indicators:
IHD
  Run using the hypervolume indicator

EPS
  Run using the epsilon (+) indicator

Mutation Operators:
UniformMutation
  Run using a uniform mutation operator
  Options:
    -p,--mutation-probability FLOAT:FLOAT in [0 - 1] REQUIRED
        = probability of occurrence of a mutation in the individual's genotype

Crossover Operators:
NPointCrossover
  Run using a n-point crossover operator
  Options:
    -p,--crossover_probability FLOAT:FLOAT in [0 - 1] REQUIRED
        = probability of occurrence of a mutation in the individual's genotype
    -n,--n-points UINT:NONNEGATIVE REQUIRED
        = number of randomly picked crossover points

UniformCrossover
  Run using a uniform crossover operator
  Options:
    -p,--crossover_probability FLOAT:FLOAT in [0 - 1] REQUIRED
        = probability of occurrence of a mutation in the individual's genotype


Selection Operators:
KWayTournament
  Run using a k-way tournament selection operator
  Options:
    -s,--matting-pool-size UINT:NONNEGATIVE REQUIRED
        = target size of the matting pool to be obtained from the selection step
    -t,--tournament-size UINT:NONNEGATIVE REQUIRED
        = size of the tournament used for individual selection
```

<details>
<summary>Examples</summary>

</details>


## API


For more information about the usage of this api please consult the available
[documentation](docs/documentation/documentation.pdf)

## Acknowledgments

This project was only made possible thanks to the support of a research grant 
at the Center for Informatics and Systems of the University of Coimbra (CISUC),
with reference number UIDB/00326/2020, awarded by the Foundation for Science 
and Technology (FCT).
 
## Contributors
- [Alexandre D. Jesus](https://adbjesus.com/)

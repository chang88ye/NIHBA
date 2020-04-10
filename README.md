# **NIHBA**

![](nihba.png)

NIHBA (Network Interdiction with Hybrid Benders Algorithm) is a software platform for identifying as many genetic modification strategies as possible for metabolic engineering. It is based on genome-scale metabolic models (GSMM), and considers the host 
cell and a metabolic engineer as two different agenets playing a game: the metabolic engineer attempts to maximally violate the host cell's primary objective. This game is called network interdiction in game theory. NIHBA employs a hybrid Bender's decomposition algorithm to solve this network interdiction problem, thereby achieving high efficiency without numeric issues in search of best design strategies.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

## Prerequisites

MATLAB2016 and later
cobratoolbox (latest version: 3.0)
Gurobi7.5 (or other software packages containing LP and ILP solvers)

## Installation

1, Download the whole NIHBA package directly from this repository or clone using
$ git clone --depth=1 https://github.com/Chang88ye/NIHBA.git NIHBA

2, Launch MATLAB and change the working directory to NIHBA as the current folder. Then run nihba_setup

>> nihba_setup

Note: nihba_setup will help you add all the NIHBA files to matlab directory. It will also test if you have all prerequisites available for NIHBA.

It is extremely important to retrieve unbounded rays for unbounded LP in NIHBA. This can be done by setting the InfUnbdInfo parameter to 1 for gurobi, and then return the unbounded ray in solveCobraLP (i.e., x=resultgurobi.unbdray in Cobratoolbox). The following shows what should look like in the unbounded case of solveCobraLP for gurobi.

        case 'UNBOUNDED'
                stat = 2; % unbounded
                x=resultgurobi.unbdray;

If you encounter that +optimality conditions are not satisfied in solveCobraQP+, it is likely caused by the numerical issue in the 'analyzeOptKnock'function in the cobratoolbox. In this case, you can replace it with the one provided in NIBA.

3, Run a toy example that maximises succinate production from the E. coli core model.

>> run_example

4, If you want to apply NIHBA to other genome-scale models for different products, just change target product in run_example

5, If you wish to visualise the results of NIHBA as done in our paper, please feel free to use or modify the code provided at 
https://github.com/chang88ye/NIHBA-supplement

## License

Released under the MIT license. All included network models provided under their respective licenses.

## Cite Us

Jiang S, Wang Y, Kaiser M, Krasnogor N.(2020) NIHBA: A network interdiction approach with Hbybrid Benders Algorithm for Metabolic Engineering. https://www.biorxiv.org/content/10.1101/752923v2.full

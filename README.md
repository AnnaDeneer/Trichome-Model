# Trichome-Model
MATLAB code to simulate trichome patterns as used in "Identification of the trichome patterning core network using data from weak ttg1 alleles to constrain the model space".

## Usage
The main file is simModel.m, from this script the model equations are solved and the pattern analyzed. To run the simulation, a parameter set is needed. Two example parameter sets are supplied in Parameterset_wt.mat and Parameterset_ttg19.mat, representing a wild-type situation and the corresponding ttg1-9 mutation, respectively. Example:<br/>
k = load('Parameterset_wt.mat');<br/>
[t,y] = simModel(k, 1);

module load gurobi/900
module load matlab/R2017a
module load perl/5.10.1
module load git/2.9.5
git config --global http.sslverify "false"
matlab < myFlux.m > run.log

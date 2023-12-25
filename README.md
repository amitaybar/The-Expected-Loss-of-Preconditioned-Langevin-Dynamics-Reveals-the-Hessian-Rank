# The-Expected-Loss-of-Preconditioned-Langevin-Dynamics-Reveals-the-Hessian-Rank
This is the code for the paper "The Expected Loss of Preconditioned Langevin Dynamics Reveals the Hessian Rank", AAAI 2024

The first file to be run is LinearNets.py. It creates the data for the linear NN experiment. The data can be saved in the working directory or moved to a 'results' directory.

The file LinearNetsLossPlot.m creates Figure 1 by using the file dim_2_MonteCarlo_1 that is created by LinearNets.py, and is also attached to this zip file for convenience.

In the U&S_method directory, there is the implementation of the U&S method according to "Fast methods for estimating the numerical rank of large matrices". In International Conference on Machine Learning, 468–477. PMLR, 2016, Ubaru, S.; and Saad, Y.
This is the implementation published by the authors. The changes we made in the code are the efficient computation of the matrix-vector product according to the structure of the Hessian of linear NN. 

The main file in the U&S_method directory is MainUnS_method.m. It creates and saves RankEstTotErrFastMAt.mat and RankEstTotFastMat.mat.

The file LinearNetsRankEst.m uses the data created by LinearNets.py and MainU&S_method.m to produce Figure 2 and Table 1. For convenience, the files RankEstTotErrFastMAt.mat and RankEstTotFastMat.mat are attached to this zip file (LinearNets.py still needs to be run).

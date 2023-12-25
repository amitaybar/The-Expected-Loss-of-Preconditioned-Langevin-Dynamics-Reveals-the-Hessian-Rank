
This is the published code by the authors of "Fast methods for estimating the numerical rank of large matrices". In International Conference on Machine Learning, 468–477. PMLR, 2016, Ubaru, S.; and Saad, Y.
The code in this folder implements the rank estimation method we refer to as U&S.

The changes we made in the code are the efficient computation of the matrix-vector product according to the structure of the Hessian of linear NN. More details are in our paper "The Expected Loss of Preconditioned Langevin Dynamics Reveals the Hessian Rank".

The main file is MainUnS_method.m. It saves RankEstTotErrFastMAt.mat and RankEstTotFastMat.mat which are used by LinearNetsRankEst.m to produce Figure 2 and the data in Table 1.


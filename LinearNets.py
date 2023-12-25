import numpy as np
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import time
import os
from math import sqrt
from HessianPapyan import Hessian
import scipy


os.environ["CUBLAS_WORKSPACE_CONFIG"]=":4096:8"

DBUG = 0
dev = 'cuda:0'
eta     = 1e-4
KIter   = int(15e3)
KIter4Mean = int(10e3)
KMonteCarlo = 100 # 50
seed = 30
sigma_sqr = 2e-5
KDepth = 5
DimVec = [8, 16, 32, 64, 128, 256]
FixedDataStatsFlag = 0
KDiffDim = len(DimVec)


total_params = KDepth*DimVec[0]**2
Factor =  sigma_sqr

np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.backends.cudnn.benchmark = True
torch.use_deterministic_algorithms(True)

t = time.time()


for dim in range(KDiffDim):


    KDimIn  = DimVec[dim]
    KDimOut = DimVec[dim]


    Sigma_x  = torch.eye(KDimIn, device=dev, dtype=torch.double)
    Sigma_yx_r = torch.randn(KDimOut,KDimIn, device=dev, dtype=torch.double)
    Sigma_yx = Sigma_yx_r.matmul(Sigma_yx_r.transpose(0,1))/sqrt(KDimIn)
    T = Sigma_yx.matmul(torch.linalg.inv(Sigma_x))
    Lambda, V = torch.linalg.eigh(T)

    T_KDepth_root = V.matmul(torch.diag(Lambda).pow(1/KDepth)).matmul(V.transpose(0,1))

    I = torch.eye(KDimOut,device=dev, dtype=torch.double)
    
    for mc in range(KMonteCarlo):
        WTensor = []
        Loss = np.zeros(KIter)
        file_name = 'dim_'+str(dim)+'_MonteCarlo_'+str(mc)
        for ii in range(KDepth):
            WTensor.append(nn.Linear(in_features=KDimIn,out_features=KDimIn,bias=False, dtype=torch.double))
            WTensor[ii].weight.data = T_KDepth_root
        
        PTmp = torch.randn(KDimIn,KDimIn, device=dev, dtype=torch.double)
        Q, R = torch.linalg.qr(PTmp)
        P = Q

        WTensor[0].weight.data = P.matmul(WTensor[0].weight.data)
        for ii in range(1,KDepth-1):
            WTensor[ii].weight.data = P.matmul(WTensor[ii].weight.data).matmul(torch.linalg.inv(P))
        WTensor[KDepth-1].weight.data = WTensor[KDepth-1].weight.data.matmul(torch.linalg.inv(P))

        net = nn.Sequential(*WTensor).to(device=dev)
        with torch.no_grad():
            OpMatrix = net(I)

        
        torch.linalg.norm(OpMatrix-T,'fro')

        
        net.requires_grad_ = True


         


        # Our approach
        for it in range(KIter):
            net.zero_grad()
            A = net(I)
            loss = torch.trace(A.matmul(Sigma_x).matmul(A.transpose(0,1)))-2*torch.trace(A.matmul(Sigma_yx.transpose(0,1)))
            loss.backward()
            for layer in WTensor:
                layer.weight.data.add_(layer.weight.grad, alpha=-eta)
                layer.weight.data.add_(torch.randn_like(layer.weight), alpha=sqrt(sigma_sqr)*sqrt(eta))
            Loss[it] = loss.detach().cpu().numpy()

        
        RankEst = (np.mean(Loss[-KIter4Mean:]) - Loss[0])/Factor*4
        EstError = (RankEst-KDimIn**2)/(KDimIn**2)*100
        print("Dim: {}, Error: {}%".format(KDimIn,EstError))
        mdict = {'T':T.cpu().numpy(),
                 'P':P.cpu().numpy(),
                 'Loss':Loss,
                 'RankEst':RankEst,
                 'DimIn':KDimIn,
                 'Depth':KDepth}
        

        scipy.io.savemat(file_name, mdict)

total_time = time.time()-t
        
print(total_time)

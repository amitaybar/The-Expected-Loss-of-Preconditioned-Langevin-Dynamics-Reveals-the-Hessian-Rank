%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code produces Figure 1 in the paper. It requires the file
% dim_2_MonteCarlo_1.mat, which is created by LinearNets.py
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc
%%
addpath(genpath([pwd]));
DataFoundFlag = 1;
%%
try
    load('dim_2_MonteCarlo_1','-mat'); % Created by LinearNets.py
catch
    disp('Missing Data files that should be created by LinearNets.py');
    disp('Please run LinearNets.py first and save ''dim_2_MonteCarlo_1s'' to the current folder');
    DataFoundFlag = 0;
end
if DataFoundFlag
    KDepth      = double(Depth);
    KDimIn      = double(DimIn);
    Sigma_x     = eye(KDimIn);
    sigma2      = 2e-5;
    Sigma       = sigma2*kron(eye(KDimIn^2),eye(KDepth));
    G           = kron(eye(KDimIn^2),eye(KDepth));
    Eta         = 1e-4;
    KIter       = 15e3;
    %%
    for ii=1:KDepth
        WTensor(:,:,ii) = T^(1/KDepth);
    end
    
    WTensor(:,:,1)      = P * WTensor(:,:,1);
    for ii=2:KDepth-1
        WTensor(:,:,ii) = P*WTensor(:,:,ii)*inv(P);
    end
    WTensor(:,:,KDepth) = WTensor(:,:,KDepth) * inv(P) ;
    
    %% Hessian computation
    PhiMat              = [];
    for k=1:KDepth
        Phi_kTmp1       = MatrixMulFuc(WTensor(:,:,1:k-1))*Sigma_x^0.5;
        Phi_kTmp2       = MatrixMulFuc(WTensor(:,:,k+1:KDepth));
        Phi_k           = kron(Phi_kTmp1, Phi_kTmp2');
        PhiMat          = [PhiMat Phi_k'];
    end
    PhiMat              = PhiMat';
    HessianTheory       = 2*(PhiMat*PhiMat');
    
    %%
    [P1,Lambda]         = eig(HessianTheory);
    J                   = diag((diag(Lambda)>1e-10)); %for numeric issues
    %%
    I                   = eye(size(G));
    tVec                =(0:KIter-1) * Eta;
    TheoreticLossTime   = 0.25*sigma2*(KDimIn^2*KDepth-sum(exp(-2*diag(Lambda)*tVec),1));
    %%
    figure; hold on
    plot((1:KIter),Loss-Loss(1));
    plot((1:KIter),TheoreticLossTime,'r')
    yline(trace(Sigma*G*(P1*J/P1))*0.25,'Linewidth',2);
    xlim([-0.05*KIter KIter]) ;
    xlabel('Iteration');
    ylabel('Loss');
    legend('Empirical Loss','Theory (Thm 1)','Theory (Prop 1)','interpreter','latex');
    box on; grid on
    width = 3.3;
    height = 3/4*width;
    FontSize = 8;
    LineWidth = 1.5;
    MarkerSize = 2;
    LegendFontSize = 6;
    Res = '-r1200';
    FileName = 'LinNetLoss';
    myPrint(FileName,width,height,FontSize,LineWidth,MarkerSize,LegendFontSize,Res)
end
%%

function MatMul = MatrixMulFuc(MatTensor)
if size(MatTensor,3)>0
    MatMul = MatTensor(:, :, 1);
else
    MatMul = eye(size(MatTensor,1));
end
for k = 2:size(MatTensor,3)
    MatMul = MatTensor(:, :, k) *  MatMul;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code produces Figure 2 in the paper. It requires data 
% which is created by LinearNets.py and the data created by U&S_method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
clc
%%
addpath(genpath(pwd));
%%
% FolderName  = fullfile(pwd,'results',filesep); 
if exist('results','dir')
    FolderName  = fullfile(pwd,'results',filesep); 
else
    FolderName  = fullfile(pwd,filesep); 
end
%% Loading the results for the U&S method
UnSDataFound = 1;
try
load('RankEstTotFastMat.mat')
load('RankEstTotErrFastMAt.mat')
catch
    disp('Please run MainU&S_method.m in U&S_method folder first, and copy the two mat files to the current folder.');
    UnSDataFound = 0;
end
%% Colors definition
BlueColDef  = [0 0.4470 0.7410];
RedColDef   = [0.8500    0.3250    0.0980];
%%
KMonteCarlo = 100;
KDims       = 6;
DimVec      = [8 16 32 64 128 256];
KIter4Mean  = 10e3;
sigma2      = 2e-5;
%%
DataFoundFlag = 1;
%%
for dim=0:KDims - 1
    for mc=0:KMonteCarlo-1    
        FileName = [FolderName 'dim_' num2str(dim) '_MonteCarlo_' num2str(mc)];
        try
            load(FileName,'-mat');
        catch
            disp('Missing Data files that are created by LinearNets.py');
            disp(['Please run LinearNets.py first and save the data to: ' FolderName]);
            DataFoundFlag = 0;
            break;
        end
        MeanLossTot(dim+1,mc+1) = mean(Loss(end-KIter4Mean :end))-Loss(1);
        RankEstTotErr(dim+1,mc+1) = RankEst - double(DimIn)^2;               
        RankEstTot(dim+1,mc+1) = RankEst;        
    end
    LossTheory(dim+1) = DimVec(dim+1)^2*sigma2/4;
end
%% 
if DataFoundFlag % This Matlab file requires data from LinearNets.py
    %% Values in Table 1 
    Ours_RMSE   = [sqrt(mean(RankEstTotErr.^2,2))./DimVec'.^2*100]'
    if UnSDataFound
        UnS_RMSE    = [sqrt(mean(RankEstTotErrFastMAt.^2,2))./DimVec'.^2*100]'
    end
    %%
    figure; hold on
    plot(10.^[1 5],10.^[1 5],'--','color',[0 0 0 0.5]);
    h1 = errorbar(DimVec.^2,mean(RankEstTot,2),std(RankEstTot,1,2),'color',BlueColDef);
    if UnSDataFound
        h2 = errorbar(DimVec.^2,mean(RankEstTotFastMat,2),std(RankEstTotFastMat,1,2),'color',RedColDef);
        legend([h1 h2],{'Ours','U&S'})
    end    
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    xlabel('Hessian rank');
    ylabel('Hessian rank estimation');
    xlim([1e1 1e5]);
    box on;
    width = 3.3;
    height = 3/4*width;
    FontSize = 8;
    LineWidth = 1.5;
    MarkerSize = 2;
    LegendFontSize = 6;
    Res = '-r1200';
    FileName = 'LinNetRankEst';
    myPrint(FileName,width,height,FontSize,LineWidth,MarkerSize,LegendFontSize,Res)
end

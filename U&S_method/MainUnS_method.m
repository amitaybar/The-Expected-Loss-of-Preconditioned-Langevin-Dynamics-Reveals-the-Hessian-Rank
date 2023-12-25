%% Computes the approximate rank of a sample matrix using Chebyshev filters
% This code was published by the authors of  "Fast methods for estimating the numerical rank of large matrices",
% In International Conference on Machine Learning, 468?477. PMLR, 2016, Ubaru, S.; and Saad, 
% We updatd the code to efficiently compute the matrix-vector product using
% the structure of the Hessian of linear NN

close all
clear
%%-------------------- NOTE: small gap requires very high degree
%%
rng(0,'twister');
%%    parameters
deg   = 50;
eps=0.01;
nvecs = 300;
AutoGap=1;
Npts=100;

%plot flags
bDo_plot = 1;
Filter_plot=1;
%%--------------------
%% load matrices and their singular values
%  load 'netz4504.mat'
% load 'netz4504_SVD.mat'
%   X = Problem.A;
%   n = size(X,2);
% d=S.s;
% ACompare = X'*X;
% A = HessianTheory;
% n = size(A,2);
% d = d .^ 2;
% d = eig(A);
%%
%%
addpath(genpath(fileparts(pwd)));
DataFoundFlag = 1;
%%
FolderName = fileparts(pwd);
if exist('results','dir')
    FolderName  = fullfile(FolderName,'results',filesep); 
else
    FolderName  = fullfile(FolderName,filesep); 
end
KMonteCarlo = 100;
KDims = 6;
DimVec = [8 16 32 64 128 256];
tic
%%
for dim=0:KDims - 1
    for mc=0:KMonteCarlo-1
        FileName = [FolderName 'dim_' num2str(dim) '_MonteCarlo_' num2str(mc)];
        try
        load(FileName,'-mat');
        catch
            disp('Missing Data files that should be created by LinearNets.py');
            disp(['Please run LinearNets.py first and save the data to: ' FileName]);    
            DataFoundFlag = 0;
            break
        end
        %%
%         sqrt(mean(RankEstTotErr.^2,2))./DimVec'.^2*100;
        
%         load('dim_0_MonteCarlo_11','-mat');
        KDepth = double(Depth);
        KDimIn = double(DimIn);
        Sigma_x  = eye(KDimIn);
        sigma2 = 2e-5;
        
        WTensor = zeros(KDimIn,KDimIn,KDepth);
        for ii=1:KDepth
            WTensor(:,:,ii) = T^(1/KDepth);
        end
        WTensor(:,:,1) = P * WTensor(:,:,1);
        for ii=2:KDepth-1
            WTensor(:,:,ii) = P*WTensor(:,:,ii)/P;
        end
        WTensor(:,:,KDepth) = WTensor(:,:,KDepth) / P ;
        Phi_kTmp1 = zeros(KDimIn,KDimIn,KDepth+1);
        Phi_kTmp2 = zeros(KDimIn,KDimIn,KDepth+1);
        for k=1:KDepth
            Phi_kTmp1(:,:,k)       = MatrixMulFuc(WTensor(:,:,1:k-1))*Sigma_x^0.5;
            Phi_kTmp2(:,:,k)       = MatrixMulFuc(WTensor(:,:,k+1:KDepth));
        end
        
        A = @(v)(HessianVecEffCalc(Phi_kTmp1,Phi_kTmp2,v,KDimIn,KDepth));
        % A = @(v)(ACompare*v);
        n = KDimIn^2*KDepth;
        % n = size(ACompare,2);
        d(2) = eigs(A,n,1);
        d(1) = 0;
        
        
        %% Find threshold and rank
        % fprintf('Finding threshold and rank using KPM of deg = %d\n',deg);
        
        % n = size(A,2);
        % d = sort(d);
        %d(1)=0;
        %%-------------------- lmin - lmax
        lmin = d(1)-0.00001;
        lmax = d(end)+0.00001;
        %%--------------------
        t1=cputime;
        [r,g,zz,z1,xx,yy]=KPM_DOS_and_Rank(A,lmin,lmax,deg,nvecs,Npts,AutoGap,n);
        RankEstTotErr(dim+1,mc+1) = r - double(DimIn)^2;
        RankEstTot(dim+1,mc+1) = r;
    end
end
TimeDur = toc/60
%eig(full(A),10000);
% t2=cputime;
% fprintf('Estimated threshold      = %7.5f \n',g) ;
% fprintf('Estimated rank      = %7.5f \n',r) ;
% num = sum(d > (g));
% fprintf('Exact count above threshold = %7.5f \n',num)
% t2-t1;

%% PLots
if DataFoundFlag
    RankEstTotFastMat = RankEstTot;
    RankEstTotErrFastMAt = RankEstTotErr;
    save('RankEstTotFastMat','RankEstTotFastMat');
    save('RankEstTotErrFastMAt','RankEstTotErrFastMAt');
end
% if (bDo_plot)
%     figure(1)
%     clf
%     hold on
%     plot(xx(1:end), yy(1:end),'b-o','linewidth',2);
%     %axis([0 4.5 0 1])
%     txt = sprintf('DOS with KPM, deg M = %3d', deg );
%     title(txt,'fontsize',20);
%     xlabel('\lambda','fontsize',20);
%     ylabel('\phi ( \lambda ) ','fontsize',20);
%     box on
%     hold off
% end
% if (Filter_plot)
%     figure(2)
%     %%-------------------- plotting -- plot line cumulative estimate
%     plot(1:nvecs,zz,'k-','linewidth',2)
%     hold on;
%     %%-------------------- plot each value (Pv,v) for each sample vect.
%     plot(1:nvecs,z1,'o','linewidth',2)
%     %%-------------------- plot exact count - straight line
%     plot([0 nvecs], [num num], 'r--', 'linewidth',2)
%     %%----------------------- plot standard deviation of z1
%     titt = sprintf('Chebyshev filter method');
%     xlab = sprintf('Number of vectors (1 -> %d)',nvecs);
%     ylab = sprintf('Estimed # eigenvalues in interval');
%     title(titt,'fontname','Helvetica','fontsize',20);
%     xlabel(xlab,'fontname','Helvetica','fontsize',20);
%     ylabel(ylab,'fontname','Helvetica','fontsize',20);
%     set(gca, 'FontSize', 20 );
%     hleg = legend('Cumulative Filter','(Pv,v)','Exact');
%     %%              'Filter + sigma','Filter - sigma');
%     set(hleg,'FontSize',22,'FontWeight','normal','location',...
%         'northeast','Color','none' );
% end

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

function ResHV = HessianVecEffCalc(Phi_kTmp1,Phi_kTmp2,v,KDimIn,KDepth)
vMat = reshape(v,KDimIn,KDimIn,KDepth);
PhiSum = 0;
for k=1:KDepth
    Phi_kVec = Phi_kTmp2(:,:,k)*vMat(:,:,k)*Phi_kTmp1(:,:,k);
    PhiSum = PhiSum + Phi_kVec;
end
ResHV = zeros(KDimIn^2*KDepth,1);
for k=1:KDepth
    Phi_kPhiSum = Phi_kTmp2(:,:,k)'*PhiSum*Phi_kTmp1(:,:,k)';
    Indx        = 1+(k-1)*KDimIn^2 : k*KDimIn^2;
    ResHV(Indx) = 2*Phi_kPhiSum(:);
end
end
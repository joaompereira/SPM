clear all
clc
addpath '../'
addpath '../helper_functions/'
%addpath(genpath('../other_packages/'))

% Test SPM against other methods for GPCA
%  This test considers the case where all subspaces have the same dimension


k = 3;    %dimension of subspaces
K = 6;    %dimension of ambient space
N = 300;  %number of points in each subspace
n = 2;    %number of groups

noisevals = 10.^linspace(-2,0,25);
nnoise = length(noisevals);

%prepare dataset
X=zeros(K,N,n);
s=zeros(N,n);
for in = 1:n
    basis = orth(randn(K,k));
    X(:,:,in) = basis*rand(k,N);
    s(:,in) = in;
end

X = reshape(X,K,[]);
s = s(:);

Algs = {'SPM', @(X,k,n) gpca_SPM(X,k,n);
        %'PDA', @(X,k,n) gpca_PDA(X,n,size(X,1)-k);
        %'GPCA-V', @(X,k,n) gpca_voting(X, k*ones(1,n) ,'postoptimization',true);
        };

nalgs = size(Algs,1);
time = zeros(nalgs, nnoise);
missrate = zeros(nalgs, nnoise);

ntries = 100;

for j = 1:nnoise
    j
    
    missratej = zeros(nalgs, ntries);
    timej = zeros(nalgs, ntries);
    
    for tries = 1:ntries
        
        sigma = noisevals(j);
        X_noisy = X + sigma*randn(size(X));

        for i = 1:nalgs
            
            timer = tic;

            groups = Algs{i,2}(X_noisy,k,n);
            
            timej(i, tries) = toc(timer);

            %compute misclassification rate
            groups = bestMap(s,groups);
            missratej(i,tries) = sum(s(:) ~= groups(:)) / length(s);
        end
        
    end
    
    missrate(:,j) = mean(missratej,2);
    time(:,j) = mean(timej,2);
end

% Save results
save results/test_gpca_vs.mat

% Plot results
hf = figure(1);
ha = semilogx(noisevals,missrate,'LineWidth',2);
ha(2).LineStyle = ':';
ha(3).LineStyle = '--';
legend(Algs(:,1),'Location','southeast');
xlabel('$\sigma$','Interpreter','latex');
ylabel('Misclassification error (7.4)',...
       'Interpreter','latex');
pdfprint('results/test_gpca_vs', hf);


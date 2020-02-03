clearvars
clc
addpath '../'
addpath '../helper_functions/'

k = 3;    %dimension of subspaces
K = 20;    %dimension of ambient space
%number of points in each group
Nval = round(10.^linspace(2,5,10));  
nvals = length(Nval);
n = 20;    %number of groups
sigma = .1;
flatrank = n*k*(k+1)/2;

%prepare dataset
subs_true = zeros(K,k,n);
for i = 1:n
    subs_true(:,:,i) = orth(randn(K,k));
end

ntries = 3;

error = zeros(nvals, 1);
momtime = zeros(1, nvals);
spmtime = zeros(1, nvals);

for j = 1:nvals
    N = Nval(j)

    errorj = zeros(ntries, 1); 
    momtimej = zeros(ntries, 1);
    spmtimej = zeros(ntries, 1);
    
    for tries = 1:ntries
        X = zeros(K,N,i);
        for i = 1:n
            X(:,:,i) = subs_true(:,:,i)*randn(k,N);
        end
        X = reshape(X,K,N*n);
        X = X(:,randperm(N*n));

        X_noisy = X + sigma*randn(K,N*n);
        
        timer = tic;
        % Estimate moments from data
        % Mest2 (Covariance) is estimated to debias Mest4
        [Mest2, Mest4] = estimate_even_moments(X_noisy);
        sigma_ = estimate_sigma_4thmoment(Mest2, Mest4);
        [Mest2, Mest4] = debias_even_moments(K, sigma_, Mest2, Mest4);
        momtimej(tries) = toc(timer);
        
        options.ntries = 3;
        
        [subs_est, ~] = subspace_power_method_general(Mest4,K,4,flatrank,k,options);
        
        spmtimej(tries) = toc(timer) - momtimej(tries);
        
        subs_est = reshape(cell2mat(subs_est),K,k,n);
        
        [~, ~, ~, errorj(tries)] = ...
            align_to_reference_multiPCA(subs_est, subs_true);
        
        errorj(tries) = sqrt(errorj(tries));

    end
    
    error(j) = mean(errorj);
    momtime(j) = mean(momtimej);
    spmtime(j) = mean(spmtimej);
end

% Save results
clear T X X_noisy
save results/test_gpca_SPM.mat

% Plot results
hf = figure(1);
ha = loglog(Nval,error,'LineWidth',2);
ha = gca;
ha.FontSize = 14;
xlabel('$N$','Interpreter','latex');
ylabel('Error (7.3)','Interpreter','latex');
pdfprint('results/test_gpca_SPM_error', hf);

hf = figure(2);
ha = loglog(Nval,[momtime;spmtime;spmtime+momtime]','LineWidth',2);
ha(2).LineStyle = ':';
ha(3).LineStyle = '--';
ha = gca;
ha.FontSize = 14;
xlabel('$N$','Interpreter','latex');
ylabel('Time (sec)','Interpreter','latex');
legend({'Estimating Moment','SPM','Total Time'});
pdfprint('results/test_gpca_SPM_time', hf);
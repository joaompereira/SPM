% (4,1) partially symmetric case
clearvars
clc
addpath('../')
addpath(pwd+"/algorithms")
addpath("../helper_functions")

nvals = 40;

rng(1,"twister")

Algs = {
    'MSPM' ,@(T, R) multiSPM(T,'rank',R,'symmetries',[1,1,1,1,2]);...
    'MSPM asym-avg',@(T,R) MSPM_11111_sym(T,R);
    'NLS',@(T,R) tensorlab_nls_41(T,R);...
    'NLS asym-avg',@(T,R) tensorlab_nls_11111_avg(T,R);...
    'MINF', @(T,R) tensorlab_minf_41(T,R);...
    'MINF asym-avg',@(T,R) tensorlab_minf_11111_avg(T,R);...
    'ALS',@(T,R) tensorlab_als_41(T,R);...
    }; 

time = zeros(nvals,size(Algs,1));
factorcos = zeros(nvals,size(Algs,1));
logerror = zeros(nvals,size(Algs,1));

for i=1:nvals
    r = 50;
    dims = [25,10];
    noise = 1/(100);

    true_factors = cell(1, 2);
    for j=1:2
        M = randn(dims(j), r);
        true_factors{j} = M ./ vecnorm(M);
    end
    
    true_lambda = exp(2*rand(1,r)-1);
    T_true = generate_lowrank_tensor(true_lambda,true_factors{:}, [4,1]);
    randT = randn(size(T_true));
    noise_const = norm(reshape(T_true,[],1))*noise/norm(reshape(randT,[],1));
    T = T_true + noise_const*randT;


    for l=1:size(Algs,1)
        tic
        [lambda,factors] = Algs{l,2}(T, r);
        time(i, l) = toc;
        T_est = generate_lowrank_tensor(lambda,factors{:}, [4,1]);
        err = norm(T-T_est,'fro');
        factorcos(i,l) = 1/2*(norm_reorder_cosine_sim(true_factors{1},factors{1}) +norm_reorder_cosine_sim(true_factors{end},factors{end})) ;
        logerror(i, l) = log10(err);

    end
end


filename = 'results/compare_41_fixed_size_25_10_50';
save(filename, "logerror","factorcos","time")


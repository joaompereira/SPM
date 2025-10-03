% three way cp decomposition test
clearvars
clc

addpath('../')
addpath(pwd+"/algorithms")
addpath("../helper_functions")


dim_vals = repmat(100,1,40);

nvals = size(dim_vals,2);


Algs = {
    'MSPM' ,@(T, R) MSPM_asym(T,'rank',R);...
    'ALS',@(T,R) tensorlab_als_111(T,R);...
    'SD',@(T,R) tensorlab_sd_111(T,R);...
    'SGSD',@(T,R) tensorlab_sgsd_111(T,R);...
    'Jennrich',@(T,R) tensorlab_gevd_111(T,R);...
    'NLS', @(T,R) tensorlab_nls_111(T,R);...
    'MINF', @(T,R) tensorlab_minf_111(T,R);
    };

time = zeros(nvals,size(Algs,1));
factorcos = zeros(nvals,size(Algs,1));
logerror = zeros(nvals,size(Algs,1));

rng(1,'twister')

for i=1:nvals
    n = dim_vals(i);
    r = 90;
    dims = [n,n,n];
    noise = 1/(100);
    
    true_factors = cell(1, 3); 
    for j=1:3
        M = randn(dims(j), r);
        true_factors{j} = M ./ vecnorm(M);
    end

    true_lambda = exp(2*rand(1,r)-1);
    [T_true, ~] = generate_lowrank_tensor(true_lambda,true_factors{:},[1,1,1]);
    randT = randn(size(T_true));
    noise_const = norm(reshape(T_true,[],1))*noise/norm(reshape(randT,[],1));
    T = T_true + noise_const*randT;

    for l=1:size(Algs,1)
        tic
        [lambda, factors] = Algs{l,2}(T, r);
        time(i, l) = toc;
        factorcos(i,l) = 1/3*(norm_reorder_cosine_sim(true_factors{1},factors{1}) +norm_reorder_cosine_sim(true_factors{2},factors{2})+norm_reorder_cosine_sim(true_factors{3},factors{3})) ;
        T_est = generate_lowrank_tensor(lambda,factors{:}, [1,1,1]);
        err = norm(T-T_est,'fro');
        logerror(i, l) = log10(err);

    end
end


filename = 'compare_111_100*3_80';
save(filename, "logerror","factorcos","time")




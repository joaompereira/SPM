% (2,1) nonorthogonal comparison 

clearvars
clc
addpath(pwd)
addpath(pwd+"/algorithms")
addpath(pwd+"/helper_functions")

rng(1,"twister")

dim_vals = repmat(100,1,40);

nvals = size(dim_vals,2);

Algs = {
    'PSSPM' ,@(T, R) spm_21sym(T, R);...
    'NLS',@(T,R) tensorlab_nls(T,R);...
    'MINF', @(T,R) tensorlab_minf(T,R);...
    'ALS',@(T,R) tensorlab_als(T,R);...
    'FFDIAG', @(T, R) ffdiag_nonortho(T);...
    'Jennrich',@(T,R) jennrich_21(T,R);...
    'HOGSVD', @(T,R) hogsvd(T,R);
    'QRJ1D', @(T,R) qrj1d(T);... 
    };

time = zeros(nvals,size(Algs,1));

logerror = zeros(nvals,size(Algs,1));

Ascore = zeros(nvals,size(Algs,1));

for i=1:nvals
    n = dim_vals(i);
    k = floor(n/2);
    r = n-20;
    dims = [n,k];
    noise = 1/(100);

    true_factors = cell(1, 2);
    for j=1:2
        M = randn(dims(j), r);
        true_factors{j} = M ./ vecnorm(M);
    end
  

    true_lambda = exp(2*rand(1,r)-1);
    [T, ~] = generate_lowrank_tensor(true_lambda,true_factors{:}, [2,1]);
    randT = randn(size(T));
    noise_const = norm(reshape(T,[],1))*noise/norm(reshape(randT,[],1));
    T = T + noise_const*randT;


    for l=1:size(Algs,1)
        tic
        [A,B, ~] = Algs{l,2}(T, r);

        lambda = vecnorm(B);
        B = B./vecnorm(B);

        time(i, l) = toc;

        Ascore(i,l) = norm_reorder_cosine_sim(true_factors{1},A);
        factors_est = {A,B};
        TF = generate_lowrank_tensor(lambda, factors_est{:}, [2,1]);
        logerror(i, l) = log10(norm(reshape(T-TF,[],1))/norm(reshape(T,[],1)));

    end
end


filename = 'compare_21_nonorthogonal_fix_size';
save(filename, "logerror","Ascore","time")





addpath '../'
addpath '../helper_functions/'
addpath(genpath('../other_packages/'))

%x_axis = 10:2:30;
dim_vals = [200, 200];
symmetries = [2, 1];
rank_vals = 100;%round(x_axis.^2 / 2)';
noise_vals = 0;

%% Broadcasting
dim_vals = dim_vals + 0*rank_vals + 0*noise_vals;
rank_vals = rank_vals + 0*dim_vals(:,1);
noise_vals = noise_vals + 0*dim_vals(:,1);

nvals = size(dim_vals,1);
order = size(dim_vals,2);

Algs = {
    'partsymSPM' ,@(T, R, symvec) partsym_SPM(T, R, symmetries=symvec);...
    'sym21SPM' ,@sym21SPM_caller;...
    %'Tensorlab', @(T, R) cpd(T, R);...
    };

time = zeros(nvals,size(Algs,1));

log10error = zeros(nvals,size(Algs,1));

for i=1:nvals
    
    dims = dim_vals(i,:);
    rank = rank_vals(i);
    noise = noise_vals(i);

    true_factors = cell(1, order);
    for k=1:order
        M = randn(dims(k), rank);
        true_factors{k} = M ./ vecnorm(M);
    end
    true_lambda = exp(2*rand(1,rank)-1);

    [T, symvec] = generate_lowrank_tensor(true_factors{:}, symmetries);
    T = T + noise*randn(size(T))/sqrt(numel(T));

    for k=1:size(Algs,1)

        tic

        [lambda, factors_est] = Algs{k,2}(T, rank, symvec);

        time(i, k) = toc;
        
        TF = generate_lowrank_tensor(lambda, factors_est{:}, symmetries);
        
        log10error(i, k) = log10(norm(reshape(T-TF,[],1))/norm(reshape(T,[],1)));

    end
end

if size(Algs,1) == 1 && nvals == 1
    dims
    rank
    time
    log10error
    
elseif nvals == 1
    
    cat = categorical(Algs(:,1));
    
    figure(1)

    stem(cat,shiftdim(time));

    figure(2)

    stem(cat,shiftdim(log10error));

else  

    figure(1)

    semilogy(x_axis, time);

    legend(Algs(:,1))

    figure(2)

    plot(x_axis, log10error,' x');

    legend(Algs(:,1))

end

function [lambda, factors_est] = sym21SPM_caller(T, R, ~) 
    [A, B] = spm_21sym(T, R);
    lambda = ones(1, R);
    factors_est = {A, B};
end 
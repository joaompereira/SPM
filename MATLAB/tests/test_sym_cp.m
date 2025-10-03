addpath '../'
addpath '../helper_functions/'
addpath(genpath('../other_packages/'))

%x_axis = 10:2:30;
dim_vals = 40;
order = 4;
rank_vals = floor(40^2/3);%round(x_axis.^2 / 2)';
noise_vals = 0;

%% Broadcasting
dim_vals = dim_vals + 0*rank_vals + 0*noise_vals;
rank_vals = rank_vals + 0*dim_vals;
noise_vals = noise_vals + 0*dim_vals;

nvals = size(dim_vals,1);

Algs = {
    'MSPM' ,MSPM_sym_wrapper(order);...
    'SPM',@(T, R) subspace_power_method(T, [], order, R);...
    };

time = zeros(nvals,size(Algs,1));

logerror = zeros(nvals,size(Algs,1));

for i=1:nvals
    
    dims = dim_vals(i);
    rank = rank_vals(i);
    noise = noise_vals(i);

    M = randn(dims, rank);
    true_factors = M ./ vecnorm(M);
    true_lambda = exp(2*rand(1,rank)-1);

    T = generate_lowrank_tensor(true_lambda, true_factors, order);
    T = T + noise*randn(size(T));

    for k=1:size(Algs,1)

        tic

        [A_est, lambda] = Algs{k,2}(T, rank);

        time(i, k) = toc;
        
        TF = generate_lowrank_tensor(lambda, A_est, order);
        
        logerror(i, k) = log10(norm(reshape(T-TF,[],1))/norm(reshape(T,[],1)));

    end
end

if size(Algs,1) == 1 && nvals == 1
    dims
    rank
    time
    logerror
    
elseif nvals == 1
    
    cat = categorical(Algs(:,1));
    
    figure(1)

    stem(cat,shiftdim(time));

    figure(2)

    stem(cat,shiftdim(logerror));

else  

    figure(1)

    semilogy(x_axis, time);

    legend(Algs(:,1))

    figure(2)

    plot(x_axis, logerror,' x');

    legend(Algs(:,1))

end

function hf = MSPM_sym_wrapper(order)
    hf = @inner;
    function [A_est, lambda] = inner(T, R)
        [lambda, A_est] = multiSPM(T, R, symmetries=ones(1, order), gradtol=1e-14);
        A_est = A_est{1};
    end
end

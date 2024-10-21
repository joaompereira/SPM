addpath '../'
addpath '../helper_functions/'
addpath(genpath('../other_packages/'))

%x_axis = 10:2:30;
dim_vals = [5, 5, 6, 6, 7, 7];
rank_vals = 10;%round(x_axis.^2 / 2)';
noise_vals = 0;

%% Broadcasting
dim_vals = dim_vals + 0*rank_vals + 0*noise_vals;
rank_vals = rank_vals + 0*dim_vals(:,1);
noise_vals = noise_vals + 0*dim_vals(:,1);

nvals = size(dim_vals,1);
order = size(dim_vals,2);

Algs = {
    'SPM v1' ,@(T, R) asym_SPM(T, R);...
    %'Tensorlab', @(T, R) cpd(T, R);...
    };

time = zeros(nvals,size(Algs,1));

logerror = zeros(nvals,size(Algs,1));

for i=1:nvals
    
    dims = dim_vals(i,:);
    rank = rank_vals(i);
    noise = noise_vals(i);

    true_factors = cell(1, order);
    for k=1:order
        true_factors{k} = randn(dims(k), rank);
    end

    T = generate_lowrank_tensor(true_factors{:});
    T = T + noise*randn(size(T));

    for k=1:size(Algs,1)

        tic

        [factors_est] = Algs{k,2}(T, rank);

        time(i, k) = toc;
        
        TF = generate_lowrank_tensor(factors_est{:});
        
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

clearvars
clc
addpath '../'
addpath '../helper_functions/'

rng_seed = 0; % for reproducibility
rng(rng_seed);

disp('Running power method iteration experiment (Figure 4)...');

L = 20;
n = 2;
vecR = 120:200;
atries = 10000;
maxiter = 5000;
f_thresholds = [0.9, 0.95, 0.99];
nfs = length(f_thresholds);
x_norm_threshs = [1e-4, 1e-10];
nxns = length(x_norm_threshs);

freq = zeros(size(vecR, 2), 4);
f_hit_it_avgs = zeros(size(vecR, 2), nfs);
f_hit_it_stds = zeros(size(vecR, 2), nfs);

xn_hit_it_avgs = zeros(size(vecR, 2), nxns);
xn_hit_it_stds = zeros(size(vecR, 2), nxns);

cn = sqrt((n-1)/n);
                        
for Riter=1:length(vecR)
    
R = vecR(Riter)

freqa = zeros(4, 1);
f_hit_it_sum = 0;
f_hit_it2_sum = 0;
xn_hit_it_sum = 0;
xn_hit_it2_sum = 0;
                        
parfor aiter=1:atries

    a_true = randn(L,R);
    a_true = a_true./vecnorm(a_true);

    decG = decomposition((a_true'*a_true).^n);

    Ak = randn(L,1);
    Ak = Ak/norm(Ak);

    f_save = 0;
    
    converged = false;
    f_ind = 1;
    xn_ind = 1;

    f_hit_it = maxiter + zeros(1, nfs);
    xn_hit_it = maxiter + zeros(1, nxns);

    for tries = 1:maxiter

        Ax = (Ak'*a_true)';

        Ak_new = a_true * (Ax.^(n-1) .* (decG \ (Ax.^n)));

        f = Ak_new'*Ak;

        if f_ind<=nfs && f > f_thresholds(f_ind)
            f_hit_it(f_ind) = tries;
            f_ind = f_ind + 1;
        end

        % Determine optimal shift
        % Sometimes due to numerical error f can be greater than 1
        clambda = 1;
        shift = cn*clambda;

        % Shifted power method
        Ak_new = Ak_new + shift*Ak;
        Ak_new = Ak_new/norm(Ak_new);
        l2norm = norm(Ak - Ak_new);

        if xn_ind <= nxns && l2norm < x_norm_threshs(xn_ind)
            xn_hit_it(xn_ind) = tries;
            xn_ind = xn_ind + 1;
        end

        if l2norm < 1e-10
            % Algorithm converged
            converged = true;
            Ak = Ak_new;
            break
        else
            Ak = Ak_new;
        end
    end

    err = 2*min(1-abs(Ak'*a_true));
    
    if ~converged
       j = 4
    elseif err < 1e-10
       j = 1;
    elseif f > 1 - 1e-10
       j = 3;
    else
       j = 2;
    end
 
    freqa = freqa + (1:4 == j)';
    f_hit_it_sum = f_hit_it_sum + f_hit_it;
    f_hit_it2_sum = f_hit_it2_sum + f_hit_it.^2;
    xn_hit_it_sum = xn_hit_it_sum + xn_hit_it;
    xn_hit_it2_sum = xn_hit_it2_sum + xn_hit_it.^2;

end

freq(Riter, :) = freqa/atries;
f_hit_it_avgs(Riter, :) = f_hit_it_sum / atries;
f_hit_it_stds(Riter, :) = sqrt(f_hit_it2_sum / atries - f_hit_it_avgs(Riter, :).^2);
xn_hit_it_avgs(Riter, :) = xn_hit_it_sum / atries;
xn_hit_it_stds(Riter, :) = sqrt(xn_hit_it2_sum / atries - xn_hit_it_avgs(Riter, :).^2);

end

clear decG
save results/test_power_method.mat
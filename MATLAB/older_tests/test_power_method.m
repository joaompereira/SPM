clearvars
clc
addpath '../'
addpath '../helper_functions/'

rng_seed = 4; % for reproducibility
rng(rng_seed);

disp('Running power method iteration experiment (Figure 4)...');

L = 20;
n = 2;
vecR = 120:200;
atries = 10000;
maxiter = 5000;
repetitions = 1;

freq = zeros(size(vecR, 1), 4);
cn = sqrt((n-1)/n);
                        
for Riter=1:length(vecR)
    
R = vecR(Riter)

freqa = zeros(4, 1);
                        
parfor aiter=1:atries

    a_true = randn(L,R);
    a_true = a_true./vecnorm(a_true);

    decG = decomposition((a_true'*a_true).^n);

    Ak = randn(L,1);
    Ak = Ak/norm(Ak);

    f_save = 0;
    
    for rep = 1:repetitions
        converged = false;
        for tries = 1:maxiter
    
            Ax = (Ak'*a_true)';
    
            Ak_new = a_true * (Ax.^(n-1) .* (decG \ (Ax.^n)));
    
            f = Ak_new'*Ak;
    
            % Determine optimal shift
            % Sometimes due to numerical error f can be greater than 1
            clambda = 1;
            shift = cn*clambda;
    
            % Shifted power method
            Ak_new = Ak_new + shift*Ak;
            Ak_new = Ak_new/norm(Ak_new);
    
            if norm(Ak - Ak_new) < 1e-13
                % Algorithm converged
                converged = true;
                Ak = Ak_new;
                break
            else
                Ak = Ak_new;
            end
        end
        
        if converged && (f > 1 - 1e-10)
            break
        elseif rep == 1 || (converged && (f > f_save))
            f_save = f
            Ak_save = Ak
            conv_save = converged            
        elseif rep == repetitions
            f = f_save
            Ak = Ak_save
            converged = conv_save
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

end

freq(Riter, :) = freqa/atries;

end

clear decG
save results/test_power_method.mat
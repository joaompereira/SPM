clearvars
clc
addpath '../'
addpath '../helper_functions/'
addpath(genpath('../other_packages/'))

rng_seed = 9; % for reproducibility
rng(rng_seed);

disp('Running noise experiment (Figure 3)...');


% L and R values
L = 15;
R = 30;
n = 4;

% Different Algorithms
Algs = {'SPM', @(T,L,R) subspace_power_method(T, L, n, R);
        'FOOBI', @(T,L,R) FOOBI1_tensor_decomposition(T, L, R);
        'Tensorlab', @(T,L,R) tensorlab_ccpd_nls(T, L, n, R);       
        };

nalgs = size(Algs,1);

ntries = 5;

sigmavals = 10.^[-4:-2, -1:.1:0];%linspace(-4, 0, 13);%
nvals = length(sigmavals);

% True rank decomposition
A_true = randn(L,R) / sqrt(L);
T_clean = generate_lowrank_tensor(A_true, n);
%gamma = norm(T_clean(:)) / sqrt(L^3 / 6);

% Time of execution and estimate error for different algorithms and (L,R)
% values

time = NaN(nalgs,nvals);
time_el = NaN(nalgs,nvals);
std_time = NaN(nalgs,nvals);
std_time_el = NaN(nalgs,nvals);
freq_converged = NaN(nalgs,nvals);
min_err = NaN(nalgs,nvals);

warning('error', 'MATLAB:rankDeficientMatrix');


for j=1:nvals
    sigma = sigmavals(j);
    fprintf('sigma =%2e\n', sigmavals(j));

    % Generate nD tensor
    T_noise = symmetrize_tensor(randn(L, L, L));
    T_noise = T_noise / norm(reshape(T_noise, [], 1));
    T = T_clean + sigma * T_noise;

    for i=1:nalgs

        if j>1 && isnan(min_err(i,j-1))
          continue
        end

        fprintf(' - %s\n',Algs{i});

        % Pre-allocate local time and error vectors for all tries
        timeij = nan(ntries, 1);
        errij = nan(ntries,1);

        for tries=1:ntries
        
            timer = tic;

            error_flag = true;
            
            try
                % Use Algorithm i to get rank decomposition
                A_est = Algs{i,2}(T,L,R);
            catch ME
                continue
            end
    
            % Get time
            timeij(tries) = toc(timer);

            % Calculate error
            errij(tries) = rderror(A_est, A_true, n);
            
        end


        % Average over all runs
        time(i,j) = mean(timeij(~isnan(timeij)));
        time_el(i,j) = mean(log(timeij(~isnan(timeij))));
        std_time(i,j) = std(timeij(~isnan(timeij)));
        std_time_el(i,j) = std(log(timeij(~isnan(timeij))));
        threshold = 8*median(errij(~isnan(errij)));
        freq_converged(i,j) = mean(errij < threshold);
        %mean(errij(errij < threshold));
        min_err(i, j) = min(errij);
        
    end
    
end

% Save results
clear T T_noise T_clean % T is cleared since it occupies a lot of space
save results/testrd_noise.mat




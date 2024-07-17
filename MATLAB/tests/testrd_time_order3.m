clearvars
clc
addpath '../'
addpath '../helper_functions/'
addpath(genpath('../other_packages/'))

rng_seed = 17; % for reproducibility
rng(rng_seed);

disp('Running computation time experiment (order 3, Figure 2a)...');

% L and R values
Lval = unique(round(exp(linspace(log(10),log(300),13))));
Rval = Lval;
nvals = length(Lval);
n = 3;

% Different Algorithms
Algs = {'SPM', @(T,L,R) subspace_power_method(T, L, n, R);
        'LRSTA', @(T,L,R) GPSTDAPRX(T, R); 
        'Tensorlab', @(T,L,R) tensorlab_ccpd_nls(T, L, n, R);
        };

nalgs = size(Algs,1);

% Number of tries for each algorithm
% The variance of the running time is bigger
% for some algorithms. We run these more
% times to produce smoother results
ntries = [5,5,100];

% Time of execution and estimate error for different algorithms and 
% (L,R) values
time = cell(nalgs, 1);
error = cell(nalgs, 1);

for i=1:nalgs
    time{i} = NaN(nvals, ntries(i));
    error{i} = NaN(nvals, ntries(i));
end

%% 1st Call
% Because of MATLAB's inner workings, the
% 1st call to a function has overhead.
% Therefore, we call the algorithms for the 
% first time without recording timings.
L = Lval(1);
R = Rval(1);

% True rank decomposition
A_true = randn(L,R) / sqrt(L);

% Generate nD tensor
T = generate_lowrank_tensor(A_true, n);

s = rng;

for i=1:nalgs

    % Use Algorithm i to get rank decomposition
    A_est = Algs{i,2}(T,L,R);

    rng(s);
end

%% The actual runs
for j=1:nvals
    L = Lval(j);
    R = Rval(j);
    fprintf('L =%3d\n',L);
   
    % True rank decomposition
    A_true = randn(L,R) / sqrt(L);

    % Generate nD tensor
    T = generate_lowrank_tensor(A_true, n);

    s = rng;

    for i=1:nalgs

        if j>1 && ~(mean(time{i}(j-1,:),'omitnan') <= 100)
          continue
        end

        fprintf(' - %s\n',Algs{i});

        for tries=1:ntries(i)
        
            timer = tic;
            
            % Use Algorithm i to get rank decomposition
            A_est = Algs{i,2}(T,L,R);
    
            % Get time
            time{i}(j, tries) = toc(timer);

            % Calculate error
            error{i}(j, tries) = rderror(A_est, A_true, n);
            
            
        end

        rng(s);
    end
    
end

% Save results
clear T % T is cleared since it occupies a lot of space
save results/testrd_time_order3.mat
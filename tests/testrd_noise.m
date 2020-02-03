clear all
clc
addpath '../'
addpath '../helper_functions/'
%addpath(genpath('../other_packages/'))

% L and R values
L = 15;
R = floor(L^2/3);

% Different Algorithms
% The tensor order has to be 4 for FOOBI to work
Algs = {'FOOBI1',@(T,L,R) FOOBI1_tensor_decomposition(T, L, R);...
        'SPM' ,@(T,L,R) subspace_power_method(T, L, 4, R);...
        %'TensorLab' ,@(T,L,R) tensorlab_ccpd_nls(T,L,4,R);...  
        };

nalgs = size(Algs,1);

noise = 10.^(linspace(-8,-2,4));
nvals = length(noise);

% Number of tries for each algorithm
% Some algorithms have more variable running times, so we run them more
% times to produce smoother results
ntries = [5,5,5];

% True rank decomposition
A_true = randn(L,R);

% Generate clean 4D tensor
T_clean = generate_lowrank_tensor(A_true, 4);

% Time of execution and estimate error for different algorithms and (L,R)
% values
time = NaN(nalgs,nvals);
error = NaN(nalgs,nvals);

s = "-";
for i = 2:45
  s = s + "-";
end

fprintf('\n%s\n',s)
fprintf('| Noise |')
for i=1:nalgs
  fprintf(' %9s |',Algs{i,1});
end
fprintf('\n%s\n',s)
   
for j=1:nvals
    
    fprintf('| %1.0e |', noise(j))
    % Add noise to entries of Tensor, making sure the tensor is still
    % symmetric
    T = T_clean + noise(j)*symmetrize_tensor(randn(L*ones(1,4)));
        
    for i=1:nalgs
        
        % Pre-allocate local time and error vectors for all tries
        timeij = zeros(ntries(i),1);
        errorij = zeros(ntries(i),1);
        
        for tries=1:ntries(i)
            
            % Start clock
            timer = tic;
            
            % Use Algorithm i to get rank decomposition
            A_est = Algs{i,2}(T,L,R);
            
            % Get time
            timeij(tries) = toc(timer);
            
            % Calculate error
            errorij(tries) = rderror(A_est, A_true, 4);
            
        end
        
        % Log Average over all runs
        time(i,j) = mean(timeij(errorij<1e-1));
        error(i,j) = mean(errorij(errorij<1e-1));
        
        fprintf(' %1.3e |',error(i,j));

    end
    
    fprintf('\n%s\n',s)
end

ratio = mean(error./noise,2);
fprintf('| Ratio |')
for i=1:nalgs
  fprintf(' %1.3e |',ratio(i));
end
fprintf('\n%s\n',s)

% Save results
clear T T_clean % T is cleared since it occupies a lot of space
save results/testrd_noise.mat



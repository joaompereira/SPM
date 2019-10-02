clear all
clc

% L and R values
Lval = 9:55;
Rval = floor(Lval.^2/3);
nvals = length(Lval);

% Different Algorithms
% The tensor order has to be 4 for FOOBI to work
Algs = {'FOOBI1',@(T,L,R) FOOBI1_tensor_decomposition(T, L);...
        'SPM' ,@(T,L,R) subspace_power_method(T, L, 4);...
        'TensorLab' ,@(T,L,R) tensorlab_ccpd_nls(T,L,4,R);...  
        };

nalgs = size(Algs,1);

% Number of tries for each algorithm
% Some algorithms have more variable running times, so we run them more
% times to produce smoother results
ntries = [3,3,11];

% Time of execution and estimate error for different algorithms and (L,R)
% values
time = NaN(nalgs,nvals);
error = NaN(nalgs,nvals);

for i=1:nalgs
    
    fprintf('\n%s:\n',Algs{i});
    
    for j=1:nvals
        L = Lval(j);
        R = Rval(j);
        fprintf('L =%3d\n',L);        
        
        % True rank decomposition
        A_true = randn(L,R);

        % Generate 4D tensor
        T = generate_lowrank_tensor(A_true, 4);
        
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
            errorij(tries) = rank_decomposition_error(A_est, A_true);
            
        end
        
        % Log Average over all runs
        time(i,j) = exp(mean(log(timeij)));
        error(i,j) = exp(mean(log(errorij)));
        
        if time(i,j) > 100
          break
        end
    end
end

% Plot results
loglog(Lval,time,'LineWidth',1);
legend(Algs(:,1));
xlabel('$L$','Interpreter','latex');
ylabel('time(s)','Interpreter','latex');

% Save results
clear T % T is cleares since it occupies a lot of space
save('rank_decomposition_time.mat')
fprintf('Execution finished');
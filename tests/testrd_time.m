clearvars
clc
addpath '../'
addpath '../helper_functions/'
%addpath(genpath('../other_packages/'))

% L and R values
Lval = round(exp(linspace(log(9),log(60),23)));
Rval = floor(Lval.^2/3);
nvals = length(Lval);

% Different Algorithms
% The tensor order has to be 4 for FOOBI to work
Algs = {'SPM' ,@(T,L,R) subspace_power_method(T, L, 4);
        'FOOBI',@(T,L,R) FOOBI1_tensor_decomposition(T, L);
        %'Tensorlab' ,@(T,L,R) tensorlab_ccpd_nls(T, L, 4, R);
        };

nalgs = size(Algs,1);

% Number of tries for each algorithm
% Some algorithms have more variable running times, so we run them more
% times to produce smoother results
ntries = [6,6,60];

% Time of execution and estimate error for different algorithms and (L,R)
% values
time = NaN(nalgs,nvals);
error = NaN(nalgs,nvals);

for j=1:nvals
    L = Lval(j);
    R = Rval(j);
    fprintf('L =%3d\n',L);        

    % True rank decomposition
    A_true = randn(L,R);

    % Generate 4D tensor
    T = generate_lowrank_tensor(A_true, 4);
        
    for i=1:nalgs
        
        if j>1 && (time(i,j-1) > 100 || isnan(time(i,j-1)))
          continue
        end

        fprintf(' - %s\n',Algs{i});
        
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
        time(i,j) = mean(timeij(errorij<1e-6));
        error(i,j) = mean(errorij(errorij<1e-6));
        
    end
end

% Save results
clear T % T is cleared since it occupies a lot of space
save results/testrd_time.mat

% Plot results
hf = figure(1);
ha = loglog(Lval(2:end),time(:,2:end),'LineWidth',2);
ha(2).LineStyle = ':';
ha(3).LineStyle = '--';
legend(Algs(:,1),'Location','southeast');
xlabel('$L$','Interpreter','latex');
ylabel('Time (sec)','Interpreter','latex');
ylim([1e-2 100]);
pdfprint('results/testrd_time', hf);


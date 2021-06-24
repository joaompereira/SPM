clearvars
clc
addpath '../'
addpath '../helper_functions/'
addpath(genpath('../other_packages/'))

% L and R values
Lval = round(exp(linspace(log(4),log(30),12)));
Rval = floor(Lval.^2/3);
nvals = length(Lval);
n = 6;

% Different Algorithms
Algs = {'SPM' ,@(T,L,R) subspace_power_method(T, L, n);
        'LRSTA' ,@(T,L,R) GPSTDAPRX(T, R);
        };

nalgs = size(Algs,1);

% Number of tries for each algorithm
% Some algorithms have more variable running times, so we run them more
% times to produce smoother results
ntries = [2,2,20];

% Time of execution and estimate error for different algorithms and (L,R)
% values
time = NaN(nalgs,nvals);
std_time = NaN(nalgs,nvals);


warning('error', 'MATLAB:rankDeficientMatrix');

for j=1:nvals
    L = Lval(j);
    R = Rval(j);
    fprintf('L =%3d\n',L);
    
    error_flag = true;
    
    while error_flag
    try
        % True rank decomposition
        A_true = randn(L,R);
        A_true = A_true./vecnorm(A_true);
        A_true = A_true.*3.^((1:R)/R-1);

        % Generate nD tensor
        T = generate_lowrank_tensor(A_true, n);

        for i=1:nalgs

            if j>1 && (time(i,j-1) > 100 || isnan(time(i,j-1)))
              continue
            end

            fprintf(' - %s\n',Algs{i});

            % Pre-allocate local time and error vectors for all tries
            timeij = zeros(ntries(i),1);
            errorij = zeros(ntries(i),1);

            for tries=1:ntries(i)
                
                error = Inf;
                while error > 1e-6
                timer = tic;

                % Use Algorithm i to get rank decomposition
                A_est = Algs{i,2}(T,L,R);

                % Get time
                timeij(tries) = toc(timer);

                % Calculate error
                error = rderror(A_est, A_true, n);
                end

            end

            % Average over all runs
            time(i,j) = mean(timeij);
            std_time(i,j) = std(timeij);

        
        end
        error_flag = false;
    end
    end
    
end

% Save results
clear T % T is cleared since it occupies a lot of space
save results/testrd_time_LRSTA.mat

% Plot results
hf = figure(1);
ha = loglog(Lval(2:end),time(:,2:end),'LineWidth',2);
ha(2).LineStyle = ':';
%ha(3).LineStyle = '--';
legend(Algs(:,1),'Location','southeast');
xlabel('$L$','Interpreter','latex');
ylabel('Time (sec)','Interpreter','latex');
ylim([1e-3 100]);
pdfprint('results/testrd_time_LRSTA', hf);


clearvars
clc
addpath '../'
addpath '../helper_functions/'

rng_seed = 9; % for reproducibility
rng(rng_seed);

shortfields = {'n',
               'L ',
               ' R ',
               'exttime';
               'powtime';
               'deftime';
               'tottime';
               'avgit';
               'nrr';
               '  error   '};

fields = {'n',
          'L',
          'R',
          'extracttime';
          'powertime';
          'deflatetime';
          'totaltime';
          'avgiter';
          'nrr';
          'error'};

nfields = length(fields);
         
disp('Running table experiment (Table 1)...');



%% WARM-UP
% calling a function for the first time in MATLAB takes longer than the
% following calls
row = 0;
fprintf('   Row  ', row);

L = 40;
R = 200;
n = 4;

% Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n, R);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

%% ROW 1 %%
row = row+1;
fprintf('\b%d', row);

% Tensor characteristics 
L = 136;     % Length
R = 136;    % Rank
n = 3;      % Order

% Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n, R);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

% Collect more stats
stat_.L = L;
stat_.R = R;
stat_.n = n;

stat(row) = stat_;

%% ROW 2 %%
row = row+1;
fprintf('\b%d', row);

% Tensor characteristics 
L = 40;     % Length
R = 200;    % Rank
n = 4;      % Order

% Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n, R);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

% Collect more stats
stat_.L = L;
stat_.R = R;
stat_.n = n;

stat(row) = stat_;

%% ROW 3 %%
row = row+1;
fprintf('\b%d', row);

% Tensor characteristics 
L = 40;     % Length
R = 400;    % Rank
n = 4;      % Order

% Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n, R);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

% Collect more stats
stat_.L = L;
stat_.R = R;
stat_.n = n;

stat(row) = stat_;

%% ROW 4 %%
row = row+1;
fprintf('\b%d', row);

% Tensor characteristics 
L = 40;     % Length
R = 600;    % Rank
n = 4;      % Order

% Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n, R);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

% Collect more stats
stat_.L = L;
stat_.R = R;
stat_.n = n;

stat(row) = stat_;

%% ROW 5 %%
row = row+1;
fprintf('\b%d', row);

% Tensor characteristics 
L = 60;     % Length
R = 400;    % Rank
n = 4;      % Order

% Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n, R);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

% Collect more stats
stat_.L = L;
stat_.R = R;
stat_.n = n;

stat(row) = stat_;

%% ROW 6 %%
row = row+1;
fprintf('\b%d', row);

% Tensor characteristics 
L = 80;     % Length
R = 400;    % Rank
n = 4;      % Order

% Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n, R);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

% Collect more stats
stat_.L = L;
stat_.R = R;
stat_.n = n;

stat(row) = stat_;

%% Row 7 %%
row = row+1;
fprintf('\b%d', row);

% Tensor characteristics 
L = 40;     % Length
R = 200;    % Rank
n = 4;      % Order

% Generate true rank decomposition
A_true = randn(L,R) + ones(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n, R);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

% Collect more stats
stat_.L = L;
stat_.R = R;
stat_.n = n;

stat(row) = stat_;

%% ROW 8 %%
row = row+1;
fprintf('\b%d', row);

% Tensor characteristics 
L = 40;     % Length
R = 200;    % Rank
n = 4;      % Order

% Generate true rank decomposition
A_true = abs(randn(L,R));
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n, R);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

%stat_.error
%end

% Collect more stats
stat_.L = L;
stat_.R = R;
stat_.n = n;

stat(row) = stat_;

%% ROW 9 %%
row = row+1;
fprintf('\b%d', row);

% Tensor characteristics 
L = 19;     % Length
R = 190;    % Rank
n = 5;      % Order

%  Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n, R);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

% Collect more stats
stat_.L = L;
stat_.R = R;
stat_.n = n;

stat(row) = stat_;

%% ROW 10 %%
row = row+1;
fprintf('\b%d', row);

% Tensor characteristics 
L = 11;     % Length
R = 250;    % Rank
n = 6;      % Order

%  Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n, R);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

% Collect more stats
stat_.L = L;
stat_.R = R;
stat_.n = n;

stat(row) = stat_;

fprintf('\n')

nrows = row;

%% Save results
clear T % T is cleared since it occupies a lot of space
save results/testrd_table.mat

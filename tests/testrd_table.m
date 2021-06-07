clearvars
clc
addpath '../'
addpath '../helper_functions/'
    
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
         
%% WARM-UP
% calling a function for the first time in MATLAB takes longer than the
% following calls
row = 0;
fprintf('Row  ', row);

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
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

%% ROW 1 %%
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
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n);

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
R = 400;    % Rank
n = 4;      % Order

% Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n);

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
R = 600;    % Rank
n = 4;      % Order

% Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n);

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
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n);

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
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n);

% Calculate decomposition error
stat_.error = rderror(A_est, A_true, lambda_est, lambda_true, n);

% Collect more stats
stat_.L = L;
stat_.R = R;
stat_.n = n;

stat(row) = stat_;

%% New ROW 5a %%
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
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n);

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
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n);

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
L = 40;     % Length
n = 4;      % Order 

% Generalized Decomposition
ranks = ones(20,1)*(2:3);
ranks = ranks(:);
nranks = length(ranks);
flatrank = 20*(nchoosek(2+1,2)+nchoosek(3+1,2));
A_true = cell(nranks,1);
W_true = cell(nranks,1);

% Generate true rank decomposition
for i=1:nranks
    R = ranks(i);
    A_true{i} = orth(randn(L,R));
    W_true{i} = symmetrize_tensor(randn(R*ones(1,n)));
end

% Generate tensor
T = generate_lowrank_tensor_general(A_true, W_true, n);

% Get rank decomposition
[A_est, W_est, stat_] = subspace_power_method_general(T, L, n, []);

% Calculate decomposition error
stat_.error = rderror_general(A_est, A_true, W_est, W_true, n);

% Collect more stats
stat_.L = L;
stat_.R = flatrank;
stat_.n = n;

stat(row) = stat_;

%% ROW 7 %%
row = row+1;
fprintf('\b%d', row);

% Tensor characteristics 
L = 16;     % Length
R = 400;    % Rank
n = 6;      % Order

%  Generate true rank decomposition
A_true = randn(L,R);
A_true = A_true./vecnorm(A_true);
lambda_true = randn(1,R);

% Generate tensor
T = generate_lowrank_tensor(A_true, lambda_true, n);

% Get rank decomposition
[A_est, lambda_est, stat_] = subspace_power_method(T, L, n);

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
L = 16;     % Length
n = 6;      % Order 

ranks = ones(8,1)*(2:3);
ranks = ranks(:);
nranks = length(ranks);
flatrank = 8*(nchoosek(2+2,3)+nchoosek(3+2,3));
A_true = cell(nranks,1);
W_true = cell(nranks,1);

% Generate true rank decomposition
for i=1:nranks
    R = ranks(i);
    A_true{i} = orth(randn(L,R));
    W_true{i} = symmetrize_tensor(randn(R*ones(1,n)));
end

% Generate tensor
T = generate_lowrank_tensor_general(A_true, W_true, n);

% Get rank decomposition
[A_est, W_est, stat_] = subspace_power_method_general(T, L, n, []);

% Calculate decomposition error
stat_.error = rderror_general(A_est, A_true, W_est, W_true, n);

% Collect more stats
stat_.L = L;
stat_.R = flatrank;
stat_.n = n;

stat(row) = stat_;

%% Show Table

% Print header
nrows = row;
fprintf('\n\n     ');
for j=1:nfields
    fprintf('|%s', shortfields{j});
end
fprintf('|\n');

% Print rest of the table
for i=1:nrows
    
    fprintf(' |T%-2d|%1d|%2d|%3d',i, stat(i).n, stat(i).L, stat(i).R);
    for j=4:7
        fprintf('| %5.2f ', stat(i).(fields{j}))
    end
    fprintf('| %3.0f ',stat(i).(fields{8}));
    fprintf('|%3d',stat(i).(fields{9}));
    fprintf('| %1.2e |\n',stat(i).(fields{10}));
end

%% Save results
clear T % T is cleared since it occupies a lot of space
%save results/testrd_table.mat

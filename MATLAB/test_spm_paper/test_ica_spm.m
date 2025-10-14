clearvars
clc
addpath '../'
addpath '../helper_functions/'

load ../../datasets/013-2015/Subject01_s1.mat

data = run{1}.eeg;

rng(0);

[p, n] = size(data);
t = (0:p-1) / 512;

means = mean(data);
data_m = data-means;

ncols = 1;
nrows = 4;

X = data_m';

X2 = reshape(X, n, 1, p) .* reshape(X, 1, n, p);
X2 = reshape(X2, [], p);
d = 3;

if d==4
    M2 = reshape((X * X') / p, [], 1);
    T = (X2 * X2') / p - 3*(M2 * M2');
    T = symmetrize_tensor(T, n, 4);
else
    T = (X2 * X') / p;
end

r = 64;

[A, lambda] = subspace_power_method(T, n, d, r, 'ntries', 5, 'ftol', 1e-8);

fprintf("L2 error: %f\n", norm(reshape(T, [], 1) - reshape(generate_lowrank_tensor(A, lambda, d), [], 1)) / norm(reshape(T, [], 1)))

SPM_icasig = ((A'*A + 0.01*eye(r)) \ A') * X;

[~, skew_inds] = sort(abs(skewness(SPM_icasig,[],2)), 'descend');

clear X2

save results\ica_spm.mat

plot_ica_spm





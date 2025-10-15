addpath '../'
addpath '../helper_functions/'

n = 10;
m = 6;
r = 150;

noise_std = 1e-6;

fprintf('Testing SPM...\nTensor specifications:\n')
fprintf('Dimension: %d   Rank: %d   Order: %d\n', n, r, m)
fprintf('Entrywise noise std: %.2e\n \n', noise_std)

A_true = randn(n,r) ./ sqrt(n);
lambda_true = ones(1, r);

T = generate_lowrank_tensor(A_true, lambda_true, m);
if noise_std>0
    T = T + noise_std * symmetrize_tensor(randn([n*ones(1, m), 1]), n, m);
end

timer = tic;

[A_est, lambda_est, stat] = subspace_power_method(T);

time = toc(timer);
error = rderror(A_est, A_true, lambda_est, lambda_true, m);

fprintf('Results\n')
fprintf('   time: %.2fs\n', time);
fprintf('  error: %.2e\n', error);
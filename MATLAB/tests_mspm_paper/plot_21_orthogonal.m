% plot the figure for test_21_orthogonal.m
clearvars
clc

hf = figure;

Algs = {
    'MSPM' ,@(T, R) MSPM_21sym(T, R);...
    'Jacobi', @(T, R) jacobi(T);...
    'FFDIAG', @(T, R) ffdiag_ortho(T);...
    'Jennrich',@(T,R) jennrich_21(T,R);...
    'NLS',@(T,R) tensorlab_nls(T,R);...
    'MINF', @(T,R) tensorlab_minf(T,R);...
    'ALS',@(T,R) tensorlab_als(T,R);...
    'SVD',@(T,R) svd_ortho(T,R);
    'HOGSVD',@(T,R) hogsvd(T,R);
    'NLS SVD-init',@(T,R) tensorlab_nls_better_initialization(T,R);
    };

data = load('results/compare_21_orthogonal_fix_size_100_50_80.mat');
Ascore = data.Ascore;
logerror = data.logerror;
time=data.time;

colors = [
    1,    0,    0;    % red
    0,    0,    1;    % blue
    0,  0.5,    0;    % green
    1,  0.5,    0;    % orange
    0.5,  0,  0.5;    % purple
    0,    1,    1;    % cyan
    1, 0.84,    0;    % gold
    0.3, 0.3, 0.3;    % dark gray
    0.6, 0.2, 0.2;    % brick red
    0.2, 0.6, 0.6     % muted teal
];

hold on;
for i = 1:size(Algs,1)
        h(i) = scatter(log10(1-real(Ascore(:,i))),log10(time(:,i)),40, colors(i,:),'filled');
end
legend(h,Algs(:,1))

xlabel('Log(1-Ascore)');
ylabel('Log Runtime');
hold off;

pdfprint('results/compare_21_orthogonal_fix_size_100_50_80', hf);
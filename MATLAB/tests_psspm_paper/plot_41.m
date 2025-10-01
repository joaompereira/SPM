% plot the figure for test_41.m
clearvars
clc

data = load("compare_41_fixed_size_25_10_50.mat");

Algs = {
    'PSSPM' ,@(T, R) PSSPM_41(T,'rank',R,'symmetries',[1,1,1,1,2]);...
    'PSSPM asym-avg',@(T,R) PSSPM_11111_sym(T,'rank',R);
    'NLS',@(T,R) tensorlab_nls_41(T,R);...
    'NLS asym-avg',@(T,R) tensorlab_nls_11111_avg(T,R);...
    'MINF', @(T,R) tensorlab_minf_41(T,R);...
    'MINF asym-avg',@(T,R) tensorlab_minf_11111_avg(T,R);...
    'ALS',@(T,R) tensorlab_als_41(T,R);...
    };

factorscore = data.factorcos;
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
];

hold on;
for i = 1:size(Algs,1)
        h(i) = scatter(logerror(:,i),log10(time(:,i)),40, colors(i,:),'filled');
end
legend(h,Algs(:,1))

xlabel('Log Reconstruction Error');
ylabel('Log Runtime');
hold off;


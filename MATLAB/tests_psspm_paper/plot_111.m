% plot the figure for test_111.m
clearvars
clc

hf = figure;

data = load("results/compare_111_100_3_80.mat");

Algs = {
    'MSPM' ,@(T, R) MSPM_asym(T,'rank',R);...
    'ALS',@(T,R) tensorlab_als_111(T,R);...
    'SD',@(T,R) tensorlab_sd_111(T,R);...
    'SGSD',@(T,R) tensorlab_sgsd_111(T,R);...
    'Jennrich',@(T,R) tensorlab_gevd_111(T,R);...
    'NLS', @(T,R) tensorlab_nls_111(T,R);...
    'MINF', @(T,R) tensorlab_minf_111(T,R);
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

pdfprint('results/compare_111_100_3_80', hf);


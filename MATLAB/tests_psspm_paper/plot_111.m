% plot the figure for test_111.m
clearvars
clc

data = load("compare_111_100*3_80.mat");

Algs = {
    'PSSPM' ,@(T, R) PSSPM_111(T,'rank',R);...
    'ALS',@(T,R) tensorlab_als_111(T,R);...
    'SD',@(T,R) tensorlab_sd_111(T,R);...
    'SGSD',@(T,R) tensorlab_sgsd_111(T,R);...
    'Jennrich',@(T,R) tensorlab_gevd_111(T,R);...
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
];

hold on;
for i = 1:size(Algs,1)
        h(i) = scatter(logerror(:,i),log10(time(:,i)),40, colors(i,:),'filled');
end
legend(h,Algs(:,1))

xlabel('Log Reconstruction Error');
ylabel('Log Runtime');
hold off;


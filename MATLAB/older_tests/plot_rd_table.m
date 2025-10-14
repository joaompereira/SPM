clearvars
addpath '../helper_functions/'

load results/testrd_table.mat

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
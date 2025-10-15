%% Create results directory if needed
if ~exist('results\', 'dir')
   mkdir('results')
end

%% Run Experiments
testrd_table        % Table 1
testrd_time_order3  % Figure 2a
testrd_time_order4  % Figure 2b
testrd_noise        % Figure 3
test_power_method   % Figure 4

%% Plot Experiments
generate_rd_table   % Table 1
plot_rd_time_order3 % Figure 2a
plot_rd_time_order4 % Figure 2b
plot_rd_noise       % Figure 3
plot_power_method   % Figure 4
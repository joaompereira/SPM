%% Create results directory if needed
if ~exist('results\', 'dir')
   mkdir('results')
end

%% Run Experiments
testrd_time_order3  % Figure 2a
testrd_time_order4  % Figure 2b
testrd_noise_ord3   % Figure 3a
testrd_noise_ord4   % Figure 3b
test_power_method   % Figure 4 & 5
test_ica_spm        % Figure 6

%% Plot Experiments
plot_rd_time_order3 % Figure 2a
plot_rd_time_order4 % Figure 2b
plot_rd_noise_ord3  % Figure 3a
plot_rd_noise_ord4  % Figure 3b
plot_power_method   % Figure 4 & 5
plot_ica_spm        % Figure 6

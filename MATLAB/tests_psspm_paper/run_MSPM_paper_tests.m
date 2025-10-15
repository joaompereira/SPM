%% Create results directory if needed
if ~exist('results\', 'dir')
   mkdir('results')
end

%% Run Experiments
test_21_nonorthogonal  
test_21_orthogonal
test_41
test_111

%% Plot Experiments
plot_21_nonorthogonal  
plot_21_orthogonal
plot_41
plot_111

clearvars
close(figure(1))
clf
addpath '../helper_functions/'

load results/testrd_noise.mat

noise_stability = mean(min_err(:, 1:4)./ sigmavals(1:4), 2);

% Plot results
hf = figure(1);
hf.Position = [100 100 450 300];
ax = gca;
ha = loglog(sigmavals, min_err,'LineWidth',2);
ha(2).LineStyle = ':';
ha(3).LineStyle = '--';
legends = cell(nalgs, 1);
for i=1:nalgs
    legends{i} = sprintf('%s (%.2f)', Algs{i,1}, noise_stability(i));
end
legend(legends, 'Interpreter','latex','Location','southeast');
xlabel('$\sigma$','Interpreter','latex');
ylabel('Err','Interpreter','latex');
ax.FontSize = 10;
ax.TickLabelInterpreter = "latex";
ylim([1e-4 10]);
pdfprint('results/testrd_noise', hf);
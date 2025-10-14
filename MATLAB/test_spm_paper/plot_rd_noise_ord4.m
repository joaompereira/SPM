clearvars
close(figure(1))
clf
addpath '../helper_functions/'

load results/testrd_noise.mat

noise_stability = mean(min_err(:, 1:7)./ sigmavals(:, 1:7), 2);

linestyles = ["-", "--", "-."];
markerstyles = ["o", "<", "diamond", ">", 'square', "^", "v"];

% Plot results
hf = figure(1);
hf.Position = [100 100 400 300];
ax = gca;
ha = loglog(sigmavals, min_err, 'LineWidth',1.5,'MarkerIndices',2:2:8);
for i=1:length(ha)
    ha(i).LineStyle = linestyles(mod(i-1,3)+1);
    ha(i).Marker= markerstyles(mod(i-1,7)+1);
    ha(i).MarkerSize = 4;
    ha(i).MarkerFaceColor = 'w';
end

legends = cell(nalgs, 1);
for i=1:nalgs
    %legends{i} = sprintf('%s (%.2f)', Algs{i,1}, noise_stability(i));
    legends{i} = sprintf('%s', Algs{i,1});
end
legend(legends, 'Interpreter','latex','Location','southeast');
xlabel('$\sigma$','Interpreter','latex');
ylabel('Err','Interpreter','latex');
ax.FontSize = 10;
ax.TickLabelInterpreter = "latex";
ylim([1e-4 100]);
pdfprint('results/testrd_noise_ord4', hf);
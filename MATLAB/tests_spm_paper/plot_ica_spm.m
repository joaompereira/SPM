clearvars
clc
addpath '../'
addpath '../helper_functions/'
load results\ica_spm.mat

close(figure(2))
hf2 = figure(2);
clf
hf2.Position = [600 250 350 400];
for i=1:ncols*nrows
    ax = subplot(nrows, ncols, i);
    plot(t, SPM_icasig(skew_inds(7*(i-1)+1),:))
    xlim([t(1), t(end)]);
    xlabel('$t$ (seconds)', Interpreter='latex')
    ax.TickLabelInterpreter = "latex";
end

peaki = skew_inds(1);

close(figure(1))
hf1 = figure(1);
hf1.Position = [200 250 350 400];
close(figure(3))
hf3 = figure(3);
hf3.Position = [1000 250 350 400];
clf
for i=1:ncols*nrows
    subplot(nrows, ncols, i)
    plot(t, [data(:, i)'; data(:, i)'-A(i, peaki)*SPM_icasig(peaki,:)]);
    yl=ylim;
    figure(hf1)
    ax = subplot(nrows, ncols, i);
    plot(t, data(:, i))
    xlim([t(1), t(end)]);
    ylim(yl)
    xlabel('$t$ (seconds)', Interpreter='latex')
    ax.TickLabelInterpreter = "latex";
    figure(hf3)
    ax = subplot(nrows, ncols, i);
    plot(t, data(:, i)'-A(i, peaki)*SPM_icasig(peaki,:), ...
        "Color", [0.8500    0.3250    0.0980]);
    xlim([t(1), t(end)]);
    ylim(yl)
    xlabel('$t$ (seconds)', Interpreter='latex')
    ax.TickLabelInterpreter = "latex";
end

pdfprint('results/ica_spm_data', hf1);
pdfprint('results/ica_spm_components', hf2);
pdfprint('results/ica_spm_denoised', hf3);



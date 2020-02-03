clear all
clc
addpath '../'
addpath '../helper_functions/'

% L and R values
Lval = 10:40;
Rval = 8:8:800;
nL = length(Lval);
nR = length(Rval);

% Number of tries for each algorithm
% Some algorithms have more variable running times, so we run them more
% times to produce smoother results
ntries = 10;

% Time of execution and estimate error for different algorithms and (L,R)
% values
freq = zeros(nL,nR);

for i=1:nL
    
    L = Lval(i);
    fprintf('L =%3d\n R =    ',L);
    
    Rmax = L*L/2;
    
    for j=1:nR
      
        R = Rval(j);
        fprintf('\b\b\b%3d',R);   
        if R > Rmax
            freq(i,j) = 0;
            continue
        end
        
        % Pre-allocate local time and error vectors for all tries
        freqij = zeros(ntries,1);
        
        for tries=1:ntries
            
            % True rank decomposition
            A_true = randn(L,R);

            % Generate 4D tensor
            T = generate_lowrank_tensor(A_true, 4);
            
            % Use SPM to get rank decomposition
            A_est = subspace_power_method(T, L, 4);
            
            error = rderror(A_est, A_true, 4);
            
            freqij(tries) = error<1e-4;
            
        end
        
        % Average over all runs
        freq(i,j) = mean(freqij);
    end
    
    fprintf('\n')
end

% Save results
clear T % T is cleared since it occupies a lot of space
save results/testrd_rank.mat

% Plot results
hf = figure(1);
colormap(linspace(0,1,11)'*[1,1,1])
hs = imagesc(Lval,Rval,freq');
xlabel('$L$','Interpreter','latex');
ylabel('$R$','Interpreter','latex');
set(gca,'YDir','normal');
hold on
plot(Lval,Lval.*(Lval-1)/2,'r','LineWidth',1);
plot(Lval,Lval.*(Lval-3)/2,'b--','LineWidth',1);
hl = legend({'$R=\frac{L(L-1)}{2}$','$R=\frac{L(L-3)}{2}$'},...
        'Interpreter','latex','Location','northwest','FontSize',12);

ylim([8,800])
    
NewHeight = hl.Position(4) * 1.5;
hl.Position(2) = hl.Position(2) - (NewHeight - hl.Position(4));
hl.Position(4) = NewHeight;

colorbar
    
pdfprint('results/testrd_rank', hf);
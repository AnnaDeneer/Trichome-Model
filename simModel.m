function [t, y] = simModel(k, figs)
% Script to run the trichome patterning model for a certain parameter set
% Input arguments:
% k:    Parameter set to solve the model for
% figs: Boolean to check whether to show figures or not. 
%       1 = show figures, 0 = no figures
%
% An example wild-type parameter set is supplied in Parameterset_wt.mat
% and the ttg1-9 mutation for this set is given in Parameterset_ttg19.mat

model = @Trichome_eqns;
NVar = 7;
label = {'TTG1','GL1','GL3', 'TRY', 'CPC','AC1', 'AC2'};

% Standard 20x20 grid:
xmax = 20; % number of cells in x-direction
ymax = 20; % number of cells in y-direction

%% Initialization

% Initialize cell grid indices
ctr = cind(1,1:ymax,1:xmax,ymax,NVar);

% Construct spatial coupling matrix for hexagonal cells and periodic
% boundary conditions, last arg: 1 = zero flux, 0 = periodic
D = diffusionMatrix(ymax,xmax,[-1 1 0 0 1 -1],[0 0 -1 1 -1 1],1);

% Set integration options including structure of the jacobian
options = odeset('Vectorized','on','JPattern',jMatrix(NVar,D));
% Set timespan, leave sufficient time for solver to reach steady state
tspan = logspace(0, 4, 1000);

% Compute single cell steady state: start from zero
% protein levels and leave sufficient time to reach steady state
% Use 1 % random perturbations of the homogeneous steady state as
% initial conditions for the grid model
[~,y] = ode15s(model, tspan, zeros(NVar,1), [], 1, 0, k);
ss = y(end,:);
y0 = repmat(ss(:),ymax*xmax,1) .* (1 + 0.01.*rand(NVar*ymax*xmax,1));

%% Integrate the model
tstart = tic;
[t,y] = ode15s(model, tspan, y0(:), options, ctr, D, k);

% Display results statistics if parameter 'figs' is 1
if figs
    % Calculate pattern quantities like trichome density, clusters and
    % coefficient of variation (CV)
    ssAC = y(end, ctr+5) + y(end,ctr+6);  % The steady state of the sum of
                                          % both active complexes
    th  = 0.5*max(ssAC);                  % The threshold to determine cell fate
    nT  = numel(find(ssAC>=th));          % Number of trichomes
    CD  = countClusters(ssAC,th,xmax,ymax,D)/nT; % Cluster density
    rnn = nn(ssAC./max(ssAC),D);          % Nearest neighbour distribution
    CV  = std(rnn)/mean(rnn);             % Coefficient of variation
    TD  = nT/(ymax*xmax);                 % Trichome density
    
    % Print results
    fprintf(1,'\nMaxima/minima and sum/var/mean of expression levels : \n');
    for k=1:NVar
        fprintf(1,'%s : \t min %.2f \t max %.2f \t sum %.2f \t var %.4f \t mean %.4f \n',...
            label{k}, ...
            min(y(end,ctr+k-1)), ...
            max(y(end,ctr+k-1)), ...
            sum(y(end,ctr+k-1)), ...
            var(y(end,ctr+k-1)), ...
            mean(y(end,ctr+k-1)));
    end
    fprintf(1, 'Cluster density: %.2f \t Trichome density: %.2f \t CV: %.2f\n', ...
            CD, TD, CV);
    toc(tstart);
end

%% Plot figures for all species

if figs
    figure();
    for i=1:NVar
        C = y(end,ctr+(i-1));
        C = C./max(C);
        if figs == 1
            subplot(2,4,i)
            plotHexagonals(C,ymax,xmax)
            title(label{i})
        end
    end
    % Finally plot the sum of AC1 and AC2, used as indicator for trichome
    % cells. When AC1+AC2 > 0.5*max(AC1+AC2) a cell is considered to be
    % a trichome.
    subplot(2,4,8)
    AC1 = y(end, ctr+5);
    AC2 = y(end, ctr+6);
    sumAC = AC1+AC2;
    plotHexagonals(sumAC./max(sumAC),ymax,xmax)
    title('AC1 + AC2')
end

end
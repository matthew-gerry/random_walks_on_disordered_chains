% data_chains.m

% Through a range of values of the proportion of A sites, and delta_gamma,
% generate a set of chains. Generate probability distributions as a
% function of time for each chain. This data can be used to calculate the
% scaled cumulants of the distributions given the properties of the chain.
% Multiple realizations are used at each correlation value to capture the
% differences between realizations.

% A matrix file is saved to the parent directory

% Matthew Gerry, May 2024

% Set the bias regime
bias_string = "nobias";

% Current date and time for filenames
time_string = string(datetime("now",Format="yy-MM-dd_HHmm"));
filename = strcat("../RWdata_",bias_string,"_",time_string);

if bias_string=="nobias"
    b = 0.0;
elseif bias_string=="lowbias"
    b = 0.2;
end % condition

% Set time values for dynamics
dt = 25.0;
tmax = 100;
time = 0:dt:tmax;

% Set parameter values required for generating a set of chains
numsites = 361; % Chain length
site_list = -floor(numsites/2):floor(numsites/2); % Site labels
epsilon = 0.05; % Tolerance for NN-correlation values
set_size = 5; % Number of realizations at each set of param values

% Set parameter values required for generating rate matrices
tau = 1.0;
ga_av = 1.0;
site0 = 0; % Initial site of walker

% Parameters to be varied
% Range of p values
dp = 0.1; p_list = 0:0.1:1;

% Range of d_ga values
ddga = 0.2; % Amount by which to step between delta gamma values
dga_list = 0:ddga:2*ga_av;
 
% List of random seeds to cycle through at each c-value
seed_list = 1:set_size;

% Pre-allocate a 3-index tensor whose 2-d slices correspond to different p-values
% Each row in a slice is to be a different chain realization
chains = zeros(length(p_list), set_size, numsites);

% Populate the chains tensor with according binary sequences
for jj=1:length(p_list)
    p = p_list(jj);

    for kk=1:set_size
        seed = seed_list(kk); % Slightly redundant, but allows us to play around with the seed list

        chains(jj, kk, :) = p_chain(p, numsites, seed);

    end % kk
end % jj

% Pre-allocate a 5-index tensor to store time-series data for each chain
% In addition to gaining an index for time, it gains another additional index
% relative to the chains tensor over which we vary the value of dga
dists = zeros(length(dga_list), length(p_list), set_size, numsites, length(time));

for ii=1:length(dga_list)
    % Set gamma values for A and B sites based on dga and ga_av
    dga = dga_list(ii);
    ga_a = ga_av + 0.5*dga; ga_b = ga_av - 0.5*dga;

    for jj=1:length(p_list)
        for kk=1:set_size
            % Call L_chain for each choice of p and each random seed
            L = L_chain(chains(jj, kk, :), b, ga_a, ga_b, tau);
    
            PDF_temp = pdf_L(L, site0, dt, tmax);
            dists(ii, jj, kk, :, :) = PDF_temp;
    
        end % kk
    end % jj
end % ii


% Save parameters, chains, and distributions to file
save(filename,"dt","tmax","numsites","epsilon", "set_size","dp","tau","ga_av","ddga","b","site0","chains","dists")

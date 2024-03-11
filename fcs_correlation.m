% fcs_correlation.m

% Through a range of values of nearest-neighbour correlation between
% site-types, generate a set of chains. Perform full counting statistics on
% each, to establish the values of the cumulants with error due to the 
% finite-sample size and finite chain length

% Matthew Gerry, March 2024

% Set time values for dynamics
dt = 1.0;
tmax = 100;
time = 0:dt:tmax;

% Set parameter values required for generating a set of chains
p = 0.5; % Proportion of A-type sites
b = 0; % Bias
numsites = 241; % Chain length
epsilon = 0.05; % Tolerance for NN-correlation values
set_size = 10; % Number of realizations at each c value

c_list = -1.0:0.1:1.0; % Range of correlation values
seed_list = 1:set_size; % List of random seeds to cycle through at each c-value

% Set parameter values required for generating rate matrices
tau = 1.0;
ga_av = 1.0; dga = 1.0;
ga_a = ga_av + 0.5*dga; ga_b = ga_av - 0.5*dga;


% Pre-allocate a 3-index tensor whose slices correspond to different c-values
% Each row in a slice is to be a different chain realization
chains = zeros(length(c_list), set_size, numsites);

% Populate the chains tensor with according binary sequences
for ii=1:length(c_list)
    c = c_list(ii);

    for jj=1:set_size
        seed = seed_list(jj); % Slightly redundant, but allows us to play around with the seed list

        chains(ii, jj, :) = corr_chain(p, numsites, c, epsilon, seed);

    end % jj
end % ii


% Pre-allocate a 4-index tensor to store time-series data for each chain
dists = zeros(length(c_list), set_size, numsites, length(time));

for ii=1:length(c_list)
    for jj=1:set_size
        L = L_chain(chains(ii, jj, :), b, ga_a, ga_b, tau);

        [PDF_temp, ~, ~, ~, ~, ~] = pdf_L(L, dt, tmax);
        dists(ii, jj, :, :) = PDF_temp;


    end % jj
end % ii


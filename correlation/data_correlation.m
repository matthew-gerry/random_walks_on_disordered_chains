% data_correlation.m

% Through a range of values of nearest-neighbour correlation between
% site-types, generate a set of chains. Generate probability distributions
% as a function of time for each chain. This data can be used to make
% statistical calculations, determining quantities such as the diffusion
% coefficient given the nearest neghbour correlations. Multiple
% realizations are used at each correlation value to capture the
% differences between realizations.

% Matthew Gerry, March 2024

% Current date and time for filenames
time_string = string(datetime("now",Format="yy-MM-dd_HHmm"));

% Set time values for dynamics
dt = 80.0;
tmax = 4000;
time = 0:dt:tmax;

% Set parameter values required for generating a set of chains
numsites = 961; % Chain length
p = 0.5; % Proportion of A-type sites
site_list = -floor(numsites/2):floor(numsites/2);
epsilon = 0.05; % Tolerance for NN-correlation values
set_size = 10; % Number of realizations at each c value

cmin = -1.0; cmax = 1.0; dc = 0.1;
c_list = cmin:dc:cmax; % Range of correlation values
seed_list = 1:set_size; % List of random seeds to cycle through at each c-value

% Set parameter values required for generating rate matrices
tau = 1.0;
ga_av = 1.0; dga = 1.0;
ga_a = ga_av + 0.5*dga; ga_b = ga_av - 0.5*dga;
b = 0.2; % Bias
site0 = -100; % Initial site of walker

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

        PDF_temp = pdf_L(L, site0, dt, tmax);
        dists(ii, jj, :, :) = PDF_temp;

    end % jj
end % ii

% Save parameters, chains, and distributions to file
save(strcat("../RWdata_lowbias_",time_string),"dt","tmax","numsites","p","epsilon", "set_size","cmin","cmax","dc","tau","ga_av","dga","b","site0","chains","dists")


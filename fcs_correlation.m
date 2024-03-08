% fcs_correlation.m

% Through a range of values of nearest-neighbour correlation between
% site-types, generate a set of chains. Perform full counting statistics on
% each, to establish the values of the cumulants with error due to the 
% finite-sample size and finite chain length

% Matthew Gerry, March 2024

% Set parameter values
p = 0.5; % Proportion of A-type sites
numsites = 121; % Chain length
epsilon = 0.01; % Tolerance for NN-correlation values
set_size = 10; % Number of realizations at each c value

c_list = -0.9:0.1:0.9; % Range of correlation values
seed_list = 1:set_size; % List of random seeds to cycle through at each c-value

% Pre-allocate a tensor whose slices correspond to different c-values
% Each row in a slice is to be a different chain realization
chains = zeros(length(c_list), set_size, numsites);
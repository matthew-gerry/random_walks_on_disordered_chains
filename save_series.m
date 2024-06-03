% Current date and time for filenames
time_string = string(datetime("now",Format="yy-MM-dd_HHmm"));

% Set time values for dynamics
dt = 0.2;
tmax = 240;
time = 0:dt:tmax;

% Set parameter values required for generating a set of chains
numsites = 161; % Chain length
p_list = [0.4,0.5,0.6,0.7];
site_list = -floor(numsites/2):floor(numsites/2);
epsilon = 0.05; % Tolerance for NN-correlation values
set_size = 10; % Number of realizations at each c value

cmin = -1.0; cmax = 1.0; dc = 0.1;
c_list = cmin:dc:cmax; % Range of correlation values
seed_list = 1:set_size; % List of random seeds to cycle through at each c-value

% Set parameter values required for generating rate matrices
tau = 1.0;
ga_av = 1.0; dga = 0.5;
b = 0.2; % Bias
site0 = -100; % Initial site of walker
% Each row in a slice is to be a different chain realization
chains = zeros(length(p_list), set_size, length(time));

% Populate the chains tensor with according binary sequences
for ii=1:length(p_list)
    p = p_list(ii);
    
    for jj=1:set_size
        seed = seed_list(jj); % Slightly redundant, but allows us to play around with the seed list
        rng(seed)
        ga_a = ga_av + 0.5*dga; ga_b = ga_av - 0.5*dga;
        [~, ~, ~, ~, D_av,~,~,~] = pdf_rand(p,b,ga_a,ga_b,tau,numsites,dt,tmax);

        chains(ii, jj, :) = D_av;

    end % jj
end % ii

% Pre-allocate a 4-index tensor to store time-series data for each chain
%dists = zeros(length(c_list), set_size, numsites, length(time));

%for ii=1:length(c_list)
%    for jj=1:set_size
%       L = L_chain(chains(ii, jj, :), b, ga_a, ga_b, tau);
%
 %       PDF_temp = pdf_L(L, site0, dt, tmax);
  %      dists(ii, jj, :, :) = PDF_temp;
%



 %   end % jj
%end % ii
save("rw_dga_biased.mat","dt","tmax","numsites","p_list", "set_size","tau","ga_av","dga","b","chains")


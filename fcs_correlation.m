% fcs_correlation.m

% Through a range of values of nearest-neighbour correlation between
% site-types, generate a set of chains. Perform full counting statistics on
% each, to establish the values of the cumulants with error due to the 
% finite-sample size and finite chain length

% Matthew Gerry, March 2024

% Set time values for dynamics
dt = 40.0;
tmax = 3200;
time = 0:dt:tmax;

% Set parameter values required for generating a set of chains
p = 0.5; % Proportion of A-type sites
b = 0.0; % Bias
numsites = 721; % Chain length
site_list = -floor(numsites/2):floor(numsites/2);
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
% including probability distributions and their time derivatives
dists = zeros(length(c_list), set_size, numsites, length(time));
dpdt = zeros(length(c_list), set_size, numsites, length(time));

for ii=1:length(c_list)
    for jj=1:set_size
        L = L_chain(chains(ii, jj, :), b, ga_a, ga_b, tau);

        PDF_temp = pdf_L(L, dt, tmax);
        dists(ii, jj, :, :) = PDF_temp;
        dpdt(ii, jj, :, :) = L*PDF_temp;

    end % jj
end % ii


% Calculate statistics
site_tensor = repmat(reshape(site_list,[1,1,numsites,1]),[length(c_list),set_size,1,length(time)]);

% Mean
n_av = sum(dists.*site_tensor,3);
v_av = sum(dpdt.*site_tensor,3);

% Diffusion coefficient
D = 0.5*sum(dpdt.*(site_tensor.^2),3) - n_av.*v_av;


% Plot diffusion coefficient as a function of time for many realizations of the same parameters
% At a few choices of correlation value

figure;
c_indices = [2,11,20];

for ii=1:length(c_indices)
    c = c_list(c_indices(ii));
    subplot(1, 3,ii); hold on; box on
    for jj=1:set_size
        plot(time, reshape(D(c_indices(ii),jj,1,:),[1,length(time)]))
    end
    xlim([0,time(end)])
    xlabel("$t$",interpreter="latex")
    ylabel("$D$",interpreter="latex")
    title(strcat("$c=",num2str(c),"$"),interpreter="latex")
    set(gca, fontsize=14)
    hold off
end

%% Animate probability distribution over time
figure;
for kk=1:length(time)
    bar(site_list, reshape(dists(20,3,:,kk),[1,numsites]))
    ylim([0,0.1])
    drawnow
    pause(0.02)
end



% animate.m

% Animates the time evolution of the probability distribution of a random
% walk on a disordered chain. Takes data from a .mat file generated using
% data_p_dga.m.

% Matthew Gerry, June 2024

load("../RWdata_nobias_24-06-06_1500.mat") % Load the data in 

% Re-define some of the parameter-dependent variables
time = 0:dt:tmax; % Time
site_list = -floor(numsites/2):floor(numsites/2); % List of site indices

% Uncomment these if we want to incorporate annotations to the plot that
% show the values of p, dga, seed
% p_list = 0:dp:1; % List of p values
% dga_list = 0:ddga:2*ga_av; % List of dga values
% seed_list = 1:set_size; % List of random seeds to cycle through at each c-value

% Animate the evolution of the probability distribution
figure;
for kk=1:length(time)
    bar(site_list, reshape(dists(6,6,1,:,kk),[1,numsites]))
    ylim([0,0.04])
    drawnow
    pause(0.05)
end
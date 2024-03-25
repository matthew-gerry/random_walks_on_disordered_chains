% fcs_correlation.m

% Carry out full counting statistics for random walks on disordered chains
% with varying nearest-neighbour correlation values. Import data generated
% by the function data_correlation.m which calculates probability
% distributions as a function of time for a range of nearest-neighbour
% correlation values.

% Matthew Gerry, March 2024

load("../RWdata_lowbias_24-03-25_1122.mat") % Load the data in 

% Re-define some of the parameter-dependent variables
time = 0:dt:tmax; % Time
site_list = -floor(numsites/2):floor(numsites/2); % List of site indices
c_list = cmin:dc:cmax; % Range of correlation values
seed_list = 1:set_size; % List of random seeds to cycle through at each c-value
ga_a = ga_av + 0.5*dga; ga_b = ga_av - 0.5*dga; % Gamma values

% Calculate time derivative of the probability distribution at each time step
dpdt = zeros(length(c_list), set_size, numsites, length(time));
for ii=1:length(c_list)
    for jj=1:set_size
      
        L = L_chain(chains(ii, jj, :), b, ga_a, ga_b, tau); % Rate matrix for the chain
        dpdt(ii, jj, :, :) = L*reshape(dists(ii, jj, :, :),[numsites,length(time)]);
    
    end % jj
end % ii

%% Calculate statistics
site_tensor = repmat(reshape(site_list,[1,1,numsites,1]),[length(c_list),set_size,1,length(time)]);

% Mean
n_av = sum(dists.*site_tensor,3);
v_av = sum(dpdt.*site_tensor,3);

% Diffusion coefficient
D = 0.5*sum(dpdt.*(site_tensor.^2),3) - n_av.*v_av;


% Plot diffusion coefficient as a function of time for many realizations of the same parameters
% At a few choices of correlation value

figure;
c_indices = [1,2,6,11];

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
    bar(site_list, reshape(dists(3,5,:,kk),[1,numsites]))
    ylim([0,0.04])
    drawnow
    pause(0.05)
end

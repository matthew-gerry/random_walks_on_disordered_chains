% stats_p_dga.m

% Carry out full counting statistics for random walks on disordered chains
% with varying proportion of A sites and varying inhomogeneity. Import data
% generated by the function data_p_dga.m which calculates probability
% distributions as a function of time for a range of these parameter values

% Matthew Gerry, May 2024

load("../RWdata_nobias_24-06-06_1201.mat") % Load the data in 

% Re-define some of the parameter-dependent variables
time = 0:dt:tmax; % Time
site_list = -floor(numsites/2):floor(numsites/2); % List of site indices
p_list = 0:dp:1; % List of p values
dga_list = 0:ddga:2*ga_av; % List of dga values
seed_list = 1:set_size; % List of random seeds to cycle through at each c-value


% Calculate time derivative of the probability distribution at each time step
dpdt = zeros(length(dga_list), length(p_list), set_size, numsites, length(time));

for ii=1:length(dga_list)
    % Recreate gamma values for A and B sites based on dga and ga_av
    dga = dga_list(ii);
    ga_a = ga_av + 0.5*dga; ga_b = ga_av - 0.5*dga;
    for jj=1:length(p_list)
        for kk=1:set_size
          
            L = L_chain(chains(jj, kk, :), b, ga_a, ga_b, tau); % Rate matrix for the chain
            dpdt(ii, jj, kk, :, :) = L*reshape(dists(ii, jj, kk, :, :),[numsites,length(time)]);
        
        end % kk
    end % jj
end % ii

%% Calculate statistics
site_tensor = repmat(reshape(site_list,[1,1,1,numsites,1]),[length(dga_list), length(p_list),set_size,1,length(time)]);

% Mean
n_av = sum(dists.*site_tensor,4);
C1 = sum(dpdt.*site_tensor,4);

% Diffusion coefficient
S = sum(dists.*(site_tensor.^2),4) - n_av.^2;
C2 = sum(dpdt.*(site_tensor.^2),4) - 2*n_av.*C1;

% Third and fourth cumulants
% Create repititions of the above results with the appropriate shapes for ease of calculations
n_av_tensor = repmat(n_av,[1, 1, 1, numsites, 1]);
C1_tensor = repmat(C1,[1, 1, 1, numsites, 1]);
C2_tensor = repmat(C2,[1 ,1, 1, numsites, 1]);

% Skewness
C3 = sum(dpdt.*(site_tensor-n_av_tensor).^3, 4) - 3*sum(dists.*C1_tensor.*(site_tensor-n_av_tensor).^2, 4);

% Kurtosis
C4 = sum(dpdt.*(site_tensor-n_av_tensor).^4, 4) - 4*sum(dists.*C1_tensor.*(site_tensor-n_av_tensor).^3, 4) - 6*sum(dists.*C2_tensor.*(site_tensor-n_av_tensor).^2, 4);


% Plot C2 as a function of time for many realizations of the same parameters
% At a few choices of p value, specific choice of dga

figure;
p_indices = [2,6,10];

for ii=1:length(p_indices)
    p = p_list(p_indices(ii));
    subplot(1, length(p_indices), ii); hold on; box on
    for jj=1:set_size
        plot(time, reshape(C2(6, p_indices(ii),jj,1,:),[1,length(time)]))
    end
    xlim([0,time(end)])
    xlabel("$t$",interpreter="latex")
%     ylabel("$D$",interpreter="latex")
    title(strcat("$p=",num2str(p),"$"),interpreter="latex")
    set(gca, fontsize=14)
    hold off
end

%% Animate probability distribution over time
figure;
for kk=1:length(time)
    bar(site_list, reshape(dists(6,6,1,:,kk),[1,numsites]))
    ylim([0,0.04])
    drawnow
    pause(0.05)
end
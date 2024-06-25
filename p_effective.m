% p_effective.m

% Testing the idea that cumulant values are instantaneously determined by
% an "effective" value of p, based on the segment of the chain the
% population distribution is sampling at any given moment.

% Matthew Gerry, June 2024


% Specify name of data file, calculate C2
filename = "RWdata_nobias_24-06-07_1547";
% filename = "RWdata_lowbias_24-06-07_1147";
[C1,C2,~,~] = stats_p_dga(filename);

% Load parameter values from data file
load(strcat("../", filename, ".mat"), "dt","tmax","numsites","set_size","dp","tau","ga_av","ddga","b","site0","chains","dists")
%%
% Re-define some of the parameter-dependent variables
time = 0:dt:tmax; % Time
site_list = -floor(numsites/2):floor(numsites/2); % List of site indices
dga_list = 0:ddga:2*ga_av; % List of dga values
p_list = 0:dp:1; % List of p values

% As a test, approximate p_effective as a function of time for one point in
% parameter space

dga_index = 6; p_index = 6;
dga = dga_list(dga_index); p = p_list(p_index);

figure; hold on; box on
for jj=1:1
    dist = transpose(reshape(dists(dga_index,p_index,jj,:,:),[numsites,length(time)]));
    chain = reshape(chains(p_index,jj,:),[numsites,1]);

    p_eff = dist*chain;

    C2_eff = (1+exp(-b))./(p_eff*(ga_av + 0.5*dga) + (1-p_eff)*(ga_av - 0.5*dga));

    plot(time, C2_eff)
    plot(time, reshape(C2(dga_index, p_index, jj, :),[1,length(time)]), 'x')
end
hold off
% chain_analysis.m

% Here we analyze the specific structure of chains that give rise to
% outlier behaviour when it comes to the steady-state values of the
% cumulants. Focus on chains 12 and 18 to start as they exhibit low and
% high values of C2 with no bias, respectively. This whole script is pretty
% rough, and created for exploratory purposes.

% Matthew Gerry, June 2024

% Specify name of data file, calculate cumulants
filename = "RWdata_nobias_24-06-07_1547";

% Load chain sequences in from data file and other useful information
load(strcat("../", filename, ".mat"), "dt", "tmax", "numsites", "set_size", "site0", "chains")
time = 0:dt:tmax;
site_list = -floor(numsites/2):floor(numsites/2); % List of site indices

chain12 = reshape(chains(6, 12, :), [1, numsites]);
chain18 = reshape(chains(6, 18, :), [1, numsites]);

% "integrals" of the chains to see if there are variations in the local
% structure
integral12 = zeros([1, numsites]);
integral18 = zeros([1, numsites]);
for ii=1:numsites
    integral12(ii) = sum(chain12(1:ii));
    integral18(ii) = sum(chain18(1:ii));
end

% Correlation analysis
bin12 = 2*chain12 - 1; bin18 = 2*chain18 - 1;
corr12 = mean(bin12(1:end-1).*bin12(2:end));
corr18 = mean(bin18(1:end-1).*bin18(2:end));


%% Testing the idea that the area around the starting point has an outsized influence

% Create a cooked up chain with a extreme local structure around the centre
% and otherwise random site types
window = 20; % Size of a window around the origin with homogeneous structure

test_chains = zeros(2, numsites);
% test_chains(:, 2:2:numsites) = 1;
% test_chains(:, randperm(numsites)) = test_chains;
% test_chains(1, floor(numsites/2)-window/2:ceil(numsites/2)+window/2) = 0;
% test_chains(2, floor(numsites/2)-window/2:ceil(numsites/2)+window/2) = 1;

test_chains(1, :) = [chain12(1:floor(numsites/2)), chain18(ceil(numsites/2):end)];
test_chains(2, :) = [chain18(1:floor(numsites/2)), chain12(ceil(numsites/2):end)];

test_dists = zeros(2, numsites, length(time));
test_dpdt = zeros(2, numsites, length(time));

% Set gamma values for A and B sites to central values for testing purposes
ga_a = 1.5; ga_b = 0.5;
tau = 1; % standard value
b = 0; % no bias

for kk=1:2
        % Call L_chain for each choice of p and each random seed
        L = L_chain(test_chains(kk, :), b, ga_a, ga_b, tau);

        PDF_temp = pdf_L(L, site0, dt, tmax);
        test_dists(kk, :, :) = PDF_temp;
        test_dpdt(kk, :, :) = L*PDF_temp;
end % kk

%% Animate the evolution of the two and compare
figure;
for kk=1:length(time)
    bar(site_list, reshape(test_dists(1,:,kk),[1,numsites]), facealpha=0.5)
    hold on
    bar(site_list, reshape(test_dists(2,:,kk),[1,numsites]), facealpha=0.5)
    hold off
    ylim([0,0.04])
    drawnow
    pause(0.05)
end

%% Calculate C2 for the two and compare

site_tensor = repmat(reshape(site_list,[1,numsites,1]),[2,1,length(time)]);

n_av = sum(test_dists.*site_tensor, 2);
C1 = sum(test_dpdt.*site_tensor, 2);

C2 = sum(test_dpdt.*(site_tensor.^2), 2) - 2*C1.*n_av;
    
figure(); hold on
plot(time, C2(1,:))
plot(time, C2(2,:))
hold off

%% Let's have a look at just the early evolution of the probability distributions

dt = 1; tmax = 500; time = 0:dt:tmax;

% Set gamma values for A and B sites to central values for testing purposes
ga_a = 1.5; ga_b = 0.5;
tau = 1; % standard value
b = 0; % no bias

short_time_chains = zeros(2, numsites);
short_time_chains(1,:) = chain12; short_time_chains(2,:) = chain18;

short_time_dists = zeros(2, numsites, length(time));
short_time_dpdt = zeros(2, numsites, length(time));

for kk=1:2
        % Call L_chain for each choice of p and each random seed
        L = L_chain(short_time_chains(kk, :), b, ga_a, ga_b, tau);

        PDF_temp = pdf_L(L, site0, dt, tmax);
        short_time_dists(kk, :, :) = PDF_temp;
        short_time_dpdt(kk, :, :) = L*PDF_temp;
end % kk

%% Animate the evolution of the two and compare
figure;
for kk=1:length(time)
    bar(site_list, reshape(short_time_dists(1,:,kk),[1,numsites]), facealpha=0.5)
    hold on
    bar(site_list, reshape(short_time_dists(2,:,kk),[1,numsites]), facealpha=0.5)
    hold off
    ylim([0,0.04])
    drawnow
    pause(0.02)
end


%% Calculate C2 for the two and compare

site_tensor = repmat(reshape(site_list,[1,numsites,1]),[2,1,length(time)]);

n_av = sum(short_time_dists.*site_tensor, 2);
C1 = sum(short_time_dpdt.*site_tensor, 2);

C2 = sum(short_time_dpdt.*(site_tensor.^2), 2) - 2*C1.*n_av;
    
figure(); hold on
plot(time, C2(1,:))
plot(time, C2(2,:))
hold off

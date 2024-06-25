% plot_C2_steady.m

% Plot the steady state value of the second scaled cumulant as a function
% of parameters (p, dga). Average over the different realizations and show
% error bars representing the standard deviation. Grab data from a .mat
% file in the parent directory.

% Matthew Gerry, June 2024

% Specify name of data file, calculate C2
% filename = "RWdata_nobias_24-06-07_1547";
filename = "RWdata_nobias_24-06-18_1524";

[~,C2,~,~] = stats_p_dga(filename);

% Load parameter values from data file
load(strcat("../", filename, ".mat"), "dt","tmax","numsites","set_size","dp","tau","ga_av","ddga","b","site0")

% Re-define some of the parameter-dependent variables
time = 0:dt:tmax; % Time
site_list = -floor(numsites/2):floor(numsites/2); % List of site indices
dga_list = 0:ddga:2*ga_av; % List of dga values
p_list = 0:dp:1; % List of p values

% Calculate theoretical value of C2 for comparison
[P, DGA] = meshgrid(p_list, dga_list);
ga_star = P.*(ga_av + 0.5*DGA) + (1 - P).*(ga_av - 0.5*DGA);
C2_star = 2*tau^2./ga_star;

epsilon = 0.05; % Tolerance in steady-state values
window_size = 10; % Number of time steps considered minimum for steady state

% Get steady state values for each time series, reducing over last index
C2_ss = zeros([length(dga_list), length(p_list), set_size]);
for ii=1:length(dga_list)
    for jj=1:length(p_list)
        for kk=1:set_size
            [ss_val, ~] = steady_state(C2(ii, jj, kk, :), epsilon, window_size);
            C2_ss(ii, jj, kk) = ss_val;

        end % kk
    end % jj
end % ii

% Slice off last row since steady states mostly not found for these
C2_ss = C2_ss(1:end-1, :, :);
C2_star = C2_star(1:end-1, :, :);

% Take mean and standard devation over different chain realizations at each
% set of parameter values
C2_ss_mean = mean(C2_ss, 3, "omitnan");
C2_ss_stddev = std(C2_ss, 0, 3, "omitnan");

% Calculate the standard error in the estimate
% For certain realizations, steady state was reached at fewer p, dga pairs
% We deal with this by excluding them from the averages, and calculating
% the standard error just based on the sample of values we got
C2_ss_samplesize = sum(~isnan(C2_ss), 3);
C2_ss_stderr = C2_ss_stddev./sqrt(C2_ss_samplesize);

%%
% Plot C2 at steady state as a function of p for a few choices of dga
dga_indices = [2, 6, 9];

figure(1)
for ii=1:length(dga_indices)
    subplot(length(dga_indices), 1, ii); hold on
    box on

    dga_index = dga_indices(ii);

    errorbar(p_list, C2_ss_mean(dga_index, :), C2_ss_stddev(dga_index, :), '.')
    plot(p_list, C2_star(dga_index, :), '--k')

    title(strcat("$\Delta\gamma=",num2str(dga_list(dga_index)),"$"), Interpreter="latex")
    xlabel("$p$", Interpreter="latex")
    ylabel("$\mathcal{C}_2$", Interpreter="latex")

    set(gca, fontsize=14)
    hold off
end % ii


% Now plot at steady state as a function of dga for a few choices of p
p_indices = [2, 6, 9];

figure(2)
for jj=1:length(p_indices)
    subplot(length(p_indices), 1, jj); hold on
    box on

    p_index = p_indices(jj);

    errorbar(dga_list(1:end-1), C2_ss_mean(:, p_index), C2_ss_stddev(:, p_index), '.')
    plot(dga_list(1:end-1), C2_star(:, p_index), '--k')

    title(strcat("$p=",num2str(p_list(p_index)),"$"), Interpreter="latex")
    xlabel("$\Delta\gamma$", Interpreter="latex")
    ylabel("$\mathcal{C}_2$", Interpreter="latex")

    set(gca, fontsize=14)
    hold off
end % jj

% stddev_disodered.m

% Probing specifically the relationship between the standard deviation in a
% cumulant (between different realizations with the same parameters), and
% the degree of inhomogeneity (i.e. delta gamma).

% Matthew Gerry, June 2024

% Specify name of data file, calculate C2
filename = "RWdata_nobias_24-06-07_1547";
[~,C2,~,~] = stats_p_dga(filename);

% Load parameter values from data file
load(strcat("../", filename, ".mat"), "dt","tmax","numsites","set_size","dp","tau","ga_av","ddga","b","site0")

% Re-define some of the parameter-dependent variables
time = 0:dt:tmax; % Time
site_list = -floor(numsites/2):floor(numsites/2); % List of site indices
dga_list = 0:ddga:2*ga_av; % List of dga values
p_list = 0:dp:1; % List of p values

epsilon = 0.02; % Tolerance in steady-state values
window_size = 4; % Number of time steps considered minimum for steady state

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

% Calculate the standard deviation over different chain realizations
C2_ss_stddev = std(C2_ss, 0, 3, "omitnan");

% Standard error is the maximum of the standard error and the tolerance in
% estimating the SS value (though arguably it really is zero for dga=0,
% p=0, p=1)
C2_ss_stderr = C2_ss_stddev; C2_ss_stderr(C2_ss_stderr<epsilon) = epsilon;


% Plot the standard deviation as a functon of dga at each value of p
figure; hold on; box on
for ii=1:length(p_list)
    plot(dga_list(1:end-1), C2_ss_stddev(:,ii), DisplayName=strcat("$p =$ ",num2str(p_list(ii))))
end
xlim([0,max(dga_list(1:end-1))])
xlabel("$\Delta\gamma$", Interpreter="latex")
ylabel("stddev($\mathcal{C}_2$)", Interpreter="latex")
legend(location="northwest", Interpreter="latex")
set(gca, fontsize=14)
hold off





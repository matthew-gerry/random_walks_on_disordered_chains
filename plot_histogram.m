% plot_histogram.m

% Plot a few selected "snapshots" of the population distribution over a
% disordered chain as it evolves in time. Reads in data from a .mat file
% generated using the script data_p_dga.m

% Matthew Gerry, June 2024

% filename = "RWdata_lowbias_24-06-07_1147";
filename = "RWdata_nobias_24-06-07_1547";
%%
load(strcat("../", filename, ".mat")) % Load the data in 

% Re-define some of the parameter-dependent variables
time = 0:dt:tmax; % Time
site_list = -floor(numsites/2):floor(numsites/2); % List of site indices
dga_list = 0:ddga:2*ga_av; % List of dga values

% Designate some kwarg values for plotting
colourlist = ["#0072BD", "#D95319", "#77AC30","#7E2F8E"];
lettlist = ["(a) ", "(b) ", "(c) ", "(d) "];
ls_list = ["-","--",":","-."];


dga_indices = [3, 6, 10];
snapshot_times = [10, 35, 90];  % Specific indices of time array at which to show PDF (arb.)



figure(1)

for ii=1:length(dga_indices)
    dga = dga_list(dga_indices(ii));
    subplot(1,3,ii); hold on; box on

    for jj=1:length(snapshot_times)
        t_snap = snapshot_times(jj)*dt;
        hom_curve = plot(site_list, reshape(dists(1,6,1,:,snapshot_times(jj)),[1,numsites]), '-k', linewidth=0.2);
        bar(site_list, reshape(dists(dga_indices(ii),6,1,:,snapshot_times(jj)),[1,numsites]), facecolor=colourlist(jj), facealpha=0.6, EdgeColor="none", DisplayName=strcat("$t=\;$",num2str(t_snap)))

        hom_curve.Annotation.LegendInformation.IconDisplayStyle = "off";

    end % jj
    ylim([0, 0.045])
    if ii==1; ylabel("$P_n$", Interpreter="latex"); legend(Interpreter="latex", Location="northeast"); end
    xlabel("$n$", Interpreter="latex")
    
    yl = ylim; xl = xlim;
    text(0.9*xl(1) + 0.1*xl(2), 0.1*yl(1) + 0.9*yl(2), strcat(lettlist(ii),"$\Delta\gamma=\;$",num2str(dga_list(dga_indices(ii)))), FontSize=14, Interpreter="latex")
    set(gca, fontsize=14)
    

    hold off;
end % ii
%Plotting PDF for different time snap shot

[PDF,sites ,~ ,~,~,~] = pdf_rand(0.5,0,1.5,0.5,1,401,0.2,580);
[PDF_hom, ~, ~, ~,~,~] = pdf_rand(1,0,1.5*0.5+0.5*0.5,1.5*0.5+0.5*0.5,1,401,0.2,580);
rng(10);
snapshot_times = [300,500,1000,1700];
colourlist = ["#0072BD", "#D95319", "#77AC30","#7E2F8E","#EDB120","#4DBEEE"];
    for jj=1:length(snapshot_times)
        t_snap = snapshot_times(jj);
        hold on
        hom_curve = plot(sites, PDF_hom(:,t_snap), '-k', linewidth=0.2);
        bar(sites, PDF(:,t_snap), facecolor=colourlist(jj), facealpha=0.6, EdgeColor="none",DisplayName=strcat("$t=\;$",num2str(0.2*t_snap)));
        hom_curve.Annotation.LegendInformation.IconDisplayStyle = "off";
    end % jj
%     xlim([-10,sites(end)]) % If high bias
    ylim([0, 1.2*max(PDF(:,snapshot_times(1)))])
    yl = ylim; xl = xlim;
    
  
    set(gca, fontsize=14)
    hold off


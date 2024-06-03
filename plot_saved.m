load('rw_vav.mat')
% Load data and plot the values
dt = 0.2;
time = 0:dt:tmax;
ga_a = ga_av + 0.5*dga; ga_b = ga_av - 0.5*dga;
% Loop through each p value

for ii = 1:length(b_list)
    figure; % Create a new figure for each p value
    hold on; % Hold on to plot multiple lines on the same figure
    [~, ~, ~, ~, D_hom,~,~,~] = pdf_rand(1,b_list(ii),ga_a*0.5+ga_b*(1-0.5),ga_a*0.5+ga_b*(1-0.5),tau,numsites,dt,tmax);
    plot(time, D_hom(:), 'r', 'DisplayName', 'Homogeneous Chain')
    title(sprintf('Chains for propability = %g', b_list(ii))); 
    xlabel('Time Step'); 
    ylabel('Chain Value'); 
end
hold off;

for ii = 1:length(b_list)
    figure; % Create a new figure for each p value
    hold on; % Hold on to plot multiple lines on the same figure
    title(sprintf('Chains for propability = %g', b_list(ii))); 
    xlabel('Time Step'); 
    ylabel('Chain Value'); 
    
    % Loop through each chain for the current p value
    for jj = 1:10
        % Use squeeze to reduce the dimensions to 2D
        plot(squeeze(chains(ii, jj, :)));
    end
    
    hold off;
    legend('Chain 1', 'Chain 2', 'Chain 3', 'Chain 4', 'Chain 5', 'Chain 6', 'Chain 7', 'Chain 8', 'Chain 9', 'Chain 10', 'Location', 'bestoutside');
end
 
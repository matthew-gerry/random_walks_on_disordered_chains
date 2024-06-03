
% Generate homogeneous chains and the biased random walk
rng(1145);
%Set seed = 1145
[~ ,~ ,~ ,v_av,~,~] = pdf_rand(0.6,0.2,1.5,0.5,1,161,0.2,120);
[~, ~, ~, v_av_hom,~,~] = pdf_rand(1,0.2,1.5*0.4+0.5*0.6,1.5*0.4+0.5*0.6,1,161,0.2,120);
rand_chain = v_av; 
hom_chain = v_av_hom; 
time = linspace(0, 120, 601);

[steady_value, index] = steady_state(rand_chain, 0.0001,5);

disp(strcat("The steady-state value is",num2str(steady_value)));

convergenceTime = time(index);

% Create a figure
figure;

% Plot the Disordered Chain
plot(time, rand_chain(:), 'b', 'DisplayName', 'Disordered Chain');
hold on; % Hold the current plot

% Plot the Homogeneous Chain
plot(time, hom_chain(:), 'r', 'DisplayName', 'Homogeneous Chain');

if ~isnan(convergenceTime)
    hold on; 
    yLimits = get(gca, 'ylim');
    line([convergenceTime convergenceTime], yLimits, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Convergence Point');
    hold off;
end

xlabel('Time');
ylabel('Average Speed');
title('Average Speed of Homogeneous and Disordered Chain');
legend('Location', 'best');
grid on;

% Release the hold
%Statistical test for average of v_av of different disordered chains
%Different probabilites for binormd
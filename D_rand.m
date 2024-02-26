
% Generate homogeneous chains and the biased random walk
%Set seed = 1145
rng(1145);
[~ ,~ ,~ ,~,D_avg,~] = pdf_rand(0.5,0,1.5,0.5,1,401,0.2,580);
[~, ~, ~, ~,D_avg_hom,~] = pdf_rand(1,0,1,1,1,401,0.2,580);
rand_D = D_avg; 
hom_D = D_avg_hom; 
time = linspace(0, 580, 2901);%longer time

[steady_value, index] = steady_state(rand_D, 0.0001,100);

disp(strcat("The steady-state value is",num2str(steady_value)));

convergenceTime = time(index);

% Create a figure
figure;

% Plot the Disordered Chain
plot(time, rand_D(:), 'b', 'DisplayName', 'Disordered Chain');
hold on; % Hold the current plot

% Plot the Homogeneous Chain
plot(time, hom_D(:), 'r', 'DisplayName', 'Homogeneous Chain');

if ~isnan(convergenceTime)
    hold on; 
    yLimits = get(gca, 'ylim');
    line([convergenceTime convergenceTime], yLimits, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Convergence Point');
    hold off;
end

xlabel('Time');
ylabel('Diffusion Coefficient');
title('Diffusion Coefficient of Homogeneous and Disordered Chain');
legend('Location', 'best');
grid on;

% Release the hold
hold off;

%Statistical test for average of v_av of different disordered chains
%Different probabilites for binormd
%Even at bias = 0, there seems existing a threshold that controls
%the behavior of value of D of the random chain. Say around 0.52 for now
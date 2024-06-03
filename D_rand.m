
% Generate homogeneous chains and the biased random walk
%Set seed = 1145
rng(1145);
[~ ,~ ,~ ,~,D_avg,~,~,~] = pdf_rand(0.6,0,1.5,0.5,1,401,0.2,120);
[~ ,~ ,~ ,~,D_avg_hom,~,~,~] = pdf_rand(1,0,1.5*0.4+0.5*0.6,1.5*0.4+0.5*0.6,1,401,0.2,120);
rand_D = D_avg; 
hom_D = D_avg_hom; 
time = 0:0.2:120;%longer time

% Create a figure
figure;

% Plot the Disordered Chain
plot(time, rand_D(:), 'b', 'DisplayName', 'Disordered Chain');
hold on; % Hold the current plot

% Plot the Homogeneous Chain
plot(time, hom_D(:), 'r', 'DisplayName', 'Homogeneous Chain');


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
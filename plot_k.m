% Generate homogeneous chains and the biased random walk
rng(1145);
%Set seed = 1145
[~ ,~ ,~ ,~,~,~,~,K] = pdf_rand(0.5,0,1.5,0.5,1,401,0.2,580);
[~, ~, ~, ~,~,~,K_hom] = pdf_rand(1,0,1.5*0.5+0.5*0.5,1.5*0.5+0.5*0.5,1,401,0.2,580);
time = linspace(0, 580, 2901);
% Plot the Disordered Chain
plot(time, K(:), 'b', 'DisplayName', 'Disordered Chain');
hold on; % Hold the current plot

% Plot the Homogeneous Chain
plot(time, K_hom(:), 'r', 'DisplayName', 'Homogeneous Chain'); 
xlabel('Time');
ylabel('Average Speed');
title('Average Speed of Homogeneous and Disordered Chain');
legend('Location', 'best');
grid on;
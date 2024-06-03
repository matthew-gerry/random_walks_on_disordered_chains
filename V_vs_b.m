% Range of bias values
biasValues = 0:0.1:3; 

convergence_V = zeros(length(biasValues), 1);

% Number of seed
numSeeds = 10;

% Loop through the bias values
for i = 1:length(biasValues)
    % Current bias value
    currentBias = biasValues(i);
    
    convergence_val = 0;
    
    % Loop through the seeds
    for seed = 1:numSeeds
        rng(seed);

        [~, ~, ~, v_av, ~, ~] = pdf_rand(0.5, currentBias, 1.5, 0.5, 1, 161, 0.2, 120);
        [steady_value, ~] = steady_state(v_av, 0.001,5);
        convergence_val = convergence_val + steady_value;
    end
    
    convergence_V(i) = convergence_val / numSeeds;
end

v = 2.*exp(-biasValues./2).*sinh(biasValues./2).*(1.5*0.5 + 0.5*0.5);
% Plot theoretical v against bias values
figure;
plot(biasValues, v, '-x', 'DisplayName', 'Theoretical v vs. Bias'); 
hold on;
plot(biasValues, convergence_V, '-o', 'DisplayName', 'V_{avg} vs. Bias'); 

xlabel('Bias');
ylabel('V_{avg}');
title('V_{avg} vs. Bias');
grid on;
legend show; %

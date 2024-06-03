% Range of bias values
biasValues = 0:0.05:0.5; 

convergence_D = zeros(length(biasValues), 1);

% Number of seed
numSeeds = 5;

% Loop through the bias values
for i = 1:length(biasValues)
    % Current bias value
    currentBias = biasValues(i);
    
    convergence_val = 0;
    
    % Loop through the seeds
    for seed = 1:numSeeds
        rng(seed);

        [~, ~, ~, ~, D_avg, ~] = pdf_rand(0.5, currentBias, 1.5, 0.5, 1, 401, 0.2, 580);
        [steady_value, ~] = steady_state(D_avg, 0.01,100);
        convergence_val = convergence_val + steady_value;
    end
    
    convergence_D(i) = convergence_val / numSeeds;
end

% Plot the bias against D_avg 
plot(biasValues, convergence_D, '-o');
xlabel('Bias');
ylabel('Average D_avg');
title('Average D_avg vs. Bias');
grid on;
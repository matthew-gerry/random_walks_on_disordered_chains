function [steady_value,time_index] = steady_state(series, tolerance, window_size)
% Find the steady state value and convergence time of time series
% Compare the difference between some consecutive observations to a tolerance
% value and determine the steady state

% Initialize the outputs
time_index = 0;
steady_value = NaN;

    % Slide window across the time series to find steady state
    for i = 1:(length(series) - window_size + 1)
        window = series(i:i+window_size-1);
        if max(window) - min(window) <= tolerance
            steady_value = mean(window);
            time_index = i;
            break; % Exit loop once a steady state window is found
        end
    end
    
    % If start index remains 0, steady state was not found
    if time_index == 0
        warning('Steady state not found within the given tolerance.');
    end
end
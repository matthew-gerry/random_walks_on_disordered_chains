% A function to generate a binary sequence whose nearest neighbour 
% correlations are within a tolerance of a desired chosen value

% Uses the function corrcoef() from the Statistics and Machine Learning
% toolbox.

% Matthew Gerry, October 2023

%%% Arguments %%%
% p - the proportion of chain sites of the first (A) type
% N - total length of the chain
% c - desired nearest-neighbour correlation value
% epsilon - tolerance for correlation value (will not be possible to get exact)
% seed - random seed to initialize random number generator

function [chain, c_calculated] = corr_chain(p, N, c, epsilon, seed)

    rng(seed) % Initialize the random number generator

    num_zeros = round(p*N);

    % Generate a chain with the appropriate number of zeros and ones
    % Start from the extreme that is closer to the desired correlation
    % value (it is easier to make the chain more disordered as we go along,
    % rather than more ordered--it'll be very difficult to "pass through"
    % the zero point)
    if c>=0
        % Correlated extreme
        chain = [zeros(1, num_zeros), ones(1,N-num_zeros)];
        if abs(c)<eps
            % Randomize the order only if the desired correlation value is
            % very close to zero
            chain = chain(randperm(length(chain)));
        end
    elseif c<0
        % Anticorrelated extreme
        chain = zeros(1,N);
        chain(2:2:N) = 1;
    end % Initialization

    % Calculate the nearest-neighbour correlations
    c_init_matrix = corrcoef(chain(1:N-1),chain(2:N));
    c_current = c_init_matrix(1,2);

    iterations = 0; % Track the number of iteration to limit runtime
    while abs(c_current-c) >= epsilon

        iterations = iterations + 1;
        new_chain = chain;
    
        % Randomly select two positions on the chain
        position1 = randi(N);
        position2 = randi(N);

        % Get a new chain by swapping the values at these positions
        new_chain([position1, position2]) = new_chain([position2, position1]);
        
        % Calculate the correlations with these entries swapped
        c_new_matrix = corrcoef(new_chain(1:N-1), new_chain(2:N));
        c_new = c_new_matrix(1,2);
        
        if abs(c_new - c) < abs(c_current - c)
            % If the new correlation is closer to the desired value,  implement this swap and continue
            chain = new_chain;
            c_current = c_new;
        end

        if iterations >= 100000
            % Break out of the loop if desired correlation isn't achieved
            error("Desired correlation not achieved after 100000 iterations.")            
        end % runtime control

    end % while loop

    % Calculate the nearest-neighbour correlations for the generated chain
    c_calculated = corrcoef(chain(1:N-1),chain(2:N));

end % function
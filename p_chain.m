% p_chain.m

% A function to generate a binary sequence without nearest neighbour
% correlations, with a proportion of 'A'-type sites given by an argument a.

% Matthew Gerry, May, 2024

%%% Arguments %%%
% p       - the proportion of chain sites of the first (A) type
% N       - total length of the chain
% seed    - random seed to initialize random number generator

function chain = p_chain(p, N, seed)
    
    if or(p<0.0, p>1.0)
        error("Choose a value between 0 and 1 for the proportion of A sites.")
    end % error

    rng(seed) % Set the random seed

    NA = floor(p*N); % Number of A-sites
    
    % Initialize a chain with the correct number of zeros and ones
    chain = [zeros([NA,1]);ones([N-NA,1])];
    
    % Shuffle the chain to a random order
    index = randperm(length(chain));
    chain(index) = chain;

end % function
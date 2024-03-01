% L_chain.m

% Function to obtain the rate matrix L given a bit string representing a
% chain. This is to allow for a little more flexibility than the L_rand
% function, which assumes the sites are randomly ordered with no
% correlations.

% Matthew Gerry, March 2024

function L = L_chain(chain,bias,ga_a,ga_b,tau)

    % Pre-allocate square matrix whose size is determined by the chain length
    numsites = length(chain);
    L = zeros(numsites);

    for ii=2:numsites-1 % Leave out edge cases for now
        % Populate the matrix with rates - rate type for the transtion
        % determined by site on the left side of the pair
        
        % Reverse rates into site ii from the right
        L(ii,ii+1) = rates(chain(ii),ga_a,ga_b,tau)*exp(-bias);

        % Forward rates into site ii from the left
        L(ii,ii-1) = rates(chain(ii-1),ga_a,ga_b,tau); 
        
        % Rates out of site ii
        L(ii,ii) = -rates(chain(ii),ga_a,ga_b,tau)...
            - rates(chain(ii-1),ga_a,ga_b,tau)*exp(-bias); 
    end % ii
    
    % Fill in remaining matrix elements
    L(1,2) = rates(chain(1),ga_a,ga_b,tau)*exp(-bias);
    L(numsites, numsites-1) = rates(chain(numsites-1),ga_a,ga_b,tau);

    L(2,1) = rates(chain(1),ga_a,ga_b,tau);
    L(numsites-1,numsites) = rates(chain(numsites-1),ga_a,ga_b,tau)*exp(-bias); 
    
    % Populate the last two diagonals
    L(1,1) = -rates(chain(1),ga_a,ga_b,tau);
    L(numsites, numsites) = -sum(L(:, numsites));

end % function
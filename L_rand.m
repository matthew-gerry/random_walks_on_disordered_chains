function [L, sites, chain] = L_rand(p,bias,ga_a,ga_b,tau,numsites)



    sites = (1:numsites) - 0.5*(numsites+1); 
    % Just an array of the site labels, from -(numsites-1)/2 to (numsites-1)/2
    chain = binorm_d(1, p, numsites) + 1; 
    % Generating a chain of 1s and 2s following a binormial distribution given probability p
    %rng(1145);
    % Set seed=1145 to make sure the results are reproducible and consitant
    % with every trial
 
    
    L = zeros(numsites); % Pre-allocate rate matrix
    
    for ii=2:numsites-1 % Leave out edge cases for now
        % Populate the matrix with rates - rate type for the transtion
        % determined by site on the left side of the pair
        L(ii,ii+1) = rates(chain(ii),ga_a,ga_b,tau)*exp(-bias)...
        ; 
        % Reverse rates into site ii from the right
        L(ii,ii-1) = rates(chain(ii-1),ga_a,ga_b,tau); 
        % Forward rates into site ii from the left
        L(ii,ii) = -rates(chain(ii),ga_a,ga_b,tau)...
            - rates(chain(ii-1),ga_a,ga_b,tau)*exp(-bias); 
        % Rates out of site ii


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
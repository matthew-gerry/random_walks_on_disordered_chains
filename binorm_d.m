function samples = binorm_d(n, p, numsites)
    % n: Number of trials
    % p: Probability of success
    % num_samples: Number of samples
    
    samples = zeros(numsites, 1);
    
    for i = 1:numsites
        successes = 0;
        for j = 1:n
            if rand() < p
                successes = successes + 1;
            end
        end
        samples(i) = successes;
    end
end
%Use this function to generate a random chain follwing binomial 
% distribution, obeying the rule for blocktypes,
%which is 1 for site A and 2 for site B
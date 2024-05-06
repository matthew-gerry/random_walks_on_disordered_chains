% pdf_L

% Function to obtain the probability distribution function for a random
% walk on a chain characterized by rate matrix L, which is taken as an
% argument. This is to allow for more flexibility than the pdf_rand 
% function, which assumes the sites are randomly ordered with no
% correlations. In addition to the probability distribution function, this
% function calculates time series of the first four scaled cumulants
% of each random walk.

% Matthew Gerry, March 2024

% Arguments:
% L     - The rate matrix characterizing the random walk
% site0 - The position of the walker at the start of the simulation
% dt    - Time-step
% tmax  - The length of time for the simulation, which begins at time zero

function PDF = pdf_L(L, site0, dt, tmax)

    % Define a variable for the length of the chain
    numsites = size(L,1);

    % Initial state (1 at central site)
%     p0 = zeros([numsites,1]); p0(floor((numsites+1)/2)) = 1;
    index0 = ceil(numsites/2) + site0;
    p0 = zeros([numsites,1]); p0(index0) = 1;

    % Solve master equation numerically
    time = 0:dt:tmax;

    PDF = zeros(numsites,length(time)); % Pre-allocate time-series of prob dist
    
    % Determine the bias to determine how to best solve the system
    bias = -log(L(2,1)/L(1,2));
%     if bias<0 % At low bias, diagonalize L to calculate PDF faster
%         [V,D] = eig(L);
%         for ii=1:length(time)
%             t = time(ii);
%             PDF(:,ii) = real(V*expm(D*t)*(V\p0)); % Exponential of L*t acting on p0
%         end
%     else % Except don't do this at high bias - leads to strange numerical effects since L is near singular
        for ii=1:length(time)
            t = time(ii);
            PDF(:,ii) = expm(L*t)*p0; % Exponential of L*t acting on p0
        end
%     end

    % Statistics of n - calculate each for the first four cumulants and
%     % scaled cumulants
%     dpdt = L*PDF;
%     numsites = size(L,1);
%     sites = -floor(0.5*numsites):floor(0.5*numsites);
%     
%     % Mean
%     n_av = sum(PDF.*repmat(sites',[1,length(time)]));
%     v_av = sum(dpdt.*repmat(sites',[1,length(time)])); % Mean "velocity"
%     
%     % Diffusion coefficient
%     [n_av_grid, sites_grid] = meshgrid(n_av,sites);
% 
%     S = sum(PDF.*(sites_grid-n_av_grid).^2);
%     D_av = 0.5*(sum(dpdt.*repmat((sites.^2)',[1,length(time)])) - 2*n_av.*v_av);
% 
%     % Skewness   
%     skw = sum(PDF.*(sites_grid-n_av_grid).^3);
%     C3 = skw./time; % Scaled skewness
% 
%     % Kurtosis
%     krt = sum(PDF.*(sites_grid-n_av_grid).^4) - 3*S.^2;
%     C4 = krt./time; % Scaled kurtosis

end % function
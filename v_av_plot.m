
% Generate homogeneous chains and the biased random walk
rng(1145);
%Set seed = 1145
[~ ,~ ,~ ,v_av,~,~] = pdf_rand(0.6,0,1.5,0.5,1,161,0.2,120);
[~, ~, ~, v_av_hom,~,~] = pdf_rand(1,0,1.5*0.6+0.5*0.4,1.5*0.6+0.5*0.4,1,161,0.2,120);
rand_chain = v_av; 
hom_chain = v_av_hom; 
time = linspace(0, 120, 601);

%Analyze convergence
diff = abs(rand_chain-hom_chain);
convergenceThreshold = 0.01;

% Initialize the no_conv flag and a variable for the convergence time
no_conv = true;
convergenceTime = NaN;


for i = 1:10:length(diff)
    if (i+9) <= length(diff)
        subset = diff(i:i+9);
        avg = mean(subset);
        if avg < convergenceThreshold
            disp(['The series steady state value is ', num2str(avg)]);
            convergenceTime = time(i);
            disp(['The series converge at time', num2str(i)]);
            no_conv = false;
            break;
            
        end
    end
end

if no_conv
    disp('The series do not converge');
end

% Create a figure
figure;

% Plot the Disordered Chain
plot(time, rand_chain(:), 'b', 'DisplayName', 'Disordered Chain');
hold on; % Hold the current plot

% Plot the Homogeneous Chain
plot(time, hom_chain(:), 'r', 'DisplayName', 'Homogeneous Chain');

if ~isnan(convergenceTime)
    hold on; 
    yLimits = get(gca, 'ylim');
    line([convergenceTime convergenceTime], yLimits, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Convergence Point');
    hold off;
end

xlabel('Time');
ylabel('Average Speed');
title('Average Speed of Homogeneous and Disordered Chain');
legend('Location', 'best');
grid on;

% Release the hold
%Statistical test for average of v_av of different disordered chains
%Different probabilites for binormd
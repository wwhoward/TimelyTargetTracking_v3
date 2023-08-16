% Written for TTT Journal by W.W.Howard in Spring 2023
% Contact: {wwhoward}@vt.edu Wireless@VT
% For TimelyTrackingNetwork v3
% 
% Task: 
% Compare an AoII where all targets at node are updated vs 

% Header
addpath(genpath(string(pwd)+'\TimelyTrackingNetwork_v3.0'))
clc;clear;close all


% Parameters
nRep = 30; 
% Strats = {"centralized_random", "distributed", "tmp_distributed"}; % Important to do tmp method last as it deletes some tracks
Strats = {"centralized_random", "distributed"}; 
% Strats = {"centralized_random", "distributed"}; 
t_step = 0.2; 
duration = 15; % #seconds to simulate
coverage = 0.1; 
nTrackers = 20; % This is the PPP density per unit km
updateRate = 0.5; % rate per node 

nStrats = length(Strats); 

stats = {}; 
parfor n = 1:nRep
    display("Rep " + string(n))

    myTargets = targetModel(30);

    myTrackers = genTrackers(nTrackers, 'Coverage', coverage, 'updateRate', updateRate);
    
    myRandomFC = fusionCenter(myTargets, myTrackers, Strats{1}, updateRate);
    % myAoiiFC = fusionCenter(myTargets, myTrackers, Strats{2}, updateRate);
    myTmpFC = fusionCenter(myTargets, myTrackers, Strats{2}, updateRate);


    % Run sim
    for t = 0:t_step:duration
        % clc; display("Rep " + string(n) + "/" + string(nRep) + ", " + string(t/duration*100)+"%")
        
        myTargets.update(t_step);
        for i = 1:length(myTrackers) % Since is PPP distributed, don't know how many
            myTrackers{i}.observe(myTargets, t);
        end
        myRandomFC.getUpdates(t);
        % myAoiiFC.getUpdates(t);
        myTmpFC.getUpdates(t);
    end

    % Get stats
    tmp_stats = {}; 
    tmp_stats{1, 1} = myRandomFC.Stats;
    % tmp_stats{2, 1} = myAoiiFC.Stats;
    tmp_stats{2, 1} = myTmpFC.Stats;    

    
    
    for i = 1:nStrats % TODO change to length Strats
        stats{i,n} = tmp_stats{i}; 
    end
end

% myTimelyFC.plotScene()


mean_stats = AverageStats(stats); 

colors = linspecer(length(Strats)); 
% display_names = {"Random", "AoII", "Target Based AoII"}; 
display_names = {"Random", "Target Based AoII"}; 
% display_names = {"AoII", "Centralized", "Random"}; 
% display_names = {"AoI Metric", "Random"}; 

% Ave selected nodes
figure; hold on
for i = 1:length(Strats)
    tmp_mean = mean(mean_stats{"nSelectedNodes"}(i,:)); 
    plot(mean_stats{"TimeSteps"}, mean_stats{"nSelectedNodes"}(i,:), 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', display_names{i}+', mean: '+string(tmp_mean))
end
legend('interpreter', 'latex', 'fontsize', 12)
xlabel('Time (s)', 'interpreter', 'latex', 'fontsize', 12)
ylabel('Nodes', 'interpreter', 'latex', 'fontsize', 12)
title('Number of Selected Nodes', 'interpreter', 'latex', 'fontsize', 16)

% Ave updated targets
figure; hold on
for i = 1:length(Strats)
    tmp_mean = mean(mean_stats{"nSelectedTargets"}(i,:)); 
    plot(mean_stats{"TimeSteps"}, mean_stats{"nSelectedTargets"}(i,:), 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', display_names{i}+', mean: '+string(tmp_mean))
end
legend('interpreter', 'latex', 'fontsize', 12)
xlabel('Time (s)', 'interpreter', 'latex', 'fontsize', 12)
ylabel('Targets', 'interpreter', 'latex', 'fontsize', 12)
title('Number of Updated Targets', 'interpreter', 'latex', 'fontsize', 16)



% Error
figure; hold on 
for i = 1:length(Strats)
    plot(mean_stats{"TimeSteps"}, mean_stats{"Error"}(i,:), 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', display_names{i})
end
legend('Interpreter', 'latex', 'fontsize', 12)
xlabel('Time (s)', 'Interpreter', 'latex', 'fontsize', 12)
ylabel('Error (m)', 'Interpreter', 'latex', 'fontsize', 12)
title('Centralized Selection Error', 'Interpreter', 'latex', 'fontsize', 16)

% More informative error CDF
figure; 
for i = 1:length(Strats)
    semilogx(mean_stats{"ECDF"}{i,2}, mean_stats{"ECDF"}{i,1}, 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', display_names{i})
    hold on
end
legend('Interpreter', 'latex', 'fontsize', 12)
xlabel('Meters', 'Interpreter', 'latex', 'fontsize', 12)
xlim([10^0, 10^3])
ylabel('Pr(Error $\leq X$)', 'Interpreter', 'latex', 'fontsize', 12)
title('Centralized Selection Error Distribution', 'Interpreter', 'latex', 'fontsize', 16)

% Steady state error CDF
figure; 
for i = 1:length(Strats)
    semilogx(mean_stats{"SS_ECDF"}{i,2}, mean_stats{"SS_ECDF"}{i,1}, 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', display_names{i})
    hold on
end
grid on
legend('Interpreter', 'latex', 'fontsize', 12, 'location', 'northwest')
xlabel('Meters', 'Interpreter', 'latex', 'fontsize', 12)
xlim([10^0, 10^3])
ylabel('Pr(Error $\leq X$)', 'Interpreter', 'latex', 'fontsize', 12)
title('Steady State Error Distribution, $\delta=0.5$', 'Interpreter', 'latex', 'fontsize', 16)


% Nodes covered
figure; hold on
plot(mean_stats{"TimeSteps"}, mean_stats{"nTotalTargets"}(1,:), '-k', 'linewidth', 2, 'DisplayName', "Total Targets")
plot(mean_stats{"TimeSteps"}, mean_stats{"nCoveredTargets"}(1,:), '--k', 'linewidth', 2, 'DisplayName', 'Covered Targets')
for i = 1:length(Strats)
    plot(mean_stats{"TimeSteps"}, mean_stats{"nTrackedTargets"}(i,:), '--', 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', "Tracked Targets, "+display_names{i})
end
legend('Interpreter', 'latex', 'fontsize', 12, 'location', 'southeast')
xlabel('Time (s)', 'Interpreter', 'latex', 'fontsize', 12)
ylabel('Number of UAVs', 'Interpreter', 'latex', 'fontsize', 12)
title('Centralized Selection Coverage', 'Interpreter', 'latex', 'fontsize', 16)


% Ages
figure; hold on
for i = 1:length(Strats)
    plot(mean_stats{"TimeSteps"}, mean_stats{"Age"}(i,:), '--', 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', display_names{i})
end
grid on
legend('Interpreter', 'latex', 'fontsize', 12, 'location', 'northeast')
xlabel('Time (s)', 'Interpreter', 'latex', 'fontsize', 12)
ylabel('Age', 'Interpreter', 'latex', 'fontsize', 12)
title('Age of Tracked Targets', 'Interpreter', 'latex', 'fontsize', 16)

% Peak Age
figure; hold on
for i = 1:length(Strats)
    plot(mean_stats{"TimeSteps"}, mean_stats{"PeakAge"}(i,:), '--', 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', display_names{i})
end
legend('Interpreter', 'latex', 'fontsize', 12, 'location', 'northeast')
xlabel('Time (s)', 'Interpreter', 'latex', 'fontsize', 12)
ylabel('Peak Age', 'Interpreter', 'latex', 'fontsize', 12)
title('Peak Age of Tracked Targets', 'Interpreter', 'latex', 'fontsize', 16)


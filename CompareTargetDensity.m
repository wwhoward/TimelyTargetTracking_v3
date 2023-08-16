% Written for TTT Journal by W.W.Howard in Summer 2023
% Contact: {wwhoward}@vt.edu Wireless@VT
% For TimelyTrackingNetwork v3
% 
% Task: 
% Compare capacities for AoII metric and random

% Header
addpath(genpath(string(pwd)+'\TimelyTrackingNetwork_v3.0'))
clc;clear;close all


% Parameters
nRep = 60; 
Strats = {"distributed", "centralized_random"}; 
t_step = 0.2; 
duration = 10; % #seconds to simulate
coverage = 0.1; 
nTrackers = 20; % This is the PPP density per unit km
updateRate = 1; % rate per node 
nUpdateRates = length(updateRate); 
targetDensity = [1, 3, 10, 30, 100, 300]; 
lengthTargetDensity = length(targetDensity); 

statsAoii = {}; 
statsRandom = {}; 

parfor n = 1:nRep
    display("Rep " + string(n))

    for r = 1:lengthTargetDensity
        myTargets = targetModel(targetDensity(r)); 
        myTrackers = genTrackers(nTrackers, 'Coverage', coverage, 'updateRate', updateRate); 

        myAoiiFC = fusionCenter(myTargets, myTrackers, Strats{1}, updateRate);
        myRandomFC = fusionCenter(myTargets, myTrackers, Strats{2}, updateRate);

        for t = 0:t_step:duration
            myTargets.update(t_step); 
            for i = 1:length(myTrackers)
                myTrackers{i}.observe(myTargets, t); 
            end
            myAoiiFC.getUpdates(t); 
            myRandomFC.getUpdates(t); 
        end

        % get stats 
        
        statsAoii{r, n} = myAoiiFC.Stats; 
        statsRandom{r, n} = myRandomFC.Stats;         
    end
end

meanStatsAoii = AverageStats(statsAoii); 
meanStatsRandom = AverageStats(statsRandom); 
colors = linspecer(lengthTargetDensity); 
display_names = string(targetDensity); 

% ECDF
figure; 
for r = 1:4
    semilogx(meanStatsAoii{"SS_ECDF"}{r,2}, meanStatsAoii{"SS_ECDF"}{r,1}, 'Color', colors(r,:), 'linewidth', 2, 'DisplayName', "AoII, $\delta=$"+display_names{r})
    hold on
    semilogx(meanStatsRandom{"SS_ECDF"}{r,2}, meanStatsRandom{"SS_ECDF"}{r,1}, '--', 'Color', colors(r,:), 'linewidth', 2, 'DisplayName', "Random, $\delta=$"+display_names{r})
end
grid on
legend('Interpreter', 'latex', 'fontsize', 12, 'location', 'northwest')
xlabel('Meters', 'Interpreter', 'latex', 'fontsize', 12)
xlim([10^0, 10^3])
ylabel('Pr(Error $\leq X$)', 'Interpreter', 'latex', 'fontsize', 12)

% Targets covered
% Nodes covered
figure; hold on
for i = 1:4
    plot(meanStatsAoii{"TimeSteps"}, meanStatsAoii{"nCoveredTargets"}(i,:), '-', 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', "Covered Targets, "+display_names{i})
    plot(meanStatsAoii{"TimeSteps"}, meanStatsAoii{"nTrackedTargets"}(i,:), '--', 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', "Tracked Targets, "+display_names{i})
    plot(meanStatsRandom{"TimeSteps"}, meanStatsRandom{"nTrackedTargets"}(i,:), ':', 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', "Tracked Targets, "+display_names{i})
end
legend('Interpreter', 'latex', 'fontsize', 12, 'location', 'southeast')
xlabel('Time (s)', 'Interpreter', 'latex', 'fontsize', 12)
ylabel('Number of UAVs', 'Interpreter', 'latex', 'fontsize', 12)
title('Centralized Selection Coverage', 'Interpreter', 'latex', 'fontsize', 16)

% Nodes updated
% Ave selected nodes
figure; hold on
for i = 1:lengthTargetDensity
    tmp_mean = mean(meanStatsAoii{"nSelectedNodes"}(i,:)); 
    plot(meanStatsAoii{"TimeSteps"}, meanStatsAoii{"nSelectedNodes"}(i,:), '-', 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', display_names{i}+", mean: "+string(tmp_mean))
    plot(meanStatsRandom{"TimeSteps"}, meanStatsRandom{"nSelectedNodes"}(i,:), ':', 'Color', colors(i,:), 'linewidth', 2, 'DisplayName', display_names{i}+", mean: "+string(tmp_mean))
end
legend('interpreter', 'latex', 'fontsize', 12)
xlabel('Time (s)', 'interpreter', 'latex', 'fontsize', 12)
ylabel('Nodes', 'interpreter', 'latex', 'fontsize', 12)
title('Number of Selected Nodes', 'interpreter', 'latex', 'fontsize', 16)
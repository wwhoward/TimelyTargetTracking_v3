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
% updateRate = [0.1, 0.3, 1, 3, 10]; % rate per node 
updateRate = [0.075, 0.25, 0.75]; % rate per node, this equals 0.3, 1, 3 I think
nUpdateRates = length(updateRate); 

statsAoii = {}; 
statsRandom = {}; 

parfor n = 1:nRep
    display("Rep " + string(n))

    for r = 1:nUpdateRates
        myTargets = targetModel(30); 
        myTrackers = genTrackers(nTrackers, 'Coverage', coverage, 'updateRate', updateRate(r)); 

        myAoiiFC = fusionCenter(myTargets, myTrackers, Strats{1}, updateRate(r));
        myRandomFC = fusionCenter(myTargets, myTrackers, Strats{2}, updateRate(r));

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
colors = linspecer(4); 
display_names = string(updateRate); 
specs = {"-", "--", ":"}; 

figure; 
for r = [3,2,1]% [4, 3, 2]
    semilogx(meanStatsAoii{"SS_ECDF"}{r,2}, meanStatsAoii{"SS_ECDF"}{r,1}, specs{r}, 'Color', colors(1,:), 'linewidth', 2, 'DisplayName', "AoII, $\delta=$"+display_names{r})
    hold on
    semilogx(meanStatsRandom{"SS_ECDF"}{r,2}, meanStatsRandom{"SS_ECDF"}{r,1}, specs{r}, 'Color', colors(4,:), 'linewidth', 2, 'DisplayName', "Random, $\delta=$"+display_names{r})
end
grid on
legend('Interpreter', 'latex', 'fontsize', 12, 'location', 'northwest')
xlabel('Meters', 'Interpreter', 'latex', 'fontsize', 12)
xlim([1e-1, 5e1])
ylabel('Pr(Error $\leq X$)', 'Interpreter', 'latex', 'fontsize', 12)
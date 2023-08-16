function [trackers] = genTrackers(nTrackers, varargin)
    % GENTRACKERS Convenience function to enable instancing of targetTrackers in parfor loops
    % Written for TTT Journal by W.W.Howard in Spring 2023
    % Version information: 
    % v3.0
    % Contact: {wwhoward}@vt.edu Wireless@VT
    
    trackers = {}; 
    % Draw the actual number of nodes from a Poisson distribution
    nTrackers = poissrnd(nTrackers); 
    for n = 1:nTrackers
        trackers{n} = targetTracker(n, nTrackers, varargin{:}); 
    end
end
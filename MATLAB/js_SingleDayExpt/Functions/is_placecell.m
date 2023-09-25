function [isPC] = is_placecell(FRmeans,FRstds,FRpeaks,sparsities,coverages)
%IS_PLACECELL Determines if neuron is a place cell.
%
% From infromation on all trajectories, determines if a neuron is a place
% cell or not. Assumes 4 trajectories of the W-track maze.
%   
% Inputs:
%   FRmeans - mean occupancy normalized firing rates across each of the 4 trajectories.
%   FRstds - standard deviations of the occupancy normalized firing rates.
%   FRpeaks - peak occupancy norm FR for each trajectory.
%   sparsities - measures of spatial sparsity for 4 trajectories.
%   coverages - measures of spatial coverage for 4 trajectories.
%
% Outputs:
%   isPC - list of logical 1s and 0s if the cell is a place cell or not for
%   each trajectory.

% Assume cell is not a place cell on any trajectory.
isPC = [0, 0, 0, 0];

for tr = 1:4
    peakRate = FRpeaks(tr);
    meanRate = FRmeans(tr);
    stdRate = FRstds(tr);
    spar = sparsities(tr);
    covg = coverages(tr);

    % Peak FR must be above 2 stds from the mean.
    if peakRate > meanRate+2*stdRate
        % Mean rate must be below 10 Hz.
        if meanRate < 100
            % Peak rate must be above 3 Hz
            if peakRate > 3
                % Sparsitty must be below 0.5
                if spar < 0.5
                    % coverage must be less than 0.5
                    if covg < 0.5

                        isPC(tr) = 1;

                    end
                end
            end
        end
    end
end




end


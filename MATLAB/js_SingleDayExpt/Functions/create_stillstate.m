function [C_stillstate] = create_stillstate(C_allpos)
%CREATE_RUNSTATE    Create cell array for times of stillness.
%
% Create a cell array variable similar to that of rem01,
% sleep01, etc. from the velocity and time information in pos01.
%   Inputs:
%       C_allpos - cell array containing pos01 file data for all rats.
%       Dimensions must be 1 x (num rats) for the outer cell array, with
%       each cell holding all epochs of data for that rat as 1 x (num epochs). Under each epoch
%       must be a struct with a field 'data' which is a matrix. The first
%       column of data must be time and the fifth column velocity in cm/s.
%
%   Outputs:
%       C_stillstate - cell array containing start and end time information
%       for the stillness of each animal. Stillness is defined as having a
%       velocity less than vel_thr. Dimensions are identical
%       to C_allpos, with structs contaning the fields 'starttime',
%       'endtime', 'timerange', 'vel_thr', and 'dur_thr'.

vel_thr = 4; % Velocity threshold, above which the rat is considered running 
% and below which stationary.
% dur_thr = 0.165; % Duration of running required to be listed as an occurance
% (in seconds). The time steps in data are approximately 0.033 seconds
% apart, so requiring 5 consecutive time steps below velocity threshold
% translates to 0.165 seconds.
C_runstate = cell(1,size(C_allpos,2));

for r = 1:size(C_allpos,2)

    epochStillState = cell(1,size(C_allpos{1,r},2));

    for e = 1:size(C_allpos{1,r},2)
        
        data = C_allpos{1,r}{1,e}.data; % Matrix of data for that epoch
        times = data(:,1); % Time indicies
        vel = data(:,5); % Velocity of the rat at those times

        velSmooth = sgolayfilt(vel,2,31); % Smooths velocity data (it needs smoothing!).

        stillStruct.timerange = [times(1), times(end)]; % Create struct and assign the time range.
        stillStruct.vel_thr = vel_thr;
        % stillStruct.dur_thr = dur_thr;
        
        % Arrays to hold start and end times of still periods.
        starttimes = []; endtimes = [];
        instill = 0; % 1 if currently still, 0 if not.

        if velSmooth(1) < vel_thr % Assign first start time if above threshold.
            starttimes(1) = times(1);
            instill = 1;
        end
        
        % Demarcate start and end times.
        for t = 1:(size(times,1)-1)
            
            if ~instill && velSmooth(t+1) < vel_thr % start times
                starttimes = [starttimes; times(t+1)];
                instill = 1;
            end

            if instill && velSmooth(t+1) >= vel_thr % end times
                endtimes = [endtimes; times(t+1)];
                instill = 0;
            end
        end

        % There should be an equal number of start and end times, or at
        % most starttimes being one index longer if the epoch ended during
        % a run period. To fix this if it occurs, add an endtime as the
        % last time in epoch.
        if abs(length(starttimes) - length(endtimes)) > 1
            error("Number of starttimes and endtimes recorded are somehow unequal for" + ...
                "Rat: %d, Epoch: %d",r,e)
        elseif (length(starttimes) - length(endtimes)) == 1
            endtimes(end+1) = times(end);
        end

        % % Finally, remove run occurances that are too short. I may not need
        % % this now that the velocity is filtered, but it doesn't hurt to
        % % leave it here.
        % islong = (endtimes-starttimes) >= dur_thr;
        % starttimes = starttimes(islong);
        % endtimes = endtimes(islong);

        % Store times for this epoch
        stillStruct.starttime = starttimes;
        stillStruct.endtime = endtimes;

        epochStillState{1,e} = stillStruct;
    end

    C_stillstate{1,r} = epochStillState;

end





end
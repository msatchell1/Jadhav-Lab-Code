function [S_ripMtcs] = calc_ripplemetrics(brainAreas, C_allriptimes, C_allspikes)
%CALC_RIPPLEMETRICS Calculates the following metrics for ripple analysis:
%...
%
% Inputs:
%
% Outputs:
% S_ripMtcs - a struct containing information for each rat and ripple.
% Fields of the struct can be 1 x (num rats) cell arrays, within which are
% matrices of (num ripples) x (num epochs), or cell arrays of the same
% dimension.








S_ripMtcs.STs_inRip = cell(length(brainAreas),size(C_allriptimes,2));


% Main loop to calculate mean firing rate for all states in all brain areas.
for a = 1:length(brainAreas)
    fprintf("Working on %s... \n",brainAreas{a})

    for r = 1:size(C_allriptimes,2) 

        for e = 1:length(C_allriptimes{1,r})
            

            epochData = C_allriptimes{1,r}{1,e}; % Data on a single rat and epoch.
            nrnsAlltets = [C_allspikes{1,r}{1,e}{:}]; % Combining nrn spike data from all tets. 
            STs_inRip = cell(size(epochData.starttime,1),length(nrnsAlltets)); % Holds spike times for each neuron and each ripple.

            for o = 1:size(epochData.starttime,1) % Loop through occurances of ripples.

                for nrn = 1:length(nrnsAlltets) % Loop through all neurons in that rat/epoch

                    % isfield returns false if the struct does not exist
                    if isfield(nrnsAlltets{nrn},'data') && ~isempty(nrnsAlltets{nrn}.data)...
                            && strcmp(nrnsAlltets{nrn}.area, brainAreas{a})
        
                        STs = nrnsAlltets{nrn}.data(:,1); % Time of all spikes for that nrn.

                        % Logical 1s and 0s determining if spikes fall within the state occurance window.
                        inOcc = (epochData.starttime(o) <= STs) & (STs <= epochData.endtime(o));
                        STs_inRip{o,nrn} = STs(inOcc);
                        
                    end
                end
            end
        end
        % NOTE THE ST DATA IS NOT SPLIT UP BY EPOCH WHICH I WILL WANT TO DO
        
        % Now that all data across the epochs has been extracted...

        S_ripMtcs.STs_inRip{a,r} = STs_inRip;
        
            
        occTimes = [epochData.starttime, epochData.endtime];
        occDurs = epochData.endtime - epochData.starttime;

        
        % I can't believe matlab lets me make a cell array insert data
        % within it all in one line like this...
        S_ripMtcs.durs{1,r} = occDurs;
        S_ripMtcs.startEndTimes{1,r} = occTimes;
           


    fprintf("       Ripple metrics for rat %d completed. \n",r)
    end
    
end












end
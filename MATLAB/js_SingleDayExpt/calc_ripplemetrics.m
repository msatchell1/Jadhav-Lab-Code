function [S_ripMtcs] = calc_ripplemetrics(brainAreas, C_allriptimes, C_allspikes)
%CALC_RIPPLEMETRICS Calculates metrics for ripple analysis.
%
% Inputs:
%   brainAreas - cell array of string, brain regions of interest.
%
%   C_allriptimes - cell array containing loaded ripple time data from each rat.
%   This is used to determine the occurance times for each ripple.
%
%   C_allspikes - cell array containing loaded spike data from each rat.
%   This is used to get the compare spike times of each neuron to ripple times.
%
% Outputs:
%   S_ripMtcs - a struct containing information for each rat and ripple.
%   Fields are:
%       STs_inRip - Nested cell arrays of dimensions: 1 x (num rats) / 1 x 
%       (num epochs) / (num ripples) x (num neurons), where each cell in
%       the (num ripples) x (num neurons) array contains the spike times
%       for that neuron during that ripple.









% Main loop to calculate mean firing rate for all states in all brain areas.
for a = 1:length(brainAreas)
    fprintf("Working on %s... \n",brainAreas{a})

    for r = 1:size(C_allriptimes,2) 

        for e = 1:length(C_allriptimes{1,r})
            

            epochData = C_allriptimes{1,r}{1,e}; % Data on a single rat and epoch.
            nrnsAlltets = [C_allspikes{1,r}{1,e}{:}]; % Combining nrn spike data from all tets. 
            STs_inRip = cell(size(epochData.starttime,1),length(nrnsAlltets)); % Holds spike times for each neuron and each ripple.
            numClstNrns = 0; % Counter for the number of clustered neurons on all the tetrodes in one brain region (for that epoch).
            
            for nrn = 1:length(nrnsAlltets) % Loop through all neurons in that rat/epoch
                 % isfield returns false if the struct does not exist
                if isfield(nrnsAlltets{nrn},'data') && ~isempty(nrnsAlltets{nrn}.data)...
                            && strcmp(nrnsAlltets{nrn}.area, brainAreas{a})

                    numClstNrns = numClstNrns + 1;
                    STs = nrnsAlltets{nrn}.data(:,1); % Time of all spikes for that nrn.

                    for o = 1:size(epochData.starttime,1) % Loop through occurances of ripples.
                   
                        % Logical 1s and 0s determining if spikes fall within the state occurance window.
                        inOcc = (epochData.starttime(o) <= STs) & (STs <= epochData.endtime(o));
                        STs_inRip{o,nrn} = STs(inOcc);
                        
                    end
                end
            end

            S_ripMtcs.spikeTimes{a,r}{1,e} = STs_inRip;
    
            occTimes = [epochData.starttime, epochData.endtime];
            occDurs = epochData.endtime - epochData.starttime;
    
            % Start and end times of ripples
            S_ripMtcs.startEndTimes{1,r}{1,e} = occTimes;
            % Duration of ripples
            S_ripMtcs.durs{1,r}{1,e} = occDurs;

            % The number of spikes per ripple and neuron
            numSpikes = cellfun(@numel,STs_inRip);
            S_ripMtcs.numSpikes{a,r}{1,e} = numSpikes;

            % The number of clustered neurons
            S_ripMtcs.numClstNrns{a,r}{1,e} = numClstNrns;

            % The number of participating cells per ripple (out of all
            % those clustered for this region).
            hasSpikes = ~cellfun(@isempty,STs_inRip);
            S_ripMtcs.fracNrns{a,r}{1,e} = sum(hasSpikes,2)/numClstNrns;
            
        end
        
        % Now that all data across the epochs has been extracted...

        
        
            

           


    fprintf("       Ripple metrics for rat %d completed. \n",r)
    end
    
end












end
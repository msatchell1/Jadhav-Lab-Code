function [C_FRs] = label_interneurons(C_allspikes,brainAreas)
%FIND_INTERNEURONS Label interneurons based on mean firing rate.
%   Detailed explanation goes here


% This is how Shantanu, Frank 2016 sorted out interneurons:
% Fast-spiking (FS) putative inhibitory interneurons (n = 10) were
% detected by a consistent spike width and firing rate criterion across the two tasks (spikewidth < 0.3 ms and firing rate > 17 Hz) and excluded from further analysis. A further 31
% cells were excluded as they did not have sufficient number of spikes (< 50) for
% quantification of SWR modulation. This analysis therefore included a total of 312 PFC
% neurons. Similarly for CA1 neurons, FS neurons (n = 22, spike-width < 0.3 ms and firing
% rate > 7 Hz).


dataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt/matclust files';

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
loadRats = {'JS14'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01','pos01'
filetypes = {};




%%
for a = 1:size(brainAreas,2)

    for r = 1:size(C_allspikes,2)
        % To hold mean firing rate data across all epochs for a single rat.
        ratFRmeans = NaN(length([C_allspikes{1,r}{1,1}{:}]), length(C_allspikes{1,r}));

        for e = 1:size(C_allspikes{1,r},2)


            nrnsAlltets = [C_allspikes{1,r}{1,e}{:}]; % Combining nrn spike data from all tets.

            epochisPyr = zeros(1,size(nrnsAlltets,2));
            epochisInh = zeros(1,size(nrnsAlltets,2));
            epochFRmeans = NaN(1,size(nrnsAlltets,2));

            for nrn = 1:length(nrnsAlltets)
        
                
                if isfield(nrnsAlltets{nrn},'meanrate') && strcmp(nrnsAlltets{nrn}.area, brainAreas{a})

                    epochFRmeans(1,nrn) = nrnsAlltets{nrn}.meanrate;

                end



                    
                
            end

            ratFRmeans(:,e) = epochFRmeans; % Turns row vector into column vector. 
            % Each column is an epoch, rows are neurons.
            



        end


        % Find neurons that fire during at least one epoch
        numNaNepochs = sum(isnan(ratFRmeans),2);
        isActiveNrn = numNaNepochs < size(C_allspikes{1,r},2); % 1 for nrns that have any 
        % non-NaN epochs, and 0 otherwise.

        count = 0;
        ratBothMeans = [];
        figure;
        hold on;
        
        for nrn = 1:size(ratFRmeans,1)
            nrnFRsNaN = ratFRmeans(nrn,:);
            nrnFRs = nrnFRsNaN(~isnan(nrnFRsNaN));
            pltEpochs = find(~isnan(nrnFRsNaN));
            % if ~isempty(nrnFRs) && any(nrnFRs > 7) && any(nrnFRs <= 7) 
            if ~isempty(nrnFRs) && all(nrnFRs < 7)
            % if ~isempty(nrnFRs)

                count = count + 1;
                ratBothMeans(nrn,:) = nrnFRsNaN;
                plot(pltEpochs,nrnFRs,'*-')

            end
        end
        xlabel("Epochs")
        ylabel("Mean FR")
        title(sprintf("Region %s Rat %d | count = %d/%d \n",brainAreas{a},r,count,sum(isActiveNrn,1)))
        % fprintf("%s Rat %d | above&below7Hz = %d/%d \n",brainAreas{a},r,count,sum(isActiveNrn,1))
        
        pause
        close all

        C_FRs{1,a}{1,r} = ratFRmeans;

    end

    

end














end


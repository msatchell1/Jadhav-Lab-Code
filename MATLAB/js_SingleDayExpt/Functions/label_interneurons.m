function [C_FRs] = label_interneurons(C_allspikes,brainAreas,loadRats)
%FIND_INTERNEURONS Label interneurons based on mean firing rate.
%   Detailed explanation goes here


% This is how Shantanu, Frank 2016 sorted out interneurons:
% For PFC: Fast-spiking (FS) putative inhibitory interneurons (n = 10) were
% detected by a consistent spike width and firing rate criterion across the 
% two tasks (spikewidth < 0.3 ms and firing rate > 17 Hz) and excluded from further analysis. A further 31
% cells were excluded as they did not have sufficient number of spikes (< 50) for
% quantification of SWR modulation. This analysis therefore included a total of 312 PFC
% neurons. Similarly for CA1 neurons, FS neurons (n = 22, spike-width < 0.3 ms and firing
% rate > 7 Hz).


dataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt/matclust files';

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01','pos01'
filetypes = {};

C_matclust = load_matclust_data(dataDir,loadRats);


%%
% Overall plan: 
% I have the data to label the cells now in the form of the matclust files,
% which are essentially the raw output of the clustering program Justin used.
% The data is organized into clusters, each representing a cell, and each 
% with indices corresponding to the spikes of that cell. The actual spike 
% waveforms are stored in “waves” files, and I use the indices of the clusters
% to access the correct waveforms for that cell. In this way I can take an 
% average of all the waveforms for each cell to get an idea of what the average
% spike looks like, but I need to make sure I know which tetrode to get the 
% waveform off of, because the waves file has the waveform on all 4 tetrodes.
% Once this is done I can calculate the average spike width and determine what
% are interneurons like Shantanu, Frank 2016: in PFC having a spikewidth < 0.3 ms
% and FR > 17 Hz. In CA1 having spike width < 0.3 ms and firing rate > 7 Hz. 

allSpikewidths = [];
for r = 1:size(C_matclust,2)
    
    for mt = 1:size(C_matclust,1) % NOTE this is not the tetrode number because these
        % are out of order.

        clustattrib = C_matclust{mt,r}.clustattrib;
        allWaveforms = C_matclust{mt,r}.waves;
        clusters = clustattrib.clusters; % Each cluster is a cell
        for nrn = 1:size(clusters,2)
            % Long column vector of spike indices for the waves file.
            spikeIndxs = clusters{1,nrn}.index;
            % The waveforms for all spikes of this cell as measured on all
            % 4 electrodes of the tetrode.
            waveforms4 = allWaveforms(:,:,spikeIndxs);
            meanWF4 = mean(waveforms4,3,'omitnan'); % Takes mean shape across all spikes
            % figure
            % plot(meanWF4) % Empty plots are clusters with no spikes.
            % pause
            % close all

            % Return peak heights for all electrode channels
            pk4 = max(meanWF4,[],1);
            [maxPk4, maxIdx4] = max(pk4); % largest peak among electrodes
            meanWF = meanWF4(:,maxIdx4);
            [WFmax,maxIdx] = max(meanWF);
            [WFmin, spikewidth] = min(meanWF(maxIdx:end)); % Note the index here
            % is spikewidth because it measures how many indices after
            % maxIdx that the trough occurs.
            spikewidth_s = spikewidth/30000; % Assumes 30,000 Hz sampling rate.
            
            allSpikewidths(end+1) = spikewidth_s;


            % NEXT to do:
            % 1) figure out how to link matclust tetrode order (mt) to the
            % actual order in C_allspikes.
            % 2) Make histograms of spikewidths in CA1 and PFC, check that
            % I agree with the decision line on what to call interneurons.
            % 3) Assign cells their spikewidth and classification as pyramidal or
            % interneuron.
        end
    end
end

figure;
histogram(allSpikewidths,25)


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

function [C_FRs] = label_interneurons(brainAreas,loadRats)
%FIND_INTERNEURONS Label interneurons based on mean firing rate.
%   Detailed explanation goes here


% This is how Shantanu, Frank 2016 sorted out interneurons:
% For PFC: "Fast-spiking (FS) putative inhibitory interneurons (n = 10) were
% detected by a consistent spike width and firing rate criterion across the 
% two tasks (spikewidth < 0.3 ms and firing rate > 17 Hz) and excluded from further analysis. A further 31
% cells were excluded as they did not have sufficient number of spikes (< 50) for
% quantification of SWR modulation. This analysis therefore included a total of 312 PFC
% neurons. Similarly for CA1 neurons, FS neurons (n = 22, spike-width < 0.3 ms and firing
% rate > 7 Hz)."




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


% Load spiking data so cell type and spikewidth attributes can be added to
% the neuron structs.
dataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
loadRats = 
% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01','pos01'
filetypes = {'spikes01'};

C_alldata = load_data(dataDir,loadRats,filetypes);

spikes_idx = find(contains(filetypes,'spikes01')); 
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end

C_allspikes = C_alldata(spikes_idx,:);


% Edit to folder of matclust files for each rat. Name of each subfolder
% should be '*Rat name*.matclust'
matclustDataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt/matclust files';

allNrnStructs = {};
for r = 1:size(loadRats,2)
    
    for mt = 1:size(C_allspikes{1,r}{1,1},2) % NOTE this is not the tetrode number because these
        % are out of order.

        % Loads matclust data from one rat and one tetrode.
        matclustTet = load_matclust_data(matclustDataDir,loadRats{1,r},mt);

        if isempty(matclustTet)

            fprintf("No clustered cells on tet: %d",mt)

        clustattrib = C_matclust{mt,r}.clustattrib;
        allWaveforms = C_matclust{mt,r}.waves;
        clusters = clustattrib.clusters; % Each cluster is a cell
        
        % Now find extract the tetrode number from the filenames in
        % clustattrib.
        startStr = "_nt"; % target string before tetrode number
        endStr = ".mat"; % string after tet num
        startIdx = strfind(clustattrib.currentfilename,startStr);
        endIdx = strfind(clustattrib.currentfilename,endStr);
        tetnum = str2double(clustattrib.currentfilename(tgtIdx+strlength(tgtStr):endIdx-1));

        for nrn = 1:size(clusters,2)
            nrnStruct = struct([]);
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
            
            nrnStruct.spikewidth = spikewidth_s;

            % Get the neuron's mean rate averaged across all epochs
            nrnFRepochs = nan(17);
            for e = 1:17
                if ~isempty(C_allspikes{1,r}{1,e}(1,tetnum))
                    if ~isempty(C_allspikes{1,r}{1,e}{1,tetnum}(1,nrn))
                        nrnEpoch= C_allspikes{1,r}{1,e}{1,tetnum}{1,nrn};
                        nrnFRepochs(e) = nrnEpoch.meanrate; % Mean rate of each epoch
                        nrnStruct.area = nrnEpoch.area;

                    end
                end
            end
            meanrate = mean(nrnFRepochs,'omitnan');
            nrnStruct.meanrate = meanrate;
            
            allNrnStructs{mt,r}{nrn} = nrnStruct;


            % NEXT to do:
            % 2) Make histograms of spikewidths in CA1 and PFC, check that
            % I agree with the decision line on what to call interneurons.
            % 3) Assign cells their spikewidth and classification as pyramidal or
            % interneuron.
        end
    end
end



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


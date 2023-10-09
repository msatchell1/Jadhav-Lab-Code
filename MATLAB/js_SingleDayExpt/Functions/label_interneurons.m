function [] = label_interneurons(dataDir,loadRats)
%FIND_INTERNEURONS Label interneurons based on mean firing rate.
%   Loads the spikes01 file for rats in loadRats as well as the spike
%   waveforms for all neurons of that rat and labels ineterneurons based on
%   spike width and mean firing rate across all epochs. Saves this
%   information to the spikes01 file.


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



filetypes = {'spikes01'};

C_alldata = load_data(dataDir,loadRats,filetypes);

spikes_idx = find(contains(filetypes,'spikes01')); 
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end

C_allspikes = C_alldata(spikes_idx,:);


% Edit to folder of matclust files for each rat. Name of each subfolder
% should be '*Rat name*.matclust'
matclustDataDir = dataDir + '/matclust files';


for r = 1:size(loadRats,2)
    fprintf("Processing matclust data for %s \n",loadRats{1,r})
    
    for tetnum = 1:size(C_allspikes{1,r}{1,1},2) % Loop through all tetrodes for rat.

        % Loads matclust data from one rat and one tetrode.
        matclustTet = load_matclust_data(matclustDataDir,loadRats{1,r},tetnum);

        % % Now find extract the tetrode number from the filenames in
        % % clustattrib.
        % startStr = "_nt"; % target string before tetrode number
        % endStr = ".mat"; % string after tet num
        % startIdx = strfind(clustattrib.currentfilename,startStr);
        % endIdx = strfind(clustattrib.currentfilename,endStr);
        % tetnum = str2double(clustattrib.currentfilename(tgtIdx+strlength(tgtStr):endIdx-1));

        if any(~isfield(matclustTet,{'clustattrib','waves'})) % Note the strings MUST be character vectors.
            % fprintf("Missing or incomplete data on tet: %d \n",tetnum)
            
        else
            clustattrib = matclustTet.clustattrib;
            allWaveforms = matclustTet.waves;
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
                
    
                % Get the neuron's mean rate averaged across all epochs
                nrnFRepochs = nan(1,17);
                for e = 1:17
                    if ~isempty(C_allspikes{1,r}{1,e}(1,tetnum))
                        if isfield(C_allspikes{1,r}{1,e}{1,tetnum}{1,nrn},'meanrate') &&...
                                isfield(C_allspikes{1,r}{1,e}{1,tetnum}{1,nrn},'area')

                            nrnEpoch= C_allspikes{1,r}{1,e}{1,tetnum}{1,nrn};
                            nrnFRepochs(1,e) = nrnEpoch.meanrate; % Mean rate of each epoch
                            nrnStruct.area = nrnEpoch.area;



                        end
                    end
                end
                
                acrossEpochFR = mean(nrnFRepochs,2,'omitnan');


                for e = 1:17
                    if ~isempty(C_allspikes{1,r}{1,e}(1,tetnum))
                        if isfield(C_allspikes{1,r}{1,e}{1,tetnum}{1,nrn},'meanrate') &&...
                                isfield(C_allspikes{1,r}{1,e}{1,tetnum}{1,nrn},'area')
                            

                            % Plotting spikewidth vs mean rate averaged across all epochs shows that
                            % the data clusters reasonably well, indicating that I can probably trust
                            % averaging FRs across all epochs to identify interneurons.
                            if strcmp(nrnEpoch.area, 'CA1')
                                
                                % Finds interneurons. The minimum necessary spike width increases as
                                % the mean FR increases with a slope of 0.01 ms/Hz, such that at 27 Hz
                                % the threshold is a 0.5 ms spikewidth.
                                if acrossEpochFR >= 7 && 1000*spikewidth_s <= (0.3+.01*(acrossEpochFR-7))
                                    nrnType = 'Int';
                                else
                                    nrnType = 'Pyr';
                                end
                            
                            elseif strcmp(nrnEpoch.area,'PFC')
                            
                                % Does the same for PFC neurons, which are less clearly defined. 
                                if acrossEpochFR >= 11 && 1000*spikewidth_s <= (0.3+.01*(acrossEpochFR-11))
                                    nrnType = 'Int';
                                else
                                    nrnType = 'Pyr';
                                end
                            end
            
                            % Save changes back to neuron struct
                            C_allspikes{1,r}{1,e}{1,tetnum}{1,nrn}.epochFRs = nrnFRepochs;
                            C_allspikes{1,r}{1,e}{1,tetnum}{1,nrn}.acrossEpochFR = acrossEpochFR;
                            C_allspikes{1,r}{1,e}{1,tetnum}{1,nrn}.spikeWidth = spikewidth_s;
                            C_allspikes{1,r}{1,e}{1,tetnum}{1,nrn}.nrnType = nrnType;

                        end
                    end
                end



            end
        end
    end

    % Save data from each rat to overwrite original file.
    spikes = C_allspikes(1,r);

    short_name = loadRats{r};
    chop_idx = strfind(loadRats{r},'_') - 1;
    if ~isempty(chop_idx)
        short_name = loadRats{r}(1:chop_idx); % Gets the first characters of the rat's name before an '_'.
        % So far this is only needed for ER1_NEW to remove the '_NEW'.
    end

    save(string(fullfile(dataDir,sprintf("%s_direct/%sspikes01.mat",loadRats{r},short_name))),"spikes")
    fprintf("Saved file: %sspikes01 \n",short_name)
end






end


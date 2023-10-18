function [] = add_structfields()
%ADD_STRUCTFIELDS   add fields to structs within a given data file across
%all rats and epochs. Overwrites the existing saved animal file and
%replaces it with the new version.
% Michael Satchell 08/30/2023

% Load data
dataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
loadRats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01','pos01'
filetypes = {'spikes01','linfields','cellinfo'};

C_alldata = load_data(dataDir,loadRats,filetypes); % Load data to struct




spikes_idx = find(contains(filetypes,'spikes01')); % Gets row index of spike data in C_alldata.
% Note that indexing like this is only supposed to give one index as a
% result.
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end

cellinfo_idx = find(contains(filetypes,'cellinfo'));
if isempty(cellinfo_idx)
    error("cellinfo data must be loaded to run this analysis.")
end

linf_idx = find(contains(filetypes,'linfields'));
if isempty(linf_idx)
    error("linfields data must be loaded to run this analysis.")
end


% % Extract spike data
% C_allspikes = C_alldata(spikes_idx,:);
% 
% % Cell info data so I can determine brain region of each nrn. I could do a
% % separate loop here and add an 'area' field to the neurons in C_allspikes.
% % I could also just treat C_allinfo like C_allspikes and flatten the cells
% % across tetrodes into one array inside the main for-loop, then grab the
% % area from C_allinfo. For now I am going to choose the first option.
% C_allinfo = C_alldata(cellinfo_idx,:);
% 
% % Linfield data
% C_alllinf = C_alldata(linf_idx,:);


% behper_idx = find(contains(filetypes,'behavperform'));
% if isempty(behper_idx)
%     error("behavperform data must be loaded to run this analysis.")
% end
% 
% 
% % Extract performance data
% C_behper = C_alldata(behper_idx,:);


% -------------------------------------------------------------------


% %% Add 'area', 'traj_x', and 'isPC' fields to spikes file.
% % Add 'area' field to the spike data for later sorting. Also add place
% % field data such as mean and std FR, sparsity, and spatial coverage for each trajectory.
% % linfields must be loaded to add the placefield data.
% for r = 1:size(C_allspikes,2)
%     for e = 1:size(C_allspikes{1,r},2)
%         for tet = 1:size(C_allspikes{1,r}{1,e},2)
%             if ~isempty(C_allspikes{1,r}{1,e}{1,tet})
%                 for nrn = 1:size(C_allinfo{1,r}{1,e}{1,tet},2)
%                     if ~isempty(C_allspikes{1,r}{1,e}{1,tet}{1,nrn}) 
% 
%                         % If the neuron exists in the spike file it should
%                         % exist in the cellinfo file.
%                         nrnCellinfo = C_allinfo{1,r}{1,e}{1,tet}{1,nrn};
% 
% 
%                         if isfield(nrnCellinfo,'area')
%                             C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.area = nrnCellinfo.area; % Create field in spike struct.       
%                         end
% 
%                         % Linfields data is only available for behavioral
%                         % epochs (even). Existence of neurons does not
%                         % quite match that of C_allspikes.
%                         if size(C_alllinf{1,r},2) >= e && ~isempty(C_alllinf{1,r}{1,e})
%                             if ~isempty(C_alllinf{1,r}{1,e}{1,tet}) && size(C_alllinf{1,r}{1,e}{1,tet},2) >= nrn
%                                 if ~isempty(C_alllinf{1,r}{1,e}{1,tet}{1,nrn})
%                                     trajsLinf = C_alllinf{1,r}{1,e}{1,tet}{1,nrn}; % Trajectories for a single neuron.
%                                     for tr = 1:size(trajsLinf,2)
% 
%                                         trData = trajsLinf{tr};
%                                         binRate = trData(:,7)./trData(:,6); % Firing rate normalized by occupancy
%                                         FRmean = mean(binRate,1,'omitnan');
%                                         FRstd = std(binRate,1,'omitnan');
% 
%                                         % calc sparsity
%                                         totOcc = sum(trData(:,6),1); % Total time spent on track
%                                         probOcc = trData(:,6)./totOcc; % Probability of occupying each spot on the trajectory
%                                         sum1 = sum(probOcc.*binRate);
%                                         sum2 = sum(probOcc.*(binRate.^2));
%                                         sparsity = (sum1^2)/sum2;
% 
%                                         % calc spatial coverage
%                                         numBins = size(trData,1); % number of spatial bins this trajectory is split into.
%                                         peakRate = max(binRate);
%                                         isAboveThr = binRate > 0.25*peakRate;
%                                         coverage = sum(isAboveThr)/numBins; % Fraction of spatial bins with at least 25%
%                                         % of the peak firing rate on that trajectory.
% 
% 
% 
%                                         % Assign to fields in struct
%                                         C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_FRmeans(tr) = FRmean;
%                                         C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_FRstds(tr) = FRstd;
%                                         C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_sparsity(tr) = sparsity;
%                                         C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_coverage(tr) = coverage;
%                                         C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.traj_FRpeaks(tr) = peakRate;
% 
%                                     end
%                                     spStruct = C_allspikes{1,r}{1,e}{1,tet}{1,nrn};
%                                     isPC = is_placecell(spStruct.traj_FRmeans, spStruct.traj_FRstds, spStruct.traj_FRpeaks,...
%                                         spStruct.traj_sparsity, spStruct.traj_coverage);
%                                     C_allspikes{1,r}{1,e}{1,tet}{1,nrn}.isPlaceCell = isPC;
% 
%                                 end
%                             end
%                         end
% 
%                     end
%                 end
%             end
%         end
%     end
% 
%     % Save data from each rat to overwrite original file.
%     spikes = C_allspikes(1,r);
% 
%     short_name = loadRats{r};
%     chop_idx = strfind(loadRats{r},'_') - 1;
%     if ~isempty(chop_idx)
%         short_name = loadRats{r}(1:chop_idx); % Gets the first characters of the rat's name before an '_'.
%         % So far this is only needed for ER1_NEW to remove the '_NEW'.
%     end
% 
%     fileName = "spikes01.mat";
%     save(sprintf(dataDir+"/%s_direct/%s"+fileName,loadRats{r},short_name),"spikes")
%     fprintf("Saved file: %s"+fileName+"\n",short_name)



% behEpochs = 2:2:17;
% restEpochs = 1:2:17;
% 
% %% Add additional performance calculations to behavperform file
% for r = 1:size(C_behper,2)
%     % Calculate performance by epoch
%     fracCorrIn = zeros(1,size(behEpochs,2));
%     fracCorrOut = zeros(1,size(behEpochs,2));
%     fracCorr = zeros(1,size(behEpochs,2));
%     numInTrials = zeros(1,size(behEpochs,2));
%     numOutTrials = zeros(1,size(behEpochs,2));
%     numTrials = zeros(1,size(behEpochs,2));
% 
%     for be = 1:size(behEpochs,2)
%         % 1s and 0s of trials for one epoch
%         isCorrOut = C_behper{1,r}.outreward(...
%             C_behper{1,r}.dayouttrials(be,1):C_behper{1,r}.dayouttrials(be,2));
%         isCorrIn = C_behper{1,r}.inreward(...
%             C_behper{1,r}.dayintrials(be,1):C_behper{1,r}.dayintrials(be,2));
% 
%         fracCorrIn(1,be) = sum(isCorrIn)/size(isCorrIn,1);
%         fracCorrOut(1,be) = sum(isCorrOut)/size(isCorrOut,1);
% 
%         isCorr = [isCorrOut; isCorrIn]; % Combine out and in trials for that epoch
%         fracCorr(1,be) = sum(isCorr)/size(isCorr,1); % Calculate fraction of correct trials
% 
%         numInTrials(1,be) = size(isCorrIn,1);
%         numOutTrials(1,be) = size(isCorrOut,1);
%         numTrials(1,be) = size(isCorr,1);
%     end
% 
%     C_behper{1,r}.eFracCorrIn = fracCorrIn;
%     C_behper{1,r}.eFracCorrOut = fracCorrOut;
%     C_behper{1,r}.eFracCorr = fracCorr;
%     C_behper{1,r}.eNumInTrials = numInTrials;
%     C_behper{1,r}.eNumOutTrials = numOutTrials;
%     C_behper{1,r}.eNumTrials = numTrials;
% 
%     % Save data from each rat to overwrite original file.
%     behper = C_behper(1,r);
% 
%     short_name = loadRats{r};
%     chop_idx = strfind(loadRats{r},'_') - 1;
%     if ~isempty(chop_idx)
%         short_name = loadRats{r}(1:chop_idx); % Gets the first characters of the rat's name before an '_'.
%         % So far this is only needed for ER1_NEW to remove the '_NEW'.
%     end
%     fileName = "behavperform.mat";
%     save(sprintf(dataDir+"/%s_direct/%s"+fileName,loadRats{r},short_name),"behper")
%     fprintf("Saved file: %s"+fileName+"\n",short_name)
% 
% end


end
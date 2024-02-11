function [] = create_nrninfo_MES(dataDir,loadRats)
% Creates cellinfo_MES.mat structs that contain information pertaining to
% each neuron across all epochs.
%
% This function saves a file to each rat's folder called nrninfo_MES.mat.
% It is a cell array with dimensions 1 x (num rats), within which is a cell
% array (tot num neurons) x 1. Each cell contains a struct with that neuron's
% information.


%% Assign cell type in the spikes files for each rat, because this is used later.
% % Only run this if the structs in spikes01.mat files don't seem to have the
% % nrnType field. Cells will have an empty "type" field in nrninfo_MES if
% % the nrnType field is missing from spikes01.mat.
% label_interneurons(dataDir,loadRats);

%% Load data
filetypes = {'spikes01','behavperform','linfields','linpos01','combstates_MES'};

C_alldata = load_data(dataDir,loadRats,filetypes);

spikes_idx = find(contains(filetypes,'spikes01')); 
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end
behper_idx = find(contains(filetypes,'behavperform'));
if isempty(behper_idx)
    error("behavperform data must be loaded to run this analysis.")
end
if ~any(contains(filetypes,'linfields'))
    error("linfields data must be loaded to run this analysis.")
end
if ~any(contains(filetypes,'linpos01'))
    error("linpos01 data must be loaded to run this analysis.")
end
if ~any(contains(filetypes,'combstates_MES'))
    error("combstates_MES data must be loaded to run this analysis.")
end

C_allspikes = C_alldata(spikes_idx,:);
C_behper = C_alldata(behper_idx,:);
C_linf = C_alldata(find(contains(filetypes,'linfields')),:);
C_linpos = C_alldata(find(contains(filetypes,'linpos01')),:);
C_combstates = C_alldata(find(contains(filetypes,'combstates_MES')),:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;



%% Meat and potatoes


fprintf("Creating nrninfo_MES files... \n \n")
for r = 1:size(loadRats,2)

    nrninfo_MES = cell(1,1);

    % All these arrays are of dims (num nrns) x (num epochs). % The first 
    % epoch can be used as a template because the number of neurons is constant
    % across epochs in spikes01. 
    epochSpikeData = cell(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochHasSpikeData = false(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2)); % empty logical array.
    tetNums = zeros(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochFR = NaN(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochSW = NaN(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochTimeRange = cell(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochArea = strings(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    meanSW = NaN(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    nrnType = strings(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochTrajFR = cell(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochTrajSTD = cell(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochTrajSparsity = cell(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochTrajCoverage = cell(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochTrajPeaks = cell(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochTrajIsPC = cell(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));
    epochLinField = cell(size([C_allspikes{1,r}{1,2}{:}],2),size(C_allspikes{1,r},2));
    epochStateFR = cell(size([C_allspikes{1,r}{1,1}{:}],2),size(C_allspikes{1,r},2));

    areaExist = ones(size([C_allspikes{1,r}{1,1}{:}],2),1); 
    typeExist = ones(size([C_allspikes{1,r}{1,1}{:}],2),1); % Variables to catch any animals where at least
    % one cell doesn't have a given struct field.

    for e = 1:size(C_allspikes{1,r},2) % I want to loop through every epoch to 
        % get the data that I need (like which epochs a neuron is clustered
        % for, etc.), but I don't want to assign values to the new
        % nrninfo_MES cell array here because that does not have a
        % dimension for epochs.
        
        
        nrnCount = 1; % To count neuron number across tetrodes
        for tet = 1:size(C_allspikes{1,r}{1,e},2)

            for n = 1:size(C_allspikes{1,r}{1,e}{1,tet},2)
                tetNums(nrnCount,e) = tet;

                if isfield(C_allspikes{1,r}{1,e}{1,tet}{1,n}, 'data') && ~isempty(C_allspikes{1,r}{1,e}{1,tet}{1,n}.data)
                    epochSpikeData{nrnCount,e} = C_allspikes{1,r}{1,e}{1,tet}{1,n}.data;
                    epochHasSpikeData(nrnCount,e) = true;
                    epochFR(nrnCount,e) = C_allspikes{1,r}{1,e}{1,tet}{1,n}.meanrate;
                    epochTimeRange{nrnCount,e} = C_allspikes{1,r}{1,e}{1,tet}{1,n}.timerange;
                    
                end

                if isfield(C_allspikes{1,r}{1,e}{1,tet}{1,n}, 'area')
                    epochArea(nrnCount,e) = C_allspikes{1,r}{1,e}{1,tet}{1,n}.area;
                end
                
                if isfield(C_allspikes{1,r}{1,e}{1,tet}{1,n}, 'nrnType') % For those cells that were
                    % processed by the label_interneurons.m function.
                    meanSW(nrnCount,e) = C_allspikes{1,r}{1,e}{1,tet}{1,n}.spikeWidth;
                    nrnType(nrnCount,e) = C_allspikes{1,r}{1,e}{1,tet}{1,n}.nrnType;
                end

                if isfield(C_allspikes{1,r}{1,e}{1,tet}{1,n}, 'spikewidth')
                    epochSW(nrnCount,e) = C_allspikes{1,r}{1,e}{1,tet}{1,n}.spikewidth;
                end

                if isfield(C_allspikes{1,r}{1,e}{1,tet}{1,n}, 'isPlaceCell')
                    epochTrajFR{nrnCount,e} = C_allspikes{1,r}{1,e}{1,tet}{1,n}.traj_FRmeans;
                    epochTrajSTD{nrnCount,e} = C_allspikes{1,r}{1,e}{1,tet}{1,n}.traj_FRstds;
                    epochTrajSparsity{nrnCount,e} = C_allspikes{1,r}{1,e}{1,tet}{1,n}.traj_sparsity;
                    epochTrajCoverage{nrnCount,e} = C_allspikes{1,r}{1,e}{1,tet}{1,n}.traj_coverage;
                    epochTrajPeaks{nrnCount,e} = C_allspikes{1,r}{1,e}{1,tet}{1,n}.traj_FRpeaks;
                    epochTrajIsPC{nrnCount,e} = C_allspikes{1,r}{1,e}{1,tet}{1,n}.isPlaceCell;
                end
                

                % Collect linfield FR info. 
                % Reminder for the LinField data: 
                % There are four columns here representing the four W-track
                % trajectories:
                % Col 1: out right, col 2: in right, col 3: out left, col 4: in left.
                % Inside each trajectory lies 7 columns of data:
                % Col 1 is the distance in cm along that linearized
                % trajectory, col 2 occupancy, col 3 spike count, col 4 occupancy normalized
                % firing rate, col 5 smoothed occupancy normalized firing rate, col 6
                % smoothed occupancy, col 7 smoothed spike count.
                if e <= size(C_linf{1,r},2) % C_linf only has 16 epochs
                    if ~isempty(C_linf{1,r}{1,e}) % Exclude sleep epochs which don't have data
                        if ~isempty(C_linf{1,r}{1,e}{1,tet}) % Not all tetrodes have data.
                            % It seems to be that inside each tet, the
                            % column of each cell's data matches the column
                            % that same cell would have in C_allspikes.
                            % However, the number of cells changes with
                            % each epoch (cells get added on), so I need to
                            % ensure I don't index out of range in the
                            % early epochs.
                            if n <= size(C_linf{1,r}{1,e}{1,tet},2)
                                if ~isempty(C_linf{1,r}{1,e}{1,tet}{1,n})
                                    epochLinField{nrnCount,e} = C_linf{1,r}{1,e}{1,tet}{1,n};
                                    
                                end
                            end
                        end
                    end
                end
         
                % Collect state FR info
                allStateOccs = C_combstates{1,r}{1,e}.sepData;
                stateFRs = NaN(size(allStateOccs));
                for s = 1:size(allStateOccs,1)    
                    stateOccs = allStateOccs{s,1};
                    % A state must have occurances and neuron must have
                    % spikes to calculate state FR
                    if ~isempty(stateOccs) && isfield(C_allspikes{1,r}{1,e}{1,tet}{1,n}, 'data') && ...
                        ~isempty(C_allspikes{1,r}{1,e}{1,tet}{1,n}.data)
                        
                        spikes = C_allspikes{1,r}{1,e}{1,tet}{1,n}.data(:,1);
                        % Holds the number of spikes for a given neuron
                        % that occur during each state occurrance.
                        numSpikesInOcc = zeros(size(stateOccs,1),1);
                        % Loop through occurances of the state,
                        % counting the number of spikes that fall into
                        % each occurrance.
                        for o = 1:size(stateOccs,1)
                            % spike times greater than state start time and less
                            % than state end time.
                            isInOcc = (stateOccs(o,1) <= spikes) & (stateOccs(o,2) > spikes);
                            numSpikesInOcc(o,1) = sum(isInOcc);
                        end
                        spikeSum = sum(numSpikesInOcc);
                        stateDur = sum(stateOccs(:,2) - stateOccs(:,1)); % total duration of state
                        stateFRs(s,1) = spikeSum/stateDur;
                    end
                end
                epochStateFR{nrnCount,e} = stateFRs;


                nrnCount = nrnCount + 1;
            end
        end

    end


    % Now loop through all the neurons of that rat, assigning struct values
    for nrn = 1:size([C_allspikes{1,r}{1,1}{:}],2)
        S_nrn = struct();
        S_nrn.ID = nrn;

        % Remove empty entries then check that all areas are the same for
        % brain area.
        emptyIdx = epochArea(nrn,:) == "";
        nonemptyAreas = epochArea(nrn, ~emptyIdx);
        if ~isempty(nonemptyAreas) && all(nonemptyAreas == nonemptyAreas(1))
            S_nrn.area = nonemptyAreas(1);
        else
            S_nrn.area = "";
            areaExist(nrn) = 0;
        end
        
        % Same thing for neuron type
        emptyIdxT = nrnType(nrn,:) == "";
        nonemptyTypes = nrnType(nrn, ~emptyIdxT);
        if ~isempty(nonemptyTypes) && all(nonemptyTypes == nonemptyTypes(1))
            S_nrn.type = nonemptyTypes(1);
        else
            S_nrn.type = "";
            typeExist(nrn) = 0;
        end

        % By taking the mean I will know when there is a problem if a
        % tetrode number is not an integer.
        S_nrn.tet = mean(tetNums(nrn,:),2);
        S_nrn.eSpikeData = epochSpikeData(nrn,:);
        S_nrn.spikeDataCols = 'time, x, y, dir, not_used, amplitude(highest variance channel), posindex';
        S_nrn.eHasSpikeData = epochHasSpikeData(nrn,:);
        S_nrn.eTimeRange = epochTimeRange(nrn,:);
        S_nrn.eFR = epochFR(nrn,:);
        S_nrn.meanFR = mean(epochFR(nrn,:),2,'omitnan');
        S_nrn.eSpikeWidth = epochSW(nrn,:);

        % Same thing for neuron type. Note that this value is not the mean
        % of eSpikeWidth because I don't know how those values were
        % calculated by Justin, and some seem  to be way off. The value of
        % meanSW here is the one derived by label_interneurons.m.
        emptyIdxSW = isnan(meanSW(nrn,:));
        nonemptySW = meanSW(nrn, ~emptyIdxSW);
        if ~isempty(nonemptySW) && all(nonemptySW == nonemptySW(1))
            S_nrn.meanSW = nonemptySW(1);
        else
            S_nrn.meanSW = NaN;
        end

        S_nrn.eTrajisPC = epochTrajIsPC(nrn,:);
        S_nrn.eTrajFR = epochTrajFR(nrn,:);
        S_nrn.eTrajSTD = epochTrajSTD(nrn,:);
        S_nrn.eTrajSparsity = epochTrajSparsity(nrn,:);
        S_nrn.eTrajCoverage = epochTrajCoverage(nrn,:);
        S_nrn.eTrajFRPeak = epochTrajPeaks(nrn,:);
        
        % Linear Field calculations
        S_nrn.eTrajLinField = epochLinField(nrn,:);
        ctrStemLen = 80; % length of center stem of W-track in cm.
        minNumSpikes = 150; % Minimum number of spikes needed to run the selectivity
        % index on an epoch
        eSelIdx2 = NaN(1,size(C_allspikes{1,r},2));
        % figure;
        % tl = tiledlayout(2,4);
        % lgdExist = 0;
        for e = find(~cellfun(@isempty,(epochLinField(nrn,:))))
            % Some cells don't have data during beh epochs
            if ~any([isempty(epochLinField{nrn,e}), isempty(S_nrn.eSpikeData{1,e})])

                % Actually finding the spikes that occur on the center
                % stem.
                % The linear distance of the animal matched with indices
                linDist = C_linpos{1,r}{1,e}.statematrix.lindist;
                linIsCtr = linDist <= ctrStemLen;
                ctrIdxs = find(linIsCtr);
                % get all spikes
                spikeData = S_nrn.eSpikeData{1,e};
                isCtrSpike = ismember(spikeData(:,7),ctrIdxs);
                ctrSpikeData = spikeData(isCtrSpike,:); % Spikes that occur on center stem
                % Identify indices of outbound trajs
                linTraj = C_linpos{1,r}{1,e}.statematrix.traj;
                isOutboundR = ismember(linTraj,[1]);
                isOutboundL = ismember(linTraj,[3]);
                outRIdxs = find(isOutboundR);
                outLIdxs = find(isOutboundL);
                % Sort spikes even further
                ctrOutRSpkData = ctrSpikeData(ismember(ctrSpikeData(:,7),outRIdxs),:);
                ctrOutLSpkData = ctrSpikeData(ismember(ctrSpikeData(:,7),outLIdxs),:);


                % Get the occupancy-normalized spike count for right and left
                % outbound trajectories on the center stem
                isCtr_right = epochLinField{nrn,e}{1,1}(:,1) <= ctrStemLen;
                isCtr_left = epochLinField{nrn,e}{1,3}(:,1) <= ctrStemLen;

                % Smoothed occupancy normalized spike count (I'm calling it FR,
                % but in reality the rate is calculated when normalizing by
                % occupancy).
                ctrRData = epochLinField{nrn,e}{1,1}(isCtr_right,:);
                ctrLData = epochLinField{nrn,e}{1,3}(isCtr_left,:);
                
                % Better selectivity index, confirmed by inspecting
                % examples.
                diffNaN = isnan(ctrRData(:,5) - ctrLData(:,5));
                
                % Only give selectivity index when there are enough spikes
                if numel(ctrOutLSpkData)+numel(ctrOutRSpkData) >= minNumSpikes
                    selIdx2 = sum(ctrRData(:,5) - ctrLData(:,5),'omitnan')/...
                        mean([sum(ctrRData(~diffNaN,5)), sum(ctrLData(~diffNaN,5))]);
                else
                    selIdx2 = NaN;
                end
                eSelIdx2(1,e) = selIdx2;

                % nexttile
                % hold on;
                % plot(ctrRData(:,1),ctrRData(:,5),'r')
                % plot(ctrLData(:,1),ctrLData(:,5),'b')
                % plot(ctrRData(:,1),ctrRData(:,5)-ctrLData(:,5),Color=[0.7,0.7,0.7])
                % yline(0,'--k')
                % plot(linDist(ctrOutRSpkData(:,7)),rand(size(ctrOutRSpkData))*3-5,'.r')
                % plot(linDist(ctrOutLSpkData(:,7)),rand(size(ctrOutLSpkData))*3-8,'.b')
                % xlabel("Linearized Distance (cm)")
                % title(sprintf("%s nrn %d Epoch %d \n selIdx2 = %.2f, numSpikes = %d", ...
                %     S_nrn.area,S_nrn.ID,e,selIdx2,numel(ctrOutLSpkData)+numel(ctrOutRSpkData)))
                % if ~lgdExist
                %     legend("Right Trajs","Left Trajs","R-L",Location="best")
                %     lgdExist = 1;
                % end
                
            end

        end
        % pause
        % close all

        S_nrn.eChoiceSelectivityIdx = eSelIdx2;

        % Calculate LMRV value and add to nrn struct
        LMRV = calc_LMRV(C_behper{1,r}.eFracCorr,S_nrn.eFR);
        LMRV_thr = 0.5; % threshold for a cell to be learning modulated.
        S_nrn.LMRV = LMRV; % Assign LMRV to struct
        S_nrn.LMRVthreshold = LMRV_thr;
        if LMRV >= LMRV_thr
            % Determine if neuron fires more during wake or sleep
            % epochs
            if sum(S_nrn.eFR(restEpochs))/length(restEpochs) >= ...
                    sum(S_nrn.eFR(behEpochs))/length(behEpochs)
                S_nrn.LMRVtype = "rest-LMRV";
            elseif sum(S_nrn.eFR(restEpochs))/length(restEpochs) < ...
                    sum(S_nrn.eFR(behEpochs))/length(behEpochs)
                S_nrn.LMRVtype = "beh-LMRV";
            end
        else
            S_nrn.LMRVtype = "non-LMRV";
        end
        

        % Assigning state-related measures
        S_nrn.stateNames = C_combstates{1,r}{1,1}.stateNames; % Should be the same for all rats/epochs
        S_nrn.eStateFR = epochStateFR(nrn,:);


        nrninfo_MES{1,1}{nrn,1} = S_nrn; % Store with other cells

    end

    % NOTE: text must be char vectors, not string, or an "invalid file
    % identifier" error occurs.
    if any(~areaExist)
        fmt = ['%s cells with empty "area" field:',repmat(' %i',1,numel(find(~areaExist))),'\n'];
        fprintf(fmt,loadRats{r},find(~areaExist))
    end
    if any(~typeExist)
        fmt = ['%s cells with empty "type" field:',repmat(' %i',1,numel(find(~typeExist))),'\n'];
        fprintf(fmt,loadRats{r},find(~typeExist))
    end


    % Save data from each rat
    short_name = loadRats{r};
    chop_idx = strfind(loadRats{r},'_') - 1;
    if ~isempty(chop_idx)
        short_name = loadRats{r}(1:chop_idx); % Gets the first characters of the rat's name before an '_'.
        % So far this is only needed for ER1_NEW to remove the '_NEW'.
    end

    save(string(fullfile(dataDir,sprintf("%s_direct/%snrninfo_MES.mat",loadRats{r},short_name))),"nrninfo_MES", '-v7.3')
    fprintf("Saved file: %snrninfo_MES \n",short_name)
    
end

fprintf("Finished creating nrninfo_MES files. \n")

end

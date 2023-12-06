function [] = create_nrninfo_MES(dataDir,loadRats)
% Creates cellinfo_MES.mat structs that contain information pertaining to
% each neuron across all epochs.
%
% This function saves a file to each rat's folder called nrninfo_MES.mat.
% It is a cell array with dimensions 1 x (num rats), within which is a cell
% array (tot num neurons) x 1. Each cell contains a struct with that neuron's
% information.

filetypes = {'spikes01','behavperform'};

C_alldata = load_data(dataDir,loadRats,filetypes);

spikes_idx = find(contains(filetypes,'spikes01')); 
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end
behper_idx = find(contains(filetypes,'behavperform'));
if isempty(behper_idx)
    error("behavperform data must be loaded to run this analysis.")
end

C_allspikes = C_alldata(spikes_idx,:);
C_behper = C_alldata(behper_idx,:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;

%% Meat and potatoes



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
        end
        
        % Same thing for neuron type
        emptyIdxT = nrnType(nrn,:) == "";
        nonemptyTypes = nrnType(nrn, ~emptyIdxT);
        if ~isempty(nonemptyTypes) && all(nonemptyTypes == nonemptyTypes(1))
            S_nrn.type = nonemptyTypes(1);
        else
            S_nrn.type = "";
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
        



        nrninfo_MES{1,1}{nrn,1} = S_nrn; % Store with other cells

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
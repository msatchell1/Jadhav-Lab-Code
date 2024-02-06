% Script to identify trajectory-choice selective cells and analyze them
% Michael Satchell 11/02/23


%% Load data

clearvars;

dataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
loadRats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01',
% 'spikes01','tetinfo','linfields01','rippletime01','pos01','nrninfo_MES','combstates_MES'
filetypes = {'nrninfo_MES','combstates_MES'};

C_alldata = load_data(dataDir,loadRats,filetypes); % Load data to struct


C_nrninfo = C_alldata(find(contains(filetypes,'nrninfo_MES')),:);
C_combstates = C_alldata(find(contains(filetypes,'combstates_MES')),:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;
brainAreas = {'CA1','PFC'};


%% Plot choice selectivity across epochs

% There are four columns in ech epoch of linfields representing the four W-track
% trajectories:
% Col 1: out right, col 2: in right, col 3: out left, col 4: in left.
% Inside each trajectory lies 7 columns of data:
% Col 1 is the distance in cm along that linearized
% trajectory, col 2 occupancy, col 3 spike count, col 4 occupancy normalized
% firing rate, col 5 smoothed occupancy normalized firing rate, col 6
% smoothed occupancy, col 7 smoothed spike count.


cellType = "All";
area = "CA1";
% Combining data from all rats in the below arrays:
trajFR = cell(size(C_combstates{1,1})); % Firing rates on all trajectories for all epochs and cells
trajSC = cell(size(C_combstates{1,1})); % Spatial coverage
trajPk = cell(size(C_combstates{1,1})); % Peak firing rate
trajisPC = cell(size(C_combstates{1,1})); % Place cell (1) or not (0)
chcIdx = cell(size(C_combstates{1,1})); % Choice selectivity index.
% abs(choice selectivity index) averaged across all epochs for every neuron
allNrnChc = [];
for r = 1:size(C_nrninfo,2)
    eChcIdx = []; % Choice selectivity index.
    % chcIdxBehLMRV = [];
    for n = 1:size(C_nrninfo{1,r},1)

        nrn = C_nrninfo{1,r}{n,1};
        
        if strcmp(nrn.area,area) && strcmp(nrn.type,"Pyr") %&& strcmp(nrn.LMRVtype,cellType)
    
            for e = 1:size(nrn.eTrajFR,2) % Loops through all epochs
                if ~isempty(nrn.eTrajFR{1,e}) % assume that if this nrn has FR data 
                    % for this epoch, then all other traj measures will
                    % exist.
                    trajFR{1,e} = [trajFR{1,e}; nrn.eTrajFR{1,e}];
                    trajSC{1,e} = [trajSC{1,e}; nrn.eTrajCoverage{1,e}];
                    trajPk{1,e} = [trajPk{1,e}; nrn.eTrajFRPeak{1,e}];
                    trajisPC{1,e} = [trajisPC{1,e}; nrn.eTrajisPC{1,e}];
                    chcIdx{1,e} = [chcIdx{1,e}; nrn.eChoiceSelectivityIdx(1,e)];
                    eChcIdx = [eChcIdx; nrn.eChoiceSelectivityIdx];
                    
                    % if strcmp(nrn.LMRVtype,cellType)
                    %     chcIdxBehLMRV = [chcIdxBehLMRV; nrn.eChoiceSelectivityIdx];
                    % end
                end
            end

        end
    end
    allNrnChc = [allNrnChc; mean(abs(eChcIdx),2,'omitnan')];
    
    
    % % if ~isempty(chcIdxBehLMRV)
    %     figure;
    %     hold on
    %     plot(behEpochs,eChcIdx(:,behEpochs).',Color=[0.4,0.4,0.4])
    %     % plot(behEpochs,chcIdxBehLMRV(:,behEpochs).',Color='cyan')
    %     plot(behEpochs,mean(eChcIdx(:,behEpochs),1,'omitnan'),LineWidth=2,Color='k')
    %     title(sprintf("%s %s %s Choice Selectivity Across Epochs",loadRats{r},cellType,area))
    %     ylabel("Choice Selectivity Index (R - L)")
    %     xlabel("Epoch")
    % 
    % % end


end


% Analyze the distribution of cells based on selectivity
figure;
histogram(allNrnChc)

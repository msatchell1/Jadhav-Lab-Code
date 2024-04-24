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

%% Plot the occupancy normalized firing rate for each cell across epochs

% There are four columns in ech epoch of linfields representing the four W-track
% trajectories:
% Col 1: out right, col 2: in right, col 3: out left, col 4: in left.
% Inside each trajectory lies 7 columns of data:
% Col 1 is the distance in cm along that linearized
% trajectory, col 2 occupancy, col 3 spike count, col 4 occupancy normalized
% firing rate, col 5 smoothed occupancy normalized firing rate, col 6
% smoothed occupancy, col 7 smoothed spike count.

area = "PFC";
type = "Pyr";


for r = 1:size(C_nrninfo,2)

    for n = 1:size(C_nrninfo{1,r},1)
        
        nrn = C_nrninfo{1,r}{n,1};

        if strcmp(nrn.area,area) && strcmp(nrn.type,type)

            figure;
            tl = tiledlayout(2,4);
            lgdExist = 0;
            title(tl, sprintf("%s %s %s",loadRats{r},area,type))

            for e = behEpochs

                if ~isempty(nrn.eTrajLinField{e})

                    epochData = nrn.eTrajLinField{e};
                    spikeCountTrajs = [sum(epochData{1}(:,3)),sum(epochData{2}(:,3)),sum(epochData{3}(:,3)),sum(epochData{4}(:,3))];
                    
                    nexttile
                    hold on;
                    plot(epochData{1}(:,1), epochData{1}(:,5),'r')
                    plot(epochData{2}(:,1), epochData{2}(:,5),'r--')
                    plot(epochData{3}(:,1), epochData{3}(:,5),'b')
                    plot(epochData{4}(:,1), epochData{4}(:,5),'b--')

                    % plot(ctrRData(:,1),ctrRData(:,5)-ctrLData(:,5),Color=[0.7,0.7,0.7])
                    % yline(0,'--k')
                    % plot(linDist(ctrOutRSpkData(:,7)),rand(size(ctrOutRSpkData))*3-5,'.r')
                    % plot(linDist(ctrOutLSpkData(:,7)),rand(size(ctrOutLSpkData))*3-8,'.b')
                    xlabel("Linearized Distance (cm)")
                    ylabel("Occ Norm Firing Rate")
                    title(sprintf("nrn %d Epoch %d \n Spike Count %s",nrn.ID,e,num2str(spikeCountTrajs)))
                    if ~lgdExist
                        legend("out R","in R","out L","in L",Location="best")
                        lgdExist = 1;
                    end
                else % If no data for that epoch
                    nexttile
                end
            end
            pause
            close all
        end
    end
end

%% Plot choice selectivity across epochs

% There are four columns in ech epoch of linfields representing the four W-track
% trajectories:
% Col 1: out right, col 2: in right, col 3: out left, col 4: in left.
% Inside each trajectory lies 7 columns of data:
% Col 1 is the distance in cm along that linearized
% trajectory, col 2 occupancy, col 3 spike count, col 4 occupancy normalized
% firing rate, col 5 smoothed occupancy normalized firing rate, col 6
% smoothed occupancy, col 7 smoothed spike count.


cellType = "All"; % Can be commented out below.
area = "PFC";
% Combining data from all rats in the below arrays:
trajFR = cell(size(C_combstates{1,1})); % Firing rates on all trajectories for all epochs and cells
trajSC = cell(size(C_combstates{1,1})); % Spatial coverage
trajPk = cell(size(C_combstates{1,1})); % Peak firing rate
trajisPC = cell(size(C_combstates{1,1})); % Place cell (1) or not (0)
chcIdx = cell(size(C_combstates{1,1})); % Choice selectivity index.
% abs(choice selectivity index) averaged across all epochs for every neuron

for r = 1:size(C_nrninfo,2)
    
    for e = behEpochs % Loops through all epochs

        for n = 1:size(C_nrninfo{1,r},1)
        
            nrn = C_nrninfo{1,r}{n,1};
            
            if strcmp(nrn.area,area) && strcmp(nrn.type,"Pyr") %&& strcmp(nrn.LMRVtype,cellType) % Comment out to include all LMRV types
    
                trajFR{1,e} = [trajFR{1,e}; nrn.eTrajFR{1,e}];
                trajSC{1,e} = [trajSC{1,e}; nrn.eTrajCoverage{1,e}];
                trajPk{1,e} = [trajPk{1,e}; nrn.eTrajFRPeak{1,e}];
                trajisPC{1,e} = [trajisPC{1,e}; nrn.eTrajisPC{1,e}];
                chcIdx{1,e} = [chcIdx{1,e}; nrn.eChoiceSelectivityIdx(1,e)];

            else
                trajFR{1,e} = [trajFR{1,e}; NaN(1,4)];
                trajSC{1,e} = [trajSC{1,e}; NaN(1,4)];
                trajPk{1,e} = [trajPk{1,e}; NaN(1,4)];
                trajisPC{1,e} = [trajisPC{1,e}; NaN(1,4)];
                chcIdx{1,e} = [chcIdx{1,e}; NaN];
      
            end
        end
    end
    
    
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

selMat = [chcIdx{:}];
goodRows = ~all(isnan(selMat),2); % Rows with at least one non-NaN value
selMat = selMat(goodRows,:);

% % Analyze the distribution of cells based on selectivity
% figure;
% h = histogram(mean(selMat,2,'omitnan'),NumBins=round(2*sqrt(numel(mean(selMat,2,'omitnan')))));
% h.FaceColor=[0.4,0.4,0.4];
% title(sprintf("%s Cells %s Mean Selectivity Across Epochs",cellType,area))
% ylabel("Num Cells")
% xlabel(" Mean Selectivity Index")
% xlim([-2,2])

% Same plot but for abs of selectivity indices before averaging.
figure;
h = histogram(mean(abs(selMat),2,'omitnan'),NumBins=round(3*sqrt(numel(mean(selMat,2,'omitnan')))));
h.FaceColor="cyan";%[0.4,0.4,0.4];
title(sprintf("%s Cells %s Mean Selectivity Across Epochs",cellType,area))
ylabel("Num Cells")
xlabel(" Mean |Selectivity Index|")
xlim([0,2])


%% Compare choice selectivity to rate in subsequent sleep states for single epoch pairs


cellType = "All";
area = "CA1";

stateNames = C_combstates{1,1}{1,1}.stateNames;
stateColors = {[4,37,204]/204, [230,10,2]/230, [133,28,252]/252, [71,255,20]/255, [20,255,255]/255};
for s = 1:numel(stateNames)

    seFRs = cell(size(C_nrninfo)); % state and epoch FRs
    eSelIdx = cell(size(C_nrninfo));

    tl = tiledlayout('flow');
    ylabel(tl,sprintf("Subsequent %s Firing Rate (Hz)",stateNames{s}))
    xlabel(tl,"On-track Choice Selectivity Index")
    title(tl, sprintf("%s %s vs Selectivity For All Individual Epochs",area,stateNames{s}))
    

    for r = 1:size(C_nrninfo,2)

        nexttile
        hold on;
        title(sprintf("%s",loadRats{r}))
        xline(0,'--k')
    
        for n = 1:size(C_nrninfo{1,r},1)
    
            nrn = C_nrninfo{1,r}{n,1};
            
            if strcmp(nrn.area,area) && strcmp(nrn.type,"Pyr") %&& strcmp(nrn.LMRVtype,cellType)
                
                rates = cellfun(@(x) x(s), nrn.eStateFR)';
                idxs = nrn.eChoiceSelectivityIdx';
                seFRs{1,r} = [seFRs{1,r}; rates];
                eSelIdx{1,r} = [eSelIdx{1,r}; idxs];
    
                plot(idxs(behEpochs),rates(restEpochs(2:end)),'o',Color=stateColors{s})
                
            end
        end
 

    end


    pause

end


%% 
dtvals = -100 : 1 : 100;
ypos = exp(-abs(dtvals)/20);
yneg = -exp(-abs(dtvals)/20);
figure;
hold on
plot(dtvals(dtvals>=0),ypos(dtvals>=0),Color=[118, 9, 173]/173)
plot(dtvals(dtvals<=0),yneg(dtvals<=0),Color=[118, 9, 173]/173)

plot(dtvals(dtvals>=0),0.5*ypos(dtvals>=0),Color=[4, 110, 15]/110)
plot(dtvals(dtvals<=0),0.5*yneg(dtvals<=0),Color=[4, 110, 15]/110)

ylabel("Norm Synaptic Change")
xlabel("Spike time post - pre (ms)")
xline(0,'--k')
yline(0,'--k')
legend(["awake, N2"," ", "N3"])

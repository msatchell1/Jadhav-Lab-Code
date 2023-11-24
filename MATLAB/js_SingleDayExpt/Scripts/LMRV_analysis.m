% Script to calculate and then analyze neurons based on their LMRV index.
% Michael Satchell 11/02/23




%% Load data

clearvars;

dataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
loadRats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01',
% 'spikes01','tetinfo','linfields01','rippletime01','pos01','nrninfo_MES'
filetypes = {'nrninfo_MES','tetinfo','pos01','sws01','rem01','rippletimes_MES','behavperform'};

C_alldata = load_data(dataDir,loadRats,filetypes); % Load data to struct


% ------------------------------------------
% ------------------------------------------
% ------------------------------------------
% List of all states to sort the spike data into. Note that adding
% ripple state adds quite a bit of calculation time. 
% WARNING: This does not change the order of how these states are stored in
% C_allstates or later FR_allStates. Thus, stateFiles MUST MATCH ORDER
% THESE FILES ARE LISTED IN filetypes.
stateFiles = {'sws','rem','ripple'};
stateNames = {'SWS','REM','ripple','run','still'}; % Must match order states
% are added to C_allstates
% ------------------------------------------
% ------------------------------------------
% ------------------------------------------


sws_idx = find(contains(filetypes,'sws'));
if isempty(sws_idx)
    error("sws data must be loaded to run this analysis.")
end

rem_idx = find(contains(filetypes,'rem'));
if isempty(rem_idx)
    error("rem data must be loaded to run this analysis.")
end

% spikes_idx = find(contains(filetypes,'spikes01')); 
% if isempty(spikes_idx)
%     error("spikes01 data must be loaded to run this analysis.")
% end

% cellinfo_idx = find(contains(filetypes,'cellinfo'));
% if isempty(cellinfo_idx)
%     error("cellinfo data must be loaded to run this analysis.")
% end

states_idx = find(contains(filetypes, stateFiles));
if isempty(states_idx)
    error("at least one of the following rest state data must be loaded: \n +" + ...
        "   sleep01, waking01, sws01, rem01, or rippletime01.")
end

% linf_idx = find(contains(filetypes,'linfields'));
% if isempty(linf_idx)
%     error("linfields data must be loaded to run this analysis.")
% end

pos_idx = find(contains(filetypes,'pos01'));
if isempty(pos_idx)
    error("pos01 data must be loaded to run this analysis.")
end

rt_idx = find(contains(filetypes,'rippletime'));
if isempty(rt_idx)
    error("rippletime data must be loaded to run this analysis.")
end

ti_idx = find(contains(filetypes,'tetinfo'));
if isempty(ti_idx)
    error("tetinfo data must be loaded to run this analysis.")
end

behper_idx = find(contains(filetypes,'behavperform'));
if isempty(behper_idx)
    error("behavperform data must be loaded to run this analysis.")
end

nrninfo_idx = find(contains(filetypes,'nrninfo_MES'));
if isempty(nrninfo_idx)
    error("nrninfo_MES data must be loaded to run this analysis.")
end

C_swsstate = C_alldata(sws_idx,:);
C_remstate = C_alldata(rem_idx,:);
% C_allspikes = C_alldata(spikes_idx,:);
% C_allinfo = C_alldata(cellinfo_idx,:);
C_allstates = C_alldata(states_idx,:);
% C_alllinf = C_alldata(linf_idx,:);
C_runstate = create_runstate(C_alldata(pos_idx,:));
C_stillstate = create_stillstate(C_alldata(pos_idx,:));
C_allstates = [C_allstates; C_runstate; C_stillstate];
C_allriptimes = C_alldata(rt_idx,:);
C_alltetinfo = C_alldata(ti_idx,:);
C_behper = C_alldata(behper_idx,:);
C_nrninfo = C_alldata(nrninfo_idx,:);
C_pos = C_alldata(pos_idx,:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;


brainAreas = {'CA1','PFC'};




%% Assign LMRV value to all neurons
for r = 1:size(C_nrninfo,2)
    for nrn = 1:size(C_nrninfo{1,r},1)
        S_nrn = C_nrninfo{1,r}{nrn,1};
        LMRV = calc_LMRV(C_behper{1,r}.eFracCorr,S_nrn.eFR);
        C_nrninfo{1,r}{nrn,1}.LMRV = LMRV; % Assign LMRV to struct
    end
end

%% Plot LMRV with the FR of each PFC cell over time. I am doing this to get an idea of 
% how cells vary firing behavior across time. Also plot probability of
% correct choice

% LMRV_thr = 0.6; % threshold for plotting cells.
% 
% for r = 1:size(C_nrninfo,2)
% 
%     for nrn = 1:size(C_nrninfo{1,r},1)
%         S_nrn = C_nrninfo{1,r}{nrn,1};
% 
%         if strcmp(S_nrn.area,"PFC")
%             if S_nrn.LMRV >= LMRV_thr
%                 f = figure;
%                 title(sprintf("%s PFC Neuron %d \n %s | LMRV = %.2f", loadRats{r},S_nrn.ID,S_nrn.type,S_nrn.LMRV))
% 
%                 xlabel("Epoch")
%                 colororder({'b','k'})
% 
%                 yyaxis left
%                 hold on
%                 plot(1:17,S_nrn.eFR,'*-')
%                 plot(behEpochs,cellfun(@sum, S_nrn.eTrajisPC(behEpochs)),"cyan")
%                 ylabel("Mean Firing Rate (Hz)")
%                 yyaxis right
%                 plot(behEpochs,C_behper{1,r}.eFracCorr,'k')
%                 ylabel("Fraction of Correct Trials")
% 
%                 legend(["nrn FR","sum isPC","frac corr"],Location="Best")
%                 pause
%                 close all
%             end
%         end
% 
%     end
% 
% end


%% Get sorted occurances of all states

C_combStates = create_combStates(C_allstates,stateNames,dataDir);

%% Plot some example epochs

stateColors = [[0 0.4470 0.7410]; [1 0 0]; [0.5 0.1 1]; [0.6 0.9 0]; [0 0.9 1]];

r = 1;
for e = 1:size(C_combStates,2)

    velData = C_pos{1,r}{1,e}.data(:,5);
    timeData = C_pos{1,r}{1,e}.data(:,1);
    velData = sgolayfilt(velData, 2, 31); % Smooth velocity

    S_occ = C_combStates{r,e};
    occs = S_occ.combSortData; % Occurances of all states in temporal order.

    % Get PFC pyramidal neurons
    nrns = [];
    for n = 1:size(C_nrninfo{1,r},1)
        S_nrn = (C_nrninfo{1,r}{n,1});
        if strcmp(S_nrn.area,"PFC") && strcmp(S_nrn.type,"Pyr") && ...
                S_nrn.eHasSpikeData(1,e)
            nrns = [nrns; S_nrn];
        end
    end

    allSpikes = []; % Put all spikes into one vector
    for n = 1:size(nrns,1)
        allSpikes = [allSpikes; nrns(n).eSpikeData{e}(:,1)];
    end
    allSpikes = sort(allSpikes);
    
    winLen = 0.1; % window of time to average over for the firing rate
    FRtimes = nrns(1).eTimeRange{1,e}(1,1) : winLen : nrns(1).eTimeRange{1,e}(1,2);
    
    FRavg = [];
    for i = 1:length(FRtimes)
        FRavg(i) = sum( allSpikes >= (FRtimes(i)-winLen/2) & ...
            allSpikes <= (FRtimes(i)+winLen/2) )/(length(nrns)*winLen); % Counts the number of spikes
        % within winLen/2 on either side of each time point.
    end

    figure;
    hold on
    % plot(timeData,velData,'k')
    plot(FRtimes,FRavg,'k','DisplayName',"avg PFC Pyr FR")
    xlabel("Time (s)")
    ylabel("Firing Rate (Hz)")

    for i = 1:size(occs,1)
        x_vertices = [occs(i,1),occs(i,2),occs(i,2),occs(i,1)];
        y_vertices = [min(FRavg),min(FRavg),max(FRavg),max(FRavg)];
        patch(x_vertices, y_vertices, stateColors(occs(i,3),:),'FaceAlpha', 0.3, 'EdgeColor','none'...
            ,'HandleVisibility','off')
    
        % if i == 1
        %     runLabels{i} = "run times";
        % else
        %     runLabels{i} = "";
        % end
    end
    
    legend()
    
    pause
    close all
end

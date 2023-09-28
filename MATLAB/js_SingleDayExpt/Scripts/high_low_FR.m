% This script will split pyramidal cells based on their activity during
% behavior epochs into high and low FR groups for CA1 and PFC, then look
% for differences in how the firing rates of these groups changes across
% learning and during sleep.
% Michael Satchell 09/25/23

%% Load data

clearvars;

dataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
loadRats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01','spikes01','tetinfo','linfields01','rippletime01','pos01'
filetypes = {'spikes01','tetinfo','pos01','sws01','rem01','rippletimes_MES','behavperform'};

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
stateNames = {'sws','rem','ripple','run','still'}; % Must match order states
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

spikes_idx = find(contains(filetypes,'spikes01')); 
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end

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

C_swsstate = C_alldata(sws_idx,:);
C_remstate = C_alldata(rem_idx,:);
C_allspikes = C_alldata(spikes_idx,:);
% C_allinfo = C_alldata(cellinfo_idx,:);
C_allstates = C_alldata(states_idx,:);
% C_alllinf = C_alldata(linf_idx,:);
C_runstate = create_runstate(C_alldata(pos_idx,:));
C_stillstate = create_stillstate(C_alldata(pos_idx,:));
C_allstates = [C_allstates; C_runstate; C_stillstate];
C_allriptimes = C_alldata(rt_idx,:);
C_alltetinfo = C_alldata(ti_idx,:);
C_allbehper = C_alldata(behper_idx,:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;


brainAreas = {'CA1','PFC'};



%% Calculates firing rates, state occurance times, and state durations
[C_stateFRs,M_stateFR,M_stateOcc,M_stateDur] = calc_meanrates(brainAreas,C_allstates,C_allspikes);





%% Test how many cells switch between < 7 Hz FR and > 7 Hz over the epochs
% Do this to determine the best way to label interneurons.




%% First thing is to split the pyramidal cells into high and low FR.
% For now I will base high and low FR on all of epoch 2.

% for a = 1:size(M_stateFR,2)
%     for s = 1:size(M_stateFR,1)
% 
% 
% 
%     end
% end

figure;
hold on;
scatter(1:size(M_stateFR{4,1}(:,2),1),M_stateFR{4,1}(:,2))
title(sprintf("Run FRs for all Cells Epoch 2 CA1"))

figure;
hold on;
pyrFRrun = M_stateFR{4,1}(:,2);
pyrFRrun = pyrFRrun(pyrFRrun<7);
scatter(1:size(pyrFRrun,1),pyrFRrun)
title(sprintf("Run FRs for Pyramidal Cells Epoch 2 CA1"))
% It looks like a good way to split into low and high FR would be around 2
% Hz... But lets check PFC

figure;
hold on;
scatter(1:size(M_stateFR{4,2}(:,2),1),M_stateFR{4,2}(:,2))
title(sprintf("Run FRs for all Cells Epoch 2 PFC"))

figure;
hold on;
pyrFRrun = M_stateFR{4,2}(:,2);
pyrFRrun = pyrFRrun(pyrFRrun<7);
scatter(1:size(pyrFRrun,1),pyrFRrun)
title(sprintf("Run FRs for Pyramidal Cells Epoch 2 PFC"))
% 3 Hz seems like a better place to split PFC.

% M_highFR = {};






























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
C_allbehper = C_alldata(behper_idx,:);
C_nrninfo = C_alldata(nrninfo_idx,:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;


brainAreas = {'CA1','PFC'};




%% Show that my labeling worked properly by plotting spike width vs mean FR
% So one problem is that many cells only fire during wake or sleep epochs,
% and so have empty C_allspikes entries for the epochs that they do not
% fire (or were just not clustered on). 



for a = 1:size(brainAreas,2)
    figure;
    hold on;
    title(sprintf("%s Interneuron Labeling",brainAreas{a}))
    xlabel("Mean FR Across all Epochs (Hz)")
    ylabel("Spike Width (ms)")
    
    for r = 1:size(C_nrninfo,2)
    
        for nrn = 1:size(C_nrninfo{1,r},1)

            S_nrn = C_nrninfo{1,r}{nrn,1};

            if strcmp(S_nrn.area,brainAreas{a})
                if strcmp(S_nrn.type,"Pyr")
                    plot(S_nrn.meanFR,S_nrn.meanSW*1000,"ob")
                elseif strcmp(S_nrn.type,"Int")
                    plot(S_nrn.meanFR,S_nrn.meanSW*1000,"or")
                end
            end
    
        end
    end
end




% I should plot spike width across epochs for cells to make sure that it
% doesn't vary too much from epoch to epoch... otherwise I might not be
% recording the same cell.




%% Firing Rates. Plot across epochs


for a = 1:size(brainAreas,2)
    figure;
    hold on;
    title(sprintf("%s Mean FRs",brainAreas{a}))
    xlabel("Epoch")
    ylabel("Firing Rate (Hz)")
    
    for r = 1:size(C_nrninfo,2)
    
        for nrn = 1:size(C_nrninfo{1,r},1)

            S_nrn = C_nrninfo{1,r}{nrn,1};

            if strcmp(S_nrn.area,brainAreas{a})
                if strcmp(S_nrn.type,"Pyr")
                    plot(1:size(S_nrn.eFR,2),S_nrn.eFR,"ob")
                elseif strcmp(S_nrn.type,"Int")
                    plot(1:size(S_nrn.eFR,2),S_nrn.eFR,"or")
                end
            end
    
        end
    end
end


%% Plot histogram of firing rates during awake epochs to identify a good
% low/high FR dividing line.

% During all epochs
for a = 1:size(brainAreas,2)
    figure;
    hold on;
    title(sprintf("%s Mean Pyr FR, Mean of All Epochs",brainAreas{a}))
    xlabel("Firing Rate (Hz)")
    ylabel("Count")
    
    FRarray = [];
    for r = 1:size(C_nrninfo,2)
    
        for nrn = 1:size(C_nrninfo{1,r},1)

            S_nrn = C_nrninfo{1,r}{nrn,1};

            if strcmp(S_nrn.area,brainAreas{a})
                if strcmp(S_nrn.type,"Pyr")
                    FRarray = [FRarray, mean(S_nrn.eFR,2,'omitnan')];
                end
            end
    
        end
    end
    histogram(FRarray,100)
end

%% During behavioral epochs (epochs individually counted)
for a = 1:size(brainAreas,2)
    figure;
    hold on;
    title(sprintf("%s Mean Pyr FR Individual Behavioral Epochs",brainAreas{a}))
    xlabel("Firing Rate (Hz)")
    ylabel("Count")
    
    FRarray = [];
    for r = 1:size(C_nrninfo,2)
    
        for nrn = 1:size(C_nrninfo{1,r},1)

            S_nrn = C_nrninfo{1,r}{nrn,1};

            if strcmp(S_nrn.area,brainAreas{a})
                if strcmp(S_nrn.type,"Pyr")
                    FRarray = [FRarray, S_nrn.eFR(behEpochs)];
                end
            end
    
        end
    end
    histogram(FRarray,5*ceil(sqrt(size(FRarray,2))))
end


% During rest epochs (epochs individually counted)
for a = 1:size(brainAreas,2)
    figure;
    hold on;
    title(sprintf("%s Mean Pyr FR Individual Rest Epochs",brainAreas{a}))
    xlabel("Firing Rate (Hz)")
    ylabel("Count")
    
    FRarray = [];
    for r = 1:size(C_nrninfo,2)
    
        for nrn = 1:size(C_nrninfo{1,r},1)

            S_nrn = C_nrninfo{1,r}{nrn,1};

            if strcmp(S_nrn.area,brainAreas{a})
                if strcmp(S_nrn.type,"Pyr")
                    FRarray = [FRarray, S_nrn.eFR(restEpochs)];
                end
            end
    
        end
    end
    histogram(FRarray,5*ceil(sqrt(size(FRarray,2))))
end



%% Plot each epoch's FRs separately
for e = 1:17

    tl = tiledlayout(2,1);
    tl.Title.String = sprintf("Mean Pyr FR Epoch %d",e);
    xlabel(tl, "Firing Rate (Hz)")
    ylabel(tl, "Count")

    for a = 1:size(brainAreas,2)

    
        FRarray = [];
        for r = 1:size(C_nrninfo,2)
        
            for nrn = 1:size(C_nrninfo{1,r},1)
    
                S_nrn = C_nrninfo{1,r}{nrn,1};
    
                if strcmp(S_nrn.area,brainAreas{a})
                    if strcmp(S_nrn.type,"Pyr")
                        FRarray = [FRarray, S_nrn.eFR(e)];
                    end
                end
        
            end
        end
        nexttile        
        histogram(FRarray,5*ceil(sqrt(size(FRarray,2))))
        if strcmp(brainAreas{a},"PFC")
            xline(3,"--k")
        elseif strcmp(brainAreas{a},"CA1")
            xline(1,"--k")
        end
        title(brainAreas{a})
    end
    pause
    close all
end

%% Plot each rat's FRs separately (individually counted epochs) in a 4-tiled
% plot of behavioral and rest epochs in CA1 and PFC
brNames = {"Behavior","Rest"};
brEpochs = {behEpochs, restEpochs};
for r = 1:size(C_nrninfo,2)
    tl = tiledlayout(2,2);
    tl.Title.String = sprintf("Mean Pyr FR Individual Epochs Rat %s",loadRats{r});
    xlabel(tl, "Firing Rate (Hz)")
    ylabel(tl, "Count")
    
    axes = [];
    for a = 1:size(brainAreas,2)
        for br = 1:size(brNames,2)
        
            FRarray = [];
           
            for nrn = 1:size(C_nrninfo{1,r},1)
        
                S_nrn = C_nrninfo{1,r}{nrn,1};
        
                if strcmp(S_nrn.area,brainAreas{a})
                    if strcmp(S_nrn.type,"Pyr")
                        FRarray = [FRarray, S_nrn.eFR(brEpochs{1,br})];
                    end
                end
        
            end
            axes(end+1) = nexttile;        
            histogram(FRarray,BinWidth=0.15)%size(behEpochs,2)*ceil(sqrt(size(C_nrninfo{1,r},1))))
            % if strcmp(brainAreas{a},"PFC")
            %     xline(3,"--k")
            % elseif strcmp(brainAreas{a},"CA1")
            %     xline(1,"--k")
            % end
            title(sprintf("%s | %s",brNames{br},brainAreas{a}))
        end
    end
    linkaxes(axes,'x')
    pause
    close all
end

%% Plot the probability of choosing the correct arm for each rat over all trials

for r = 1:size(C_allbehper,2)
    
    figure;
    hold on;
    title(sprintf("Outbound Trial Performance %s",loadRats{r}))
    ylabel("Probability of Correct Choice")
    xlabel("Trial")
    plot(C_allbehper{1,r}.outprobcorrect(:,1),'k')
    pause
    close all


end







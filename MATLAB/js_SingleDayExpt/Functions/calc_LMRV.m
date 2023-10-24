function [LMRV] = calc_LMRV(ePerf,eFR)
%CALC_LMRV Calculate the Learning-Modulated Rate Variance index (LMRV) for a neuron.
%   The LMRV measures how the firing rate between subsequent rest and
%   behavioral epochs varies as a function of increases in task
%   performance. The first and last sleep epochs must be left out of this
%   analysis as there is no value of change in performance that can be
%   associated with them.
%
% Inputs:
%   ePerf - row vector of performance across behavioral epochs
%   eFR - row vector of firing rates across all epochs. The first and last
%   values will not be considered if they fall outside the range of
%   the behavioral epochs. eFR cannot have NaNs.
%
% Outputs:
%   LMRV - LMRV index


% eFR = [1.91170631902480	0.460790767733273	2.16597534192263	0.567374705559981	2.31557429673496	0.520032482028878	2.59307070812742	0.385433084314785	2.37179643736526	0.957470127453648	2.69534326785954	0.756453591186236	2.20775870606551	0.580041257236653	2.41107012655311	0.522692611983724	3.13850014769339];
% ePerf = [0.481481481481481	0.769230769230769	0.866666666666667	0.870370370370370	0.933333333333333	0.806451612903226	0.842105263157895	0.758620689655172];


eFR = eFR(2:end-1); % Removes first and last epoch.

if any(isnan(eFR))
    fprintf("Warning: NaN values exist in firing rate vector. LMRV cannot be calculated. \n")
    LMRV = NaN;
else
        
    
    % % First we need to interpolate to get values of performance during the rest
    % % epochs. This will correspond to the learning achieved in the previous
    % % beh epoch. This puts NaNs into intrPerf at the edges of ePerf.
    % intrPerf = interp1(2:2:2*length(ePerf),ePerf,1:length(eFR),'linear');
    
    % Get slope of performance. The slope represents the amount "learned"
    % in that behavioral and subsequent rest epoch, so those are the epochs
    % that slope should be associated with.
    eSlope = diff(ePerf);
    % Interpolate to get slope values for the rest periods. The slope in
    % the previous behavioral epoch is used.
    intrSlope = interp1(2:2:2*length(eSlope),eSlope,2:length(eFR),'previous','extrap');

    % Getting change in rates
    eCR = diff(eFR);
    % Squared change in rates
    eSCR = eCR.^2;
    
    Smax = max(intrSlope);
    LMRVsum = 0;
    % Calculate numerator of LMRV
    for i = 1:length(eCR)
        S = intrSlope(i); % Change in performance at epoch
        if S < 0 % Disallow negative slopes because they show an equal lack
            % of learning as slope = 0.
            S = 0;
        end

        dr2 = eSCR(i); % Change in firing rate squared at epoch
        temp = dr2 * ((S/Smax)^0.5 + (S-Smax)/(S+Smax)); 
        LMRVsum = LMRVsum + temp;
    end

    LMRV = LMRVsum/(sum(eSCR)); % Normalize by the sum of squared differences in rates

     

    % figure;
    % hold on;
    % % plot(2:2:2*length(eSlope)+1,eSlope)
    % plot(2:length(intrSlope)+1,intrSlope)
    % plot(2:2:2*length(ePerf), ePerf)
    % plot(2:length(eSCR)+1,eSCR/max(eSCR))
    % % plot(2:length(eFR)+1, eFR)
    
    

end

end


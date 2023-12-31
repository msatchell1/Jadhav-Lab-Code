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


% eFR = eFR(2:end-1); % Removes first and last epoch.
eFR = eFR/max(eFR); % Normalize to 1 to eliminate differences between high and low FR neurons.

if any(isnan(eFR))
    %fprintf("Warning: NaN values exist in firing rate vector. LMRV cannot be calculated. \n")
    LMRV = NaN;
else

    % First we need to interpolate to get values of performance during the rest
    % epochs. This will correspond to the learning achieved in the previous
    % beh epoch. This puts NaNs into intrpPerf at the edges of ePerf. This
    % array starts at epoch 2 and ends at epoch 16, so should have 15
    % elements.
    intrpPerf = interp1(2:2:2*length(ePerf),ePerf,2:2*length(ePerf),'linear');
    % plot(2:length(intrpPerf)+1,intrpPerf)

    % Get slope of performance. The slope represents the amount "learned"
    % in that behavioral and subsequent rest epoch, so those are the epochs
    % that slope should be associated with.
    ePSlope = diff(intrpPerf);
    % % Interpolate to get slope values for the rest periods. The slope in
    % % the previous behavioral epoch is used.
    % intrSlope = interp1(2:2:2*length(eSlope),eSlope,2:length(eFR),'previous','extrap');
    
    % Getting change in rates
    eCR = diff(eFR);

    eCR(eCR<0) = 0;

    % The absolute value in the change of rates is what matters
    % eCR = abs(eCR);
    
    % Also, normalize change in rates to 1 to prevent high-FR and low-FR
    % neurons from being different.
    % eCR = eCR/max(eCR);

    % % And the sum of adjacent rates to normalize the difference by
    % eSR = eFR(1:end-1) + eFR(2:end);

    % To make the dp vector equal in length to the dr vector, simply extend
    % the end of the dp vector using the end values.
    ePSlope = [ePSlope(1,1), ePSlope, ePSlope(1,end)];

    % Negative slopes have no meaning, so these are set to zero to avoid
    % the correlation being affected by them.
    ePSlope(ePSlope<0) = 0;

    % Now normalize slopes to 1 to make the LMRV index simpler.
    ePSlope = ePSlope/max(ePSlope);
    
    
    
    Smax = max(ePSlope);
    LMRVsum = 0;
    % Calculate numerator of LMRV
    for i = 1:length(eCR)
        S = ePSlope(i); % Change in performance at epoch. Is greater than or equal to zero.
        dr = eCR(i); % Change in firing rate at epoch.

        temp = dr^2 * ((S/Smax)^0.1 + ((S-Smax)/(S+Smax))); 
        LMRVsum = LMRVsum + temp;
    end

    LMRV = LMRVsum/(sum(eCR.^2)); % Normalize to a max of 1.

     

    % figure;
    % 
    % tl = tiledlayout(1,2);
    % nexttile
    % plot(eFR)
    % hold on
    % plot(1.5:length(eCR)+0.5,eCR)
    % legend("eFR","eCR >= 0",Location="best")
    % title("Firing Rate")
    % 
    % nexttile
    % plot(2:2:2*length(ePerf),ePerf)
    % hold on
    % plot(1.5:length(ePSlope)+0.5,ePSlope)
    % legend("ePerf","ePSlope >= 0",Location="best")
    % title("Performance")
    % 
    % title(tl,sprintf("LMRV = %.2f",LMRV))

    % figure;
    % plot(1.5:length(eCR)+0.5,eCR/max(eCR))
    % hold on
    % plot(1.5:length(ePSlope)+0.5,ePSlope/max(ePSlope))
    % legend("dr","dp")
    % title("dr and dp Normalized to 1")
    
    

end

end


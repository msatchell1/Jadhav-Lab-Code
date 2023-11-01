function [corrVal] = dr_dp_corr(ePerf,eFR)
%DR_DP_CORR Returns correlation between change in FR dr and change in
%performance dp.
%   This measure does not seem to work as well as LMRV for quantifying how
%   much the change in firing rate correlates with the change in
%   performance.


% eFR = [1.91170631902480	0.460790767733273	2.16597534192263	0.567374705559981	2.31557429673496	0.520032482028878	2.59307070812742	0.385433084314785	2.37179643736526	0.957470127453648	2.69534326785954	0.756453591186236	2.20775870606551	0.580041257236653	2.41107012655311	0.522692611983724	3.13850014769339];
% ePerf = [0.481481481481481	0.769230769230769	0.866666666666667	0.870370370370370	0.933333333333333	0.806451612903226	0.842105263157895	0.758620689655172];

% Note, the first and last epochs are removed because the change in firing
% rates between these epochs and their adjacent epochs cannot be compared
% with a change in performance. This is because the performance vector is shorter due to
% sleep being the first and last epochs.
eFR = eFR(2:end-1); % Removes first and last epoch.


if any(isnan(eFR))
    %fprintf("Warning: NaN values exist in firing rate vector. Correlation cannot be calculated. \n")
    corrVal = NaN;
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

    % Negative slopes have no meaning, so these are set to zero to avoid
    % the correlation being affected by them.
    ePSlope(ePSlope<0) = 0;


    corrMat = corr([abs(eCR)', ePSlope']);
    corrVal = corrMat(1,2); % Extract an off-diagonal element of the correlation
    % matrix. This is the cross-correlation between the change in FR and change
    % in performance.

    figure;
    
    tl = tiledlayout(1,2);
    nexttile
    plot(2:length(eFR)+1,eFR)
    hold on
    plot(2.5:length(eCR)+1.5,abs(eCR))
    legend("eFR","abs(eCR)",Location="best")
    title("Firing Rate")

    nexttile
    plot(2:2:2*length(ePerf),ePerf)
    hold on
    plot(2.5:length(ePSlope)+1.5,ePSlope)
    legend("ePerf","ePSlope >= 0",Location="best")
    title("Performance")
    
    title(tl,sprintf("corr = %.2f",corrVal))
    

end
end


function [initialRate, pts_used] = compInitialRate(stoichCoeff, time_pts, concs, rxn_percent)
    % Computes the initial rate of reaction.
    % time_pts(1) should be 0.
    % concs(1) should be the initial concentration.
    % concs and time_pts are both n by 1 arrays of nonnegative floats.
    if nargin == 3
        rxn_percent = 10;
    end
    init_conc = concs(1);
    total_conc_variation = max(concs)-min(concs);
    restricted_variation = rxn_percent/100*total_conc_variation;
    conc_diffs = concs(2:end) - init_conc;
    conc_variations = abs(conc_diffs);
    restricted_diffs = conc_diffs(conc_variations <= restricted_variation);
    pts_used = size(restricted_diffs, 1);
    
    if isempty(restricted_diffs)
        initialRate = conc_variations(1)/(time_pts(2)*stoichCoeff);
    else
        A = [ones(pts_used, 1) time_pts(2:(pts_used+1))];
        b = restricted_diffs./time_pts(2:(pts_used+1));
        x = A\b;
        initialRate = abs(x(1))/stoichCoeff;
    end
end

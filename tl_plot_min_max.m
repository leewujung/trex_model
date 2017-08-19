function [tlmin,tlmax] = tl_plot_min_max(p)    
    p_abs = abs(p);
    icount = find(p_abs>1e-7);         % for stats, only these values count <-- taken from Kraken
    p_abs(p_abs<1e-7) = 1e-7;
    p_abs(isnan(p_abs)) = 1e-6;
    p_abs(isinf(p_abs)) = 1e-6;
    
    tl = -20*log10(p_abs);
    
    % Get color scale for full field
    tlmed = median(tl(icount));	% median value
    tlstd = std(tl(icount));     % standard deviation
    tlmax = tlmed+0.75*tlstd;       % max for colorbar
    tlmax = 10*round(tlmax/10);     % make sure the limits are round numbers
    tlmin = tlmax-50;               % min for colorbar

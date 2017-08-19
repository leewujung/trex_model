function [tlt,tlmax,tlmin] = clean_pressure_data(pressure)

pressure = squeeze(pressure(1,1,:,:));

tlt = abs( pressure );

tlt( isnan( tlt ) ) = 1e-6;   % remove NaNs
tlt( isinf( tlt ) ) = 1e-6;   % remove infinities

icount = find( tlt > 1e-7 );         % for stats, only these values count
tlt( tlt < 1e-7 ) = 1e-7;            % remove zeros
tlt = -20.0 * log10( tlt );          % so there's no error when we take the log

% compute some statistics to automatically set the color bar

tlmed = median( tlt( icount ) );    % median value
tlstd = std( tlt( icount ) );       % standard deviation
tlmax = tlmed + 0.75 * tlstd;       % max for colorbar
tlmax = 10 * round( tlmax / 10 );   % make sure the limits are round numbers
tlmin = tlmax - 50;                 % min for colorbar
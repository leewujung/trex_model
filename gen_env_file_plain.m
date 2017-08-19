% 2016 03 30   Use write_env.m to generate environment file
% function gen_env_file(env_param,freq)
% Function for generating ENV file for Kraken
%

env_fname = 'test_env';
env_title = 'test gen profile';
freq = 3450;
nmesh_delta = 0.01;  % nmesh resolution

ssp{1}.z = [0, 19.0];
ssp{1}.cp = [1525, 1525];
ssp{1}.rho = 1;
ssp{1}.sigma = 0;
ssp{1}.attn_p = 0;

ssp{2}.z = [19.0, 114.0];
ssp{2}.cp = [1630, 1630];
ssp{2}.rho = 2;
ssp{2}.sigma = 0;
ssp{2}.attn_p = 0.545750541536736;   % attenuation [dB/lambda]

rootFinderSwitch = 1;  % 0-slow, 1-fast
top_opt = 'NVW ';
bot_opt = 'V';
bot_sigma = 0;
clim = [0, 3000];  % search limit of phase speed [m/s]
rmax = 0;
src_z = 17.8;
rcv_z_range = [0 ssp{end}.z(end)];
rcv_z_delta = 0.1;
rcv_z_num = rcv_z_range(2)/rcv_z_delta+1;

if rootFinderSwitch==1
    top_opt = [top_opt, '.'];
end

% Start writing into file
if (strcmp(env_fname,'ENVFIL')==0 && ~strcmp(env_fname(end-3:end), '.env'))
    env_fname = [env_fname '.env']; % append extension
end
fid = fopen(env_fname, 'w');   % create new envfil

fprintf(fid,'''%s'' ! TITLE \r\n', env_title );
fprintf(fid,'%d  \t \t \t ! FREQ (Hz) \r\n', freq );
fprintf(fid,'%d    \t \t \t ! NMEDIA \r\n', length(ssp) );
fprintf(fid,'''%s'' \t \t \t ! TOP OPTIONS \r\n', top_opt);

% SSP
for iS = 1:length(ssp)
    nmesh = ssp{iS}.z(end)/nmesh_delta+1;  % number of mesh points
    if ~isfield(ssp{iS},'cs')
        ssp{iS}.cs = zeros(length(ssp{iS}.cp));
    end
    if ~isfield(ssp{iS},'attn_s')
        ssp{iS}.attn_s = zeros(length(ssp{iS}.cp));
    end
    
    fprintf( fid, '%d %4.1f %6.1f \t ! MEDIUM INFO [NMESH SIGMA Z(NSSP)] \r\n', nmesh, ssp{iS}.sigma, ssp{iS}.z(end));
    
    % first line
    fprintf(fid, '\t %6.1f %6.1f %6.1f %6.1f %6.15f %6.1f \t ! SOUND SPEED PROFILE [Z CP CS RHO AP AS] \r\n', ...
                    ssp{iS}.z(1), ssp{iS}.cp(1), ssp{iS}.cs(1), ssp{iS}.rho(1), ssp{iS}.attn_p(1), ssp{iS}.attn_s(1));

    for ii = 2:length(ssp{iS}.cp)
        fprintf( fid, '\t %6.1f %6.1f / \r\n', ...
            ssp{iS}.z(ii), ssp{iS}.cp(ii));
    end
end

% Lower halfspace
fprintf(fid,'''%s'' %6.1f  \t \t ! BOTTOM BOUNDARY CONDITION  [TYPE SIGMA] \r\n', bot_opt, bot_sigma);

% Phase speed limit
fprintf(fid, '%6.1f %6.1f \t \t ! CMIN CMAX (m/s)  PHASE SPEED LIMIT \r\n', clim(1), clim(2));

% Range max
fprintf(fid, '%6.1f \t \t ! RMAX (km) \r\n', rmax);

% Source
fprintf(fid, '%d \t \t ! NSD \r\n', length(src_z));
fprintf(fid, '%6.1f \t \t ! SD(1:NSD) (m) \r\n', src_z);

% Receiver
fprintf(fid, '%d \t \t ! NRD \r\n', rcv_z_num);
fprintf(fid, '%6.1f %6.1f / \t \t ! RD(1:NSD) (m) \r\n', rcv_z_range(1), rcv_z_range(2));

fclose( fid );





function gen_env_file(freq,env_param,env_file,env_title,root_finder)
% Generate ENV file for Kraken
%
% INPUT
%   freq                frequency [Hz]
%   env_param           env parameters such as env_param.SSP
%   env_file            file name/path of the ENV file to be created
%   env_title           title of the ENV file
%   root_finder         root finding scheme
%                       0-fast but can drop lots of important modes
%                       1-slow but capture most modes

if root_finder==1
    env_param.top_opt = [env_param.top_opt, '.'];
end

% Start writing into file
if (strcmp(env_file,'ENVFIL')==0 && ~strcmp(env_file(end-3:end), '.env'))
    env_file = [env_file '.env']; % append extension
end
fid = fopen(env_file, 'w');   % create new envfil

fprintf(fid,'''%s'' ! TITLE \r\n', env_title );
fprintf(fid,'%d  \t \t \t ! FREQ (Hz) \r\n', freq );
fprintf(fid,'%d    \t \t \t ! NMEDIA \r\n', length(env_param.SSP) );
fprintf(fid,'''%s'' \t \t \t ! TOP OPTIONS \r\n', env_param.top_opt);

% env_param.SSP
for iS = 1:length(env_param.SSP)
    if iS==1
        nmesh = env_param.SSP{iS}.z(end)/env_param.nmesh_delta+1;  % number of mesh points
    else
        nmesh = (env_param.SSP{iS}.z(end)-env_param.SSP{iS-1}.z(end))/env_param.nmesh_delta+1;  % number of mesh points
    end
    if ~isfield(env_param.SSP{iS},'cs')
        env_param.SSP{iS}.cs = zeros(length(env_param.SSP{iS}.cp));
    end
    if ~isfield(env_param.SSP{iS},'attn_s')
        env_param.SSP{iS}.attn_s = zeros(length(env_param.SSP{iS}.cp));
    end
    
    fprintf( fid, '%d %4.1f %6.1f \t ! MEDIUM INFO [NMESH SIGMA Z(Nenv_param.SSP)] \r\n', nmesh, env_param.SSP{iS}.sigma, env_param.SSP{iS}.z(end));
    
    % first line
    fprintf(fid, '\t %6.1f %6.1f %6.1f %6.1f %6.15f %6.1f \t ! SOUND SPEED PROFILE [Z CP CS RHO AP AS] \r\n', ...
                    env_param.SSP{iS}.z(1), env_param.SSP{iS}.cp(1), env_param.SSP{iS}.cs(1), env_param.SSP{iS}.rho(1), env_param.SSP{iS}.attn_p(1), env_param.SSP{iS}.attn_s(1));

    for ii = 2:length(env_param.SSP{iS}.cp)
        fprintf( fid, '\t %6.1f %6.1f / \r\n', ...
            env_param.SSP{iS}.z(ii), env_param.SSP{iS}.cp(ii));
    end
end

% Lower halfspace
fprintf(fid,'''%s'' %6.1f  \t \t ! BOTTOM BOUNDARY CONDITION  [TYPE SIGMA] \r\n', env_param.bot_opt, env_param.bot_sigma);

% Phase speed limit
fprintf(fid, '%6.1f %6.1f \t \t ! CMIN CMAX (m/s)  PHASE SPEED LIMIT \r\n', env_param.clim(1), env_param.clim(2));

% Range max
fprintf(fid, '%6.1f \t \t ! RMAX (km) \r\n', env_param.rmax);

% Source
fprintf(fid, '%d \t \t ! NSD \r\n', length(env_param.src_z));
fprintf(fid, '%6.1f \t \t ! SD(1:NSD) (m) \r\n', env_param.src_z);

% Receiver
fprintf(fid, '%d \t \t ! NRD \r\n', env_param.rcv_z_num);
fprintf(fid, '%6.1f %6.1f / \t \t ! RD(1:NSD) (m) \r\n', env_param.rcv_z_range(1), env_param.rcv_z_range(2));

fclose( fid );





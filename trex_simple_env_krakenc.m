% 2017 08 17  Run KrakenC for simple trex environment

addpath 'F:\Dropbox\0_APL_normal_mode\kraken\Matlab\wjlee_codes'

clear

freq_all = 1800:200:3600;  % [Hz]
root_finder = 1;
param_path = pwd;
param_file = 'env_param_trex13.mat';
env_path = pwd;
env_file_pre = 'trex13_env';
env_title = 'TREX13 water depth 0-19m, sediment 19-57m';

kraken_exe = which( 'krakenc.exe' );
env_param = load(fullfile(param_path,param_file));

tic
parfor iF=1:length(freq_all)
    freq = freq_all(iF);
    disp(['freq = ',num2str(freq),' Hz']);
    
    % generate ENV file
    env_file = sprintf('%s_%06.1fHz',env_file_pre,freq);
    gen_env_file(freq,env_param,fullfile(env_path,env_file),env_title,root_finder);
    
    % run Kraken
    cmdstr = sprintf('%s %s',kraken_exe,env_file);
    system(cmdstr);
end
toc


% 2016 03 31  Multi-freq kraken runs for time domain prediction

clear

usrn = getenv('username');

if strcmp(usrn,'Wu-Jung')   % APL computer name
    addpath(genpath('F:\Dropbox\0_APL_normal_mode\kraken'));
    env_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160516_mode_calc_100m';
    param_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160516_mode_calc_100m';
else
    addpath(genpath(['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken']));
    env_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160516_mode_calc_100m'];
    param_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160516_mode_calc_100m'];
end

if ~exist(env_path,'dir')
    mkdir(env_path);
end

freq_all = [0.1:0.1:0.9,1.1:0.2:50];  % [Hz]
root_finder = 1;
param_file = 'env_param_100m.mat';
env_file_pre = 'mfenv';
env_title = 'Multi-freq 100m env';

ori_path = pwd;
cd(env_path);

kraken_exe = which( 'krakenc.exe' );
env_param = load(fullfile(param_path,param_file));

tic
parfor iF=1:length(freq_all)
    freq = freq_all(iF);
    disp(['freq = ',num2str(freq),' Hz']);
    
    % generate ENV file
    env_file = sprintf('%s_%06.1f',env_file_pre,freq);
    gen_env_file(freq,env_param,fullfile(env_path,env_file),env_title,root_finder);
    
    % run Kraken
    cmdstr = sprintf('%s %s',kraken_exe,env_file);
    system(cmdstr);
end
toc


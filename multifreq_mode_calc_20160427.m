% 2016 03 31  Multi-freq kraken runs for time domain prediction

clear 
usrn = getenv('username');

% Set up various paths
kraken_path = 'F:\Dropbox\0_APL_normal_mode\kraken\Matlab';
path_cell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
    kraken_on_path = any(strcmpi(kraken_path, path_cell));
else
    kraken_on_path = any(strcmp(kraken_path, path_cell));
end

if strcmp(usrn,'Wu-Jung')   % APL computer name
    if ~kraken_on_path   % add path if Kraken is not already on path
        addpath(genpath('F:\Dropbox\0_APL_normal_mode\kraken'));
    end
    base_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
else
    if ~kraken_on_path   % add path if Kraken is not already on path
        addpath(genpath(['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken']));
    end
    base_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests'];
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
env_path = fullfile(base_path,script_name);
if ~exist(env_path,'dir')
    mkdir(env_path);
end


freq_all = 4:4:10000;  % [Hz]
root_finder = 1;
param_path = env_path;
param_file = 'env_param_trex13.mat';
env_file_pre = 'mfenv';
env_title = 'Multi-freq TREX13 env';

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


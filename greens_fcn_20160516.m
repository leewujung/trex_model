% 2016 05 16  Get Green's function for f-k diagram

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
        addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    end
    base_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
    mode_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160516_mode_calc_100m';
    param_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160516_mode_calc_100m';
else
    if ~kraken_on_path   % add path if Kraken is not already on path
        addpath(genpath(['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken']));
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    end
    base_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests'];
    mode_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160516_mode_calc_100m'];
    param_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160516_mode_calc_100m'];
end
env_file_pre = 'mfenv';

% Modeling parameters
freq_step = 0.1;
freq_all = freq_step:freq_step:50;  % [Hz]
param_file = 'env_param_100m.mat';
env_param = load(fullfile(param_path,param_file));

z0 = 17.8;  % source depth [m]
z = 15;   % receiver depth [m]

% Store parameters
param.env_param_path = param_path;
param.env_param_file = param_file;
param.env_param = env_param;
param.mode_path = mode_path;

param.z0 = z0;
param.z = z;
param.freq_step = freq_step;
param.freq_max = freq_all(end);

% Save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

cw = param.env_param.SSP{1}.cp(1);
kr = 0:0.001:0.2;

% Assemble across frequency
g = nan(length(freq_all),length(kr));
for iF=1:length(freq_all)
    freq = freq_all(iF);
    disp(['freq = ',num2str(freq),' Hz']);
    %     kmax = 2*pi*freq/cw;
    
    % load modes
    env_file = sprintf('%s_%06.1f',env_file_pre,freq);
    mod_file = [env_file,'.mod'];
    modes = read_modes(fullfile(mode_path,mod_file));
    mode_incl = length(modes.k);
    
    if mode_incl~=0
        phi_incl = modes.phi(:,1:mode_incl);
        phi_z = interp1(modes.z,phi_incl,z);
        phi_z0 = interp1(modes.z,phi_incl,z0);
        krm = modes.k;
        
        gm = nan(length(modes.k),length(kr));
        for iM=1:length(modes.k)
            gm(iM,:) = phi_z(iM).*phi_z0(iM)./(kr.^2-krm(iM).^2);
        end
        g(iF,:) = nansum(gm,1)/(2*pi);
    end
end



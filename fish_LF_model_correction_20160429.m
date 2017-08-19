% 2016 04 29  Fish LF scattering, comparison of new and old code
%             old code has an error on how the radiation damping is
%             incorporated with the viscous and thermal damping

clear
usrn = getenv('username');

% Set up various paths
fish_scat_model_path = 'F:\Dropbox\0_CODE\fish_scattering\fish_scat_model_new';
path_cell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
    fish_on_path = any(strcmpi(fish_scat_model_path, path_cell));
else
    fish_on_path = any(strcmp(fish_scat_model_path, path_cell));
end

if strcmp(usrn,'Wu-Jung')   % APL computer name
    if ~fish_on_path
        addpath('F:\Dropbox\0_CODE\fish_scattering\fish_scat_model_new');
    end
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    base_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
else
    if ~fish_on_path
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\fish_scattering\fish_scat_model_new']);
    end
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    base_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests'];
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

param.gamma_alpha = 1.4;  % Nero et al. eq(4)
param.rho = 1050;         % from Chu's code
param.cw = 1500;
param.c_a = 340;
param.tau = 200;          % [N/m] surface tension
param.kappa_a = 5.5e-3;   % thermal conductivity (cal/cm s ^oC)

% 1dbar = *1e5-> dynes/cm^2; *1e4-> Pascals
param.P0 = 1.013e5;      % Pa (1 atm at sea surface)

%   L         fish length in [m]
%   shape     'prosph' 1 - prolate spheroid
%             'swb'    2 - swimbladder
%             'cyl'    3 - cylinder
%   opt       1 - match V00 for 25 cm fish
%             2 - use default radius or swimbladder value
%   V00       swimbladder volume at surface [cm^3]

L = 0.05;
shape = 'prosph';
V00 = 5;
D = 10;
opt = 1;
xi = 5;
freq = 1000:10:10000;

[swb_shape.xg,swb_shape.ar] = swb_depth_compress(L,D,shape,opt,V00);  % [cm]

[f0,sigma_bs,Vz] = swb_LF_model_old(xi,D,freq,'chu',swb_shape,param);  % old model
[f0,fbs,sigma_bs,Vz] = swb_LF_model(xi,D,freq,'chu',swb_shape,param);  % new model

% Compare old and newly corrected model
fig = figure;
plot(freq,10*log10(sigma_bs));
hold on
plot(freq,10*log10(sigma_bs_new),'-.');
xlabel('Frequency (Hz)');
ylabel('TS (dB)');
ll = legend('Old','New, corrected');
set(ll,'fontsize',12);
title(sprintf('LF swimbladder model, depth %d m',D));
grid on
% set(gca,'xscale','log')
ylim([-55 -30])

% Save figure
save_name = sprintf('%s_old_new_model_cmp.png',script_name);
saveSameSize(fig,'file',fullfile(save_path,save_name),...
    'format','png','renderer','painters');


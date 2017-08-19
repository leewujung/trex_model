% 2016 05 06  Fish LF scattering across depth

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
    addpath('F:\Dropbox\0_CODE\MATLAB\brewermap');
    base_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
else
    if ~fish_on_path
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\fish_scattering\fish_scat_model_new']);
    end
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    addpath('F:\Dropbox\0_CODE\MATLAB\brewermap');
    base_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests'];
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

param.gamma_alpha = 1.4;  % Nero et al. eq(4)
param.rho = 1050;         % from Chu's code
param.cw = 1525;
param.c_a = 343;
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
shape = 'swb';
V00 = 5;
opt = 2;
xi = 10;
freq = 1000:10:10000;
D_all = 1:2:19;

fig_ts = figure;
% colorset = parula(length(D_all));
colorset = brewermap(length(D_all),'Paired');

f0_all = zeros(size(D_all));
Vz_all = zeros(size(D_all));
sigma_bs_all = zeros(length(freq),length(D_all));
for iD=1:length(D_all)
    D = D_all(iD);
    [swb_shape.xg,swb_shape.ar] = swb_depth_compress(L,D,shape,opt,V00);  % [cm]
    [f0,fbs,sigma_bs,Vz] = swb_LF_model(xi,D,freq,'chu',swb_shape,param);  % new model
    
    f0_all(iD) = f0;
    Vz_all(iD) = Vz;
    sigma_bs_all(:,iD) = sigma_bs;
    
    figure(fig_ts)
    plot(freq,10*log10(sigma_bs),'color',colorset(length(D_all)-iD+1,:));
    hold on
end

figure(fig_ts)
ll = legend(num2str([D_all]'));
set(ll,'fontsize',10);
xlabel('Frequency (Hz)');
ylabel('TS (dB)');
ylim([-60 -35])
grid on
title('Swimbladder response vs depth')

fig_res = figure;
plot(D_all,f0_all,'-o');
xlabel('Depth (m)')
ylabel('Frequency (Hz)')
grid on
title('Resonant frequency')

% Save figure
save_ts = sprintf('%s_TS_xi%d.png',script_name,xi);
saveSameSize_150(fig_ts,'file',fullfile(save_path,save_ts),...
    'format','png','renderer','painters');
save_res = sprintf('%s_res_freq_xi%d.png',script_name,xi);
saveSameSize_150(fig_res,'file',fullfile(save_path,save_res),...
    'format','png','renderer','painters');


% 2016 04 27  Fish LF scattering across depth

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
    base_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
else
    if ~fish_on_path
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\fish_scattering\fish_scat_model_new']);
    end
    base_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests'];
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
D_all = 1:19;
opt = 1;
xi = 5;
freq = 0:4:5000;

for iD=1:length(D_all)
    D = D_all(iD);
[swb_shape.xg,swb_shape.ar] = swb_depth_compress(L,D,shape,opt,V00);  % [cm]

% figure;

% [f0,sigma_bs,Vz] = swb_LF_model(xi,D,freq,'nero',swb_shape,param);
% plot(freq,10*log10(sigma_bs));
% hold on

% [f0,sigma_bs,Vz] = swb_LF_model(xi,D,freq,'ye',swb_shape,param);
% plot(freq,10*log10(sigma_bs));

[f0,sigma_bs,Vz] = swb_LF_model(xi,D,freq,'chu',swb_shape,param);
% plot(freq,10*log10(sigma_bs));

% legend('Nero','Ye','Chu')
% title(sprintf('Depth %d m',D));

    f0_all(iD) = f0;
end

figure
plot(D_all,f0_all,'o-')
grid
xlabel('Depth (m)')
ylabel('Resonant frequency (Hz)')


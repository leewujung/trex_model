% 2016 04 25  Calculate time series at different ranges using the time
%             shited trick in FFT and compare with brute force solution

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
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

env_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_modes_calc';
env_file_pre = 'mfenv';


% Modeling parameters
tau = 1e-3;  % [sec]
fc = 2000;


%% Original time
freq_step = 2;
freq_all = freq_step:freq_step:4000;  % [Hz]
param_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_modes_calc';
param_file = 'env_param_trex13.mat';

env_param = load(fullfile(param_path,param_file));
z0 = 17.8;
r = 500;  % distance [m]
% delta_z = 1;
% z = 0:delta_z:30;
z = 5;

pp_freq = zeros(length(freq_all),length(z),length(r));
for iF=1:length(freq_all)
    freq = freq_all(iF);
    disp(['freq = ',num2str(freq),' Hz']);
    
    % load modes
    env_file = sprintf('%s_%04d',env_file_pre,freq);
    mod_file = [env_file,'.mod'];
    modes = read_modes(fullfile(env_path,mod_file));
    mode_incl = length(modes.k);
    
    if mode_incl~=0
        % get field from modes
        pp_freq(iF,:,:) = get_pfield(modes,mode_incl,r,z0,z);  % direct arrival
    end
end


% Time series at a particular location
pp_rz = pp_freq;
pp_rz_cmplx = [zeros(1,size(pp_rz,2));pp_rz;flipud(conj(pp_rz))];
tt_rz = ifft(pp_rz_cmplx);
tt_vec = linspace(0,1/freq_step,length(tt_rz));

fig_ir = figure;
corder = get(gca,'colororder');  % get line colors
plot(tt_vec,tt_rz,'color',corder(1,:));
hold on
axis ij
title('Impulse response')

% Transmit signal
S = exp(-pi^2*(freq_all(:)-fc).^2*tau^2);
S_cmplx = [0;S;flipud(conj(S))];
% s = ifftshift(ifft(S_cmplx));

tt_rz_mult = zeros(size(tt_rz,1),length(z));
for iZ=1:length(z)
    tt_rz_mult(:,iZ) = ifft(pp_rz_cmplx(:,iZ).*S_cmplx);
%     tt_rz_mult(:,iZ) = ifft(pp_rz_cmplx(:,iZ).*S');
end
tt_vec_mult = linspace(0,1/freq_step,length(tt_rz_mult));

fig_sig = figure;
plot(tt_vec_mult,tt_rz_mult,'color',corder(1,:));
hold on
axis ij
title('Time series with transmit signal')


%% Reduced time
freq_step = 4;
freq_all = freq_step:freq_step:4000;  % [Hz]
param_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_modes_calc';
param_file = 'env_param_trex13.mat';

env_param = load(fullfile(param_path,param_file));
z0 = 17.8;
r = 500;  % distance [m]
z = 5;

pp_freq = zeros(length(freq_all),length(z),length(r));
for iF=1:length(freq_all)
    freq = freq_all(iF);
    disp(['freq = ',num2str(freq),' Hz']);
    
    % load modes
    env_file = sprintf('%s_%04d',env_file_pre,freq);
    mod_file = [env_file,'.mod'];
    modes = read_modes(fullfile(env_path,mod_file));
    mode_incl = length(modes.k);
    
    if mode_incl~=0
        % get field from modes
        pp_freq(iF,:,:) = get_pfield(modes,mode_incl,r,z0,z);  % direct arrival
    end
end


% Time series at a particular location
pp_rz = pp_freq;
pp_rz_cmplx = [zeros(1,size(pp_rz,2));pp_rz;flipud(conj(pp_rz))];
tt_rz = ifft(pp_rz_cmplx);
tt_vec = linspace(0,1/freq_step,length(tt_rz));

figure(fig_ir)
plot(tt_vec,tt_rz,'color',corder(2,:));

% Apply reduced time shift
t_min = 0.15;   % start of reduced time [sec]
freq_all_cmplx_reduce = [0,freq_all,-fliplr(freq_all)];
phase_shift = exp(1j*2*pi*freq_all_cmplx_reduce(:)*t_min);
pp_rz_cmplx_reduce = pp_rz_cmplx.*repmat(phase_shift,1,size(pp_rz_cmplx,2));
tt_rz_reduce = ifft(pp_rz_cmplx_reduce);
tt_vec_reduce = linspace(0,1/freq_step,length(tt_rz_reduce));

figure(fig_ir)
plot(tt_vec_reduce+t_min,tt_rz_reduce,'color',corder(3,:));
xlabel('Time (sec)');
ll = legend('Brute force','Undersampled','Reduced time','location','northwest');
set(ll,'fontsize',12);

% Transmit signal
S = exp(-pi^2*(freq_all(:)-fc).^2*tau^2);
S_cmplx = [0;S;flipud(conj(S))];

tt_rz_mult_reduce = zeros(size(tt_rz,1),length(z));
for iZ=1:length(z)
    tt_rz_mult_reduce(:,iZ) = ifft(pp_rz_cmplx_reduce(:,iZ).*S_cmplx);
%     tt_rz_mult(:,iZ) = ifft(pp_rz_cmplx(:,iZ).*S');
end
tt_vec_mult_reduce = linspace(0,1/freq_step,length(tt_rz_mult_reduce));

figure(fig_sig);
plot(tt_vec_mult_reduce+t_min,tt_rz_mult_reduce,'color',corder(2,:));
axis ij
xlabel('Time (sec)');
title('Time series with transmit signal')
ll = legend('Brute force','Reduced time','location','northwest');
set(ll,'fontsize',12);




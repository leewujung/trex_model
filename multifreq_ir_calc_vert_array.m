% 2016 03 31  Multi-freq kraken runs for time domain prediction


clear 

usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    addpath(genpath('F:\Dropbox\0_APL_normal_mode\kraken'))
    base_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
else
    addpath(genpath(['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken']));
    base_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests'];
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

env_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_modes_calc';
env_file_pre = 'mfenv';

addpath(genpath('F:\Dropbox\0_APL_normal_mode\kraken'))

clear 

env_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_time_ir';
env_file_pre = 'mfenv';

freq_step = 5;
freq_all = 1:freq_step:4000;  % [Hz]
param_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_time_ir';
param_file = 'env_param_trex13.mat';

env_param = load(fullfile(param_path,param_file));
z0 = 17.8;
r = 20;  % distance [m]
delta_z = 1;
z = 0:delta_z:30;

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
        pp_freq(iF,:,:) = get_pfield(modes,mode_incl,r,z0,z);
    end
end

plot_shift = 0.01;

% Time series at a particular location
pp_rz = pp_freq;
pp_rz_cmplx = [zeros(1,size(pp_rz,2));pp_rz;flipud(conj(pp_rz))];
tt_rz = ifft(pp_rz_cmplx);
tt_vec = linspace(0,1/freq_step,length(tt_rz));

figure
plot(tt_vec,tt_rz+repmat(z,size(tt_rz,1),1)*plot_shift);
set(gca,'ytick',z*plot_shift,'yticklabel',z)
axis ij
title('Impulse response')

% Convolution in time
tt_delta = diff(tt_vec(1:2));
tau = 1e-3;  % [sec]
s_len = (size(pp_rz_cmplx,1)-1)/2;
s_t = (-s_len:s_len)*tt_delta;
fc = 1000;
s = exp(-s_t.^2/tau^2).*(exp(1i*2*pi*fc*s_t)+exp(-1i*2*pi*fc*s_t))*1/2;

tt_rz_conv = zeros(size(tt_rz,1)+length(s)-1,length(z));
for iZ=1:length(z)
    tt_rz_conv(:,iZ) = conv(tt_rz(:,iZ),s);
end
tt_vec_conv = (0:size(tt_rz_conv,1)-1)*tt_delta;

figure
plot(tt_vec_conv+s_t(1),tt_rz_conv+repmat(z,size(tt_rz_conv,1),1)*plot_shift);
set(gca,'ytick',z*plot_shift,'yticklabel',z)
axis ij
title('Convolution in time')

% Multiplication in freq
S = fft(s);
s_freq = linspace(0,1/tt_delta,length(S));

tt_rz_mult = zeros(size(tt_rz,1),length(z));
for iZ=1:length(z)
    tt_rz_mult(:,iZ) = ifftshift(ifft(pp_rz_cmplx(:,iZ).*S.'));
end
tt_vec_mult = linspace(0,1/freq_step,length(tt_rz_mult));

figure
plot(tt_vec_mult,tt_rz_mult+repmat(z,size(tt_rz_mult,1),1)*plot_shift);
set(gca,'ytick',z*plot_shift,'yticklabel',z)
axis ij
title('Mulitplication in frequency')




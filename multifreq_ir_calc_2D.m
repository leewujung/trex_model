% 2016 03 31  Multi-freq kraken runs for time domain prediction

% addpath(genpath('F:\Dropbox\0_APL_normal_mode\kraken'))

clear 

env_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_time_ir';
env_file_pre = 'mfenv';

freq_step = 1;
freq_all = 1:freq_step:2000;  % [Hz]
param_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_time_ir';
param_file = 'env_param_trex13.mat';

env_param = load(fullfile(param_path,param_file));
z0 = 17.8;
r = 0:2:200;  % distance [m]
delta_z = 0.5;
z = 0:delta_z:30;

pp_freq = zeros(length(freq_all),length(z),length(r));
for iF=1:length(freq_all)
    freq = freq_all(iF);
    disp(['freq = ',num2str(freq),' Hz']);
    
    % load modes
    env_file = sprintf('%s_%04d',env_file_pre,freq);
%     env_file = sprintf('%s_%d',env_file_pre,freq);
    mod_file = [env_file,'.mod'];
    modes = read_modes(fullfile(env_path,mod_file));
    mode_incl = length(modes.k);
    
    if mode_incl~=0
        % get field from modes
        pp_freq(iF,:,:) = get_pfield(modes,mode_incl,r,z0,z);
    end
end

% Time series at a particular location
% r_idx = 251;
z_idx = 21;
pp_rz = squeeze(pp_freq(:,z_idx,:));
pp_rz_cmplx = [zeros(1,size(pp_rz,2));pp_rz;flipud(conj(pp_rz))];
% pp_rz = squeeze(pp_freq(:,z_idx,r_idx));
% pp_rz_cmplx = [0;pp_rz;flipud(conj(pp_rz))];
tt_rz = ifft(pp_rz_cmplx);
tt_vec = linspace(0,1/10,length(tt_rz));


% Time series across the plane
tt_rz = zeros(size(pp_freq,1)*2+1,length(z),length(r));
for iZ=1:length(z)
    for iR=1:length(r)
        r_idx = iR;
        z_idx = iZ;
        pp_rz = squeeze(pp_freq(:,z_idx,r_idx));
        pp_rz_cmplx = [0;pp_rz;flipud(conj(pp_rz))];
        tt_rz(:,iZ,iR) = ifft(pp_rz_cmplx);
    end
end

tt_rz_real = tt_rz;
idx = imag(tt_rz_real)~=0;
tt_rz_real(idx) = real(tt_rz_real(idx));

% tt_delta = 1/(2*max(freq_all));
tt_vec = linspace(0,1/2,size(tt_rz,1));

% Convolution in time
tt_delta = diff(tt_vec(1:2));
tau = 1e-3;  % [sec]
s_t = -250e-3:tt_delta:250e-3;
fc = 500;
s = exp(-s_t.^2/tau^2).*(exp(1i*2*pi*fc*s_t)+exp(-1i*2*pi*fc*s_t))*1/2;

tt_rz_conv = zeros(size(tt_rz,1)+length(s)-1,length(z),length(r));
for iZ=1:length(z)
    for iR=1:length(r)
        tt_rz_conv(:,iZ,iR) = conv(squeeze(tt_rz(:,iZ,iR)),s);
    end
end

tt_rz_conv_real = tt_rz_conv;
idx = imag(tt_rz_conv_real)~=0;
tt_rz_conv_real(idx) = real(tt_rz_conv_real(idx));


% Plot
fig = figure;
for iT=2000:1:size(tt_rz,1)
    figure(fig)
    p = squeeze(tt_rz_conv_real(iT,:,:));
    imagesc(r,z,p);
    colormap(jet)
    hold on
    plot([r(1) r(end)],[19,19],'w--');
    caxis([-.03 .03])
    colorbar
    title(sprintf('%.2f ms',tt_vec(iT)*1e3));
    hold off
    drawnow
    xlabel('Range (m)');
    ylabel('Depth (m)');
end



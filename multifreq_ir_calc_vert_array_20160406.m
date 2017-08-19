% 2016 03 31  Multi-freq kraken runs for time domain prediction

% addpath(genpath('F:\Dropbox\0_APL_normal_mode\kraken'))

clear 

env_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_modes_calc';
env_file_pre = 'mfenv';

freq_step = 1;
freq_all = 1:freq_step:4000;  % [Hz]
param_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_modes_calc';
param_file = 'env_param_trex13.mat';

env_param = load(fullfile(param_path,param_file));
z0 = 17.8;
r = 2000;  % distance [m]
delta_z = 0.1;
z = 0:delta_z:22;

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
        pp_freq(iF,:,:) = get_pfield(modes,mode_incl,r,z0,z).^2;
    end
end


% Impulse response at a particular location
pp_freq_cmplx = [zeros(1,size(pp_freq,2));pp_freq;flipud(conj(pp_freq))];
tt_rz = ifft(pp_freq_cmplx);
tt_vec = linspace(0,1/freq_step,length(tt_rz));

% Transmit signal
tt_delta = diff(tt_vec(1:2));
tau = 1e-3;  % [sec]
s_len = (size(pp_freq_cmplx,1)-1)/2;
s_t = (-s_len:s_len)*tt_delta;
fc = 3000;
s = exp(-s_t.^2/tau^2).*(exp(1i*2*pi*fc*s_t)+exp(-1i*2*pi*fc*s_t))*1/2;
S = fft(s);
s_freq = linspace(0,1/tt_delta,length(S));
S_mtx = repmat(S.',1,length(z));

% % Convolution in time
% tt_rz_conv = zeros(size(tt_rz,1)+length(s)-1,length(z));
% for iZ=1:length(z)
%     tt_rz_conv(:,iZ) = conv(tt_rz(:,iZ),s);
% end
% tt_vec_conv = (0:size(tt_rz_conv,1)-1)*tt_delta;
% 
% figure
% plot(tt_vec_conv+s_t(1),tt_rz_conv+repmat(z,size(tt_rz_conv,1),1)*plot_shift);
% set(gca,'ytick',z*plot_shift,'yticklabel',z)
% axis ij
% title('Convolution in time')

% Multiplication in freq
% tt_rz_mult = fftshift(ifft(pp_freq_cmplx.*S_mtx));
tt_rz_mult = zeros(size(tt_rz,1),length(z));
for iZ=1:length(z)
    tt_rz_mult(:,iZ) = ifftshift(ifft(pp_freq_cmplx(:,iZ).*S.'));
end
tt_vec_mult = linspace(0,1/freq_step,length(tt_rz_mult));

plot_shift = max(max(tt_rz_mult,[],1))*0.8;


figure
plot(tt_vec,tt_rz+repmat(z,size(tt_rz,1),1)*plot_shift);
set(gca,'ytick',z*plot_shift,'yticklabel',z)
axis ij
title('Impulse response')

figure
subplot(211)
plot(s_t,s);
xlabel('Time (sec)')
title(['Center freq = ',num2str(fc),'Hz'])
subplot(212)
plot(s_freq,abs(S));
xlabel('Frequency (Hz)')

figure
plot(tt_vec_mult,tt_rz_mult+repmat(z,size(tt_rz_mult,1),1)*plot_shift);
set(gca,'ytick',(0:5:30)*plot_shift,'yticklabel',num2str([0:5:30]'))
axis ij
% title('Mulitplication in frequency')
title(['Center freq = ',num2str(fc),'Hz, Dist = ',num2str(r),'m'])
xlabel('Time (sec)')
ylabel('Depth (m)')

figure
imagesc(tt_vec_mult,z,10*log10(tt_rz_mult'.^2));
hold on
plot(tt_vec_mult([1 end]),[19 19],'w--','linewidth',2)
caxis([-150 -80])
xlabel('Time (sec)')
ylabel('Depth (m)')

% Integrating energy
e = trapz(tt_vec_mult,tt_rz_mult.^2);
figure
plot(e,z);
axis ij
title(['Center freq = ',num2str(fc),'Hz, Dist = ',num2str(r),'m'])



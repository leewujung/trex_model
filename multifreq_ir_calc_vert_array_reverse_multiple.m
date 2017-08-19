% 2016 03 31  Multi-freq kraken runs for time domain prediction

% addpath(genpath('F:\Dropbox\0_APL_normal_mode\kraken'))

clear

usrn = getenv('username');

save_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160406_ensemble';
save_file = 'ensemble_test.mat';
if ~exist(save_path,'dir')
    mkdir(save_path);
end

env_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_modes_calc';
env_file_pre = 'mfenv';

freq_step = 2;
freq_all = 1:freq_step:2000;  % [Hz]
param_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_modes_calc';
param_file = 'env_param_trex13.mat';

env_param = load(fullfile(param_path,param_file));
z0 = 17.8;
% r = 1000;  % distance [m]
% delta_z = 1;
% z = 0:delta_z:30;

% Fish ball param
r_mean = 500;
r_sigma = 20;
% r_sigma = [10:10:50];
z_mean = [3:3:18];
% z_mean = 10;
z_sigma = 2;
fish_num = 100;
rep_num = 100;

pp_freq_fish_rep = zeros(length(freq_all),fish_num,length(z_mean),rep_num);
for iRP=1:rep_num
    disp(['Rep ',num2str(iRP)]);
    for iZ=1:length(z_mean)
        disp(['Depth = ',num2str(z_mean(iZ)),' m']);
        rr = randn(1,fish_num)*r_sigma+r_mean;
        zz = randn(1,fish_num)*z_sigma+z_mean(iZ);
        parfor iF=1:length(freq_all)
            freq = freq_all(iF);
%             disp(['freq = ',num2str(freq),' Hz']);
            
            % load modes
            env_file = sprintf('%s_%04d',env_file_pre,freq);
            mod_file = [env_file,'.mod'];
            modes = read_modes(fullfile(env_path,mod_file));
            mode_incl = length(modes.k);
            
            if mode_incl~=0
                for iFISH=1:fish_num
                    pp_freq_fish_rep(iF,iFISH,iZ,iRP) = get_pfield(modes,mode_incl,rr(iFISH),z0,zz(iFISH)).^2;  % roundtrip
                end
            end
        end
    end
end
save(fullfile(save_path,save_file));
pp_freq = squeeze(nansum(pp_freq_fish,2));


% Impulse response at a particular location
pp_freq_cmplx = [zeros(1,size(pp_freq,2));pp_freq;flipud(conj(pp_freq))];
tt_rz = ifft(pp_freq_cmplx);
tt_vec = linspace(0,1/freq_step,length(tt_rz));

% figure
% plot(tt_vec,tt_rz+repmat(z,size(tt_rz,1),1)*plot_shift);
% set(gca,'ytick',z*plot_shift,'yticklabel',z)
% axis ij
% title('Impulse response')

% Transmit signal
tt_delta = diff(tt_vec(1:2));
tau = 1e-3;  % [sec]
s_len = (size(pp_freq_cmplx,1)-1)/2;
s_t = (-s_len:s_len)*tt_delta;
fc = 1000;
s = exp(-s_t.^2/tau^2).*(exp(1i*2*pi*fc*s_t)+exp(-1i*2*pi*fc*s_t))*1/2;
S = fft(s);
s_freq = linspace(0,1/tt_delta,length(S));

% figure
% subplot(211)
% plot(s_t,s);
% xlabel('Time (sec)')
% title(['Center freq = ',num2str(fc),'Hz'])
% subplot(212)
% plot(s_freq,abs(S));
% xlabel('Frequency (Hz)')

% Multiplication in freq
% tt_rz_mult = zeros(size(tt_rz,1),length(z_mean));
% tt_rz_mult_env = zeros(size(tt_rz,1),length(z_mean));
% for iZ=1:length(z_mean)
%     tt_rz_mult(:,iZ) = ifftshift(ifft(pp_freq_cmplx(:,iZ).*S.'));
%     tt_rz_mult_env(:,iZ) = abs(hilbert(tt_rz_mult(:,iZ)));
% end
tt_rz_mult = zeros(size(tt_rz,1),length(r_sigma));
tt_rz_mult_env = zeros(size(tt_rz,1),length(r_sigma));
for iR=1:length(r_sigma)
    tt_rz_mult(:,iR) = ifftshift(ifft(pp_freq_cmplx(:,iR).*S.'));
    tt_rz_mult_env(:,iR) = abs(hilbert(tt_rz_mult(:,iR)));
end
tt_vec_mult = linspace(0,1/freq_step,length(tt_rz_mult));

plot_shift = max(max(tt_rz_mult,[],1))*0.8;


%%%%% z_mean variation
figure
plot(tt_vec_mult,tt_rz_mult+repmat(z_mean,size(tt_rz_mult,1),1)*plot_shift);
set(gca,'ytick',z_mean*plot_shift,'yticklabel',z_mean)
axis ij
% title('Mulitplication in frequency')
title(['Center freq = ',num2str(fc),'Hz, Dist = ',num2str(rr_mean),'m'])
xlabel('Time (sec)')
ylabel('Depth (m)')

figure
plot(tt_vec_mult,-tt_rz_mult_env+repmat(z_mean,size(tt_rz_mult,1),1)*plot_shift);
set(gca,'ytick',z_mean*plot_shift,'yticklabel',z_mean)
axis ij
% title('Mulitplication in frequency')
title(['Center freq = ',num2str(fc),'Hz, Dist = ',num2str(rr_mean),'m'])
xlabel('Time (sec)')
ylabel('Depth (m)')

% Integrating energy
e = trapz(tt_vec_mult,tt_rz_mult.^2);
figure
plot(e,z_mean);
axis ij
title(['Center freq = ',num2str(fc),'Hz, Dist = ',num2str(rr_mean),'m'])


%%%%% r_sigma variation
figure
plot(tt_vec_mult,tt_rz_mult+repmat(r_sigma/5,size(tt_rz_mult,1),1)*plot_shift);
set(gca,'ytick',z_mean*plot_shift,'yticklabel',z_mean)
axis ij
% title('Mulitplication in frequency')
title(['Center freq = ',num2str(fc),'Hz, Dist = ',num2str(rr_mean),'m'])
xlabel('Time (sec)')
ylabel('Depth (m)')

figure
plot(tt_vec_mult,-tt_rz_mult_env+repmat(r_sigma/5,size(tt_rz_mult,1),1)*plot_shift);
set(gca,'ytick',z_mean*plot_shift,'yticklabel',z_mean)
axis ij
% title('Mulitplication in frequency')
title(['Center freq = ',num2str(fc),'Hz, Dist = ',num2str(rr_mean),'m'])
xlabel('Time (sec)')
ylabel('Depth (m)')

figure
imagesc(tt_vec_mult,r_sigma/5,tt_rz_mult_env')
xlabel('Time (sec)')
ylabel('\sigma_{range} (m)')

% Integrating energy
e = trapz(tt_vec_mult,tt_rz_mult.^2);
figure
plot(e,r_sigma/5);
axis ij
title(['Center freq = ',num2str(fc),'Hz, Dist = ',num2str(rr_mean),'m'])


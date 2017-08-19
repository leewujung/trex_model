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
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    base_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
    mode_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\\multifreq_mode_calc_20160427';
    param_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\\multifreq_mode_calc_20160427';
else
    if ~kraken_on_path   % add path if Kraken is not already on path
        addpath(genpath(['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken']));
    end
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    base_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests'];
    mode_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\\multifreq_mode_calc_20160427'];
    param_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\\multifreq_mode_calc_20160427'];
end
env_file_pre = 'mfenv';

% Modeling parameters
freq_step = 4;
freq_all = freq_step:freq_step:8000;  % [Hz]
param_file = 'env_param_trex13.mat';
env_param = load(fullfile(param_path,param_file));

rt = 0;     % direct reception (0) or round trip (1)
z0 = 5;  % source depth [m]
z = 0:0.1:21;   % receiver depth [m]

% Transmit signal parameters
tau = 1e-3;  % [sec]
fc = 3000;   % [Hz]

% Store parameters
param.env_param_path = param_path;
param.env_param_file = param_file;
param.env_param = env_param;
param.mode_path = mode_path;
if rt
    param.scenario = 'round trip';
else
    param.scenario = 'direct reception';
end
param.z0 = z0;
param.z = z;
param.tx_tau = tau;
param.tx_fc = fc;
param.freq_step = freq_step;
param.freq_max = freq_all(end);

% Save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Loop through different distance
% r_all = [100:100:1000,1200:200:3000,3500:500:10000];
r_all = 11000:1000:20000;

fig_2d_lin = figure('units','normalized','outerposition',[0 0 1 1]);
fig_2d_log = figure('units','normalized','outerposition',[0 0 1 1]);
% v = VideoWriter(fullfile(save_path,'video.mp4'),'mpeg-4');
% v.FrameRate = 1;
% v.Quality = 90;
% open(v);

for iR=1:length(r_all)
    
r = r_all(iR);   % distance [m]
if rt
    t_1st_arr = r/env_param.SSP{1}.cp(1)*2;  % *2 for roundtrip
else
    t_1st_arr = r/env_param.SSP{1}.cp(1);
end
T = 1/freq_step;
% t_min = round((t_1st_arr-T/8)*100)/100;  % start of reduced time [sec]
t_min = t_1st_arr-0.01;  % start of reduced time [sec]
if t_min<0
    t_min = 0;
end

param.r = r;    % range of receiver/scatterer [m]
param.T = T;    % time window length [sec]
param.t_min = t_min; % reduced time start [sec]

% Assemble across frequency
pp_freq = zeros(length(freq_all),length(z),length(r));
for iF=1:length(freq_all)
    freq = freq_all(iF);
    disp(['freq = ',num2str(freq),' Hz']);
    
    % load modes
    env_file = sprintf('%s_%06.1f',env_file_pre,freq);
    mod_file = [env_file,'.mod'];
    modes = read_modes(fullfile(mode_path,mod_file));
    mode_incl = length(modes.k);
    
    if mode_incl~=0
        % get field from modes
        if rt
            pp_freq(iF,:,:) = get_pfield(modes,mode_incl,r,z0,z).^2;  % echoes
        else
            pp_freq(iF,:,:) = get_pfield(modes,mode_incl,r,z0,z);  % receiving only
        end
    end
end

% Time series at a particular location
pp_rz = pp_freq;
pp_rz_cmplx = [zeros(1,size(pp_rz,2));pp_rz;flipud(conj(pp_rz))];
tt_rz = ifft(pp_rz_cmplx);
tt_vec = linspace(0,1/freq_step,length(tt_rz));

% Apply reduced time shift
freq_all_cmplx_reduce = [0,freq_all,-fliplr(freq_all)];
phase_shift = exp(1j*2*pi*freq_all_cmplx_reduce(:)*t_min);
pp_rz_cmplx_reduce = pp_rz_cmplx.*repmat(phase_shift,1,size(pp_rz_cmplx,2));
% tt_rz_reduce = ifft(pp_rz_cmplx_reduce);
% tt_vec_reduce = linspace(0,1/freq_step,length(tt_rz_reduce));

% Transmit signal
S = exp(-pi^2*(freq_all(:)-fc).^2*tau^2);
S_cmplx = [0;S;flipud(conj(S))];

% Calculate time series
tt_rz_mult_reduce = ifft(pp_rz_cmplx_reduce.*repmat(S_cmplx,1,size(pp_rz_cmplx_reduce,2)));
tt_vec_mult_reduce = linspace(0,1/freq_step,length(tt_rz_mult_reduce));

% Get envelope
tt_rz_mult_reduce_env = zeros(size(tt_rz_mult_reduce));
for iZ=1:length(z)
    tt_rz_mult_reduce_env(:,iZ) = abs(hilbert(tt_rz_mult_reduce(:,iZ)));
end

% % Raw time series
% figure;
% corder = get(gca,'colororder');  % get line colors
% plot_shift = max(tt_rz_mult_reduce_env(:))*0.8;
% [~,idx] = sort(-z*plot_shift);
% z_new = (z-min(z));
% z_min = sort(z_new);
% z_min = z_min(2);
% z_new = -z_new/z_min*plot_shift;
% z_bottom = 19;
% z_bottom_new = -(z_bottom-min(z))/z_min*plot_shift;
% 
% plot(tt_vec_mult_reduce+t_min,...
%     tt_rz_mult_reduce+repmat(z_new,size(tt_rz_mult_reduce,1),1),...
%     'color',corder(1,:));
% % hold on
% % plot(tt_vec_mult_reduce([1 end])+t_min,z_bottom_new*[1 1],'k');
% set(gca,'ytick',sort(z_new),'yticklabel',z(idx));  % label yaxis correctly
% xlabel('Time (sec)');
% ylabel('Depth (z)');
% title('Time series with transmit signal')
% 
% % Envelope
% figure;
% plot(tt_vec_mult_reduce+t_min,...
%     tt_rz_mult_reduce_env+repmat(z_new,size(tt_rz_mult_reduce,1),1),...
%     'color',corder(1,:));
% % hold on
% % plot(tt_vec_mult_reduce([1 end])+t_min,z_bottom_new*[1 1],'k');
% set(gca,'ytick',sort(z_new),'yticklabel',z(idx))  % label yaxis correctly
% xlabel('Time (sec)');
% ylabel('Depth (z)');
% title('Time series with transmit signal')

if rt
    title_txt = sprintf('Range=%dm, fc=%dHz, z0=%3.1fm, round trip',r,fc,z0);
else
    title_txt = sprintf('Range=%dm, fc=%dHz, z0=%3.1fm, one-way',r,fc,z0);
end

plot_time_span = 0.12;

% 2D color plot: linear scale
figure(fig_2d_lin);
cla
imagesc(tt_vec_mult_reduce+t_min,z,tt_rz_mult_reduce_env');
hold on
plot(tt_vec_mult_reduce([1 end])+t_min,[19 19],'w--','linewidth',1)
xlim(tt_vec_mult_reduce(1)+t_min+[0 plot_time_span])
xlabel('Time (sec)');
ylabel('Depth (m)');
title(title_txt);

% 2D color plot: log scale
mm = 20*log10(tt_rz_mult_reduce_env(:)');
mm_med = median(mm);
mm_std = std(mm);
mm_max = 10*round(max(mm)/10);
mm_min = mm_max-50;                 % min for colorbar

figure(fig_2d_log);
cla
imagesc(tt_vec_mult_reduce+t_min,z,20*log10(tt_rz_mult_reduce_env'));
hold on
plot(tt_vec_mult_reduce([1 end])+t_min,[19 19],'w--','linewidth',1)
xlim(tt_vec_mult_reduce(1)+t_min+[0 plot_time_span])
xlabel('Time (sec)');
ylabel('Depth (m)');
caxis([mm_min mm_max])
colorbar
title(title_txt);

% Save figure
save_name_lin = sprintf('%s_lin_r%05dm.png',script_name,r);
save_name_log = sprintf('%s_log_r%05dm.png',script_name,r);
saveSameSize(fig_2d_lin,'file',fullfile(save_path,save_name_lin),...
    'format','png','renderer','painters');
saveSameSize(fig_2d_log,'file',fullfile(save_path,save_name_log),...
    'format','png','renderer','painters');

% Save parameters
save_name_param = sprintf('%s_param_r%05dm.mat',script_name,r);
save(fullfile(save_path,save_name_param),'param');

% % Grab frame into video
%   writeVideo(v,getframe(fig_2d_lin));
  
end  % loop through all ranges

% close(v)
close(fig_2d_lin)
close(fig_2d_log)
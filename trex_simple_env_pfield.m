% 2017 08 17  Assemble pressure field/TL from all modes
%             Compare across TREX13 frequencies

addpath 'F:\Dropbox\0_APL_normal_mode\kraken\Matlab\wjlee_codes'
addpath 'F:\Dropbox\0_APL_normal_mode\kraken\Matlab\ReadWrite'



clear

% % Save path
% [~,script_name,~] = fileparts(mfilename('fullpath'));
% save_path = fullfile(base_path,script_name);
% if ~exist(save_path,'dir')
%     mkdir(save_path);
% end

% Set up path and params
freq_all = 1800:200:3600;  % [Hz]
r = 0:1:5000;  % distance [m]
z = 0:0.1:21;   % receiver depth [m]
z0 = 17.8;       % source depth [m]
p_all = zeros(length(freq_all),length(z),length(r));

% Individual frequency pressure field
for iF=1:length(freq_all)
    mode_path = pwd;
    mode_file_pre = 'trex13_env';
    mode_file = sprintf('%s_%06.1fHz.mod',mode_file_pre,freq_all(iF));
    modes = read_modes(fullfile(mode_path,mode_file));
    phi = get_component(modes,'N');
    p = get_pfield(modes,length(modes.k),r,z0,z);
    [tlmin,tlmax] = tl_plot_min_max(p);
    tl = -20*log10(abs(p));

    p_all(iF,:,:) = p;
    
    figure;
    imagesc(r,z,tl);
    hold on
    plot(r([1 end]),[19 19],'w--','linewidth',2)
    xlabel('Range (m)')
    ylabel('Depth (m)')
    colormap(flipud(parula))
    caxis([tlmin tlmax])
    title(sprintf('freq=%dHz, z0=%3.1fm',freq_all(iF),z0));
end

% Mean pressure field
p_mean = squeeze(mean(p_all,1));
[tl_mean_min,tl_mean_max] = tl_plot_min_max(p_mean);
tl_mean = -20*log10(abs(p_mean));
    
figure;
imagesc(r,z,tl_mean);
hold on
plot(r([1 end]),[19 19],'w--','linewidth',2)
xlabel('Range (m)')
ylabel('Depth (m)')
colormap(flipud(parula))
caxis([tl_mean_min tl_mean_max])
title(sprintf('Mean pressure field, z0=%3.1fm',freq_all(iF),z0));

% Calculate depth pressure field average
depth_idx = 1:20:211;
p_mean_depth = zeros(length(depth_idx),length(r));
tl_mean_depth = zeros(length(depth_idx),length(r));
for iD=1:length(depth_idx)-1
    p_mean_depth(iD,:) = smooth(mean(p_mean(depth_idx(iD):depth_idx(iD+1),:),1),100);
    tl_mean_depth(iD,:) = smooth(mean(p_mean(depth_idx(iD):depth_idx(iD+1),:).^2,1),100);
end
%tl_mean_depth = 20*log10(abs(p_mean_depth));




%     saveSameSize_100(fig_all_mode,'file',fullfile(save_path,save_all_mode),...
%         'format','png','renderer','painters');

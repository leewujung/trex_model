% 2016 03 25  Calculate pressure field

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

% Save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% mode_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
fc = 1000;  % cetner frequency for the transmit signal [Hz]

mode_file = sprintf('mfenv_%06.1f.mod',fc);
modes = read_modes(fullfile(mode_path,mode_file));
phi = get_component(modes,'N');

kr = modes.k;

r = 0:1:3000;  % distance [m]
z = 0:0.1:21;   % receiver depth [m]
z0 = 15;       % source depth [m]

fig_all_mode = figure('units','normalized','outerposition',[0 0 1 1]);
% fig_one_mode = figure('units','normalized','outerposition',[0 0 1 1]);
mode_incl_all = 1:length(modes.k);
for iM=1:30
    mode_incl = mode_incl_all(iM);
    p = get_pfield(modes,mode_incl,r,z0,z);
    p_one_mode = get_pfield_one_mode(modes,mode_incl,r,z0,z);

    p_abs = abs(p);
    icount = find(p_abs>1e-7);         % for stats, only these values count <-- taken from Kraken
    p_abs(p_abs<1e-7) = 1e-7;
    p_abs(isnan(p_abs)) = 1e-6;
    p_abs(isinf(p_abs)) = 1e-6;
    
    tl = -20*log10(p_abs);
    
    % Get color scale for full field
    tlmed = median(tl(icount));	% median value
    tlstd = std(tl(icount));     % standard deviation
    tlmax = tlmed+0.75*tlstd;       % max for colorbar
    tlmax = 10*round(tlmax/10);     % make sure the limits are round numbers
    tlmin = tlmax-50;               % min for colorbar
    if iM==1
        tlmin_const = tlmin;
        tlmax_const = tlmax;
    elseif exist('tlmax_const','var')
        if tlmax<tlmax_const
            tlmax_const = tlmax;
            tlmin_const = tlmin;
        end
    end
    
    % Get color scale for only one mode
    p_one_mode_abs = abs(p_one_mode);
    icount = find(p_one_mode_abs>1e-7);         % for stats, only these values count <-- taken from Kraken
    p_one_mode_abs(p_one_mode_abs<1e-7) = 1e-7;
    p_one_mode_abs(isnan(p_one_mode_abs)) = 1e-6;
    p_one_mode_abs(isinf(p_one_mode_abs)) = 1e-6;
    
    tl_one_mode = -20*log10(p_one_mode_abs);
    
    tlmed_one_mode = median(tl_one_mode(icount));	% median value
    tlstd_one_mode = std(tl_one_mode(icount));     % standard deviation
    tlmax_one_mode = tlmed_one_mode+0.75*tlstd_one_mode;       % max for colorbar
    tlmax_one_mode = 10*round(tlmax_one_mode/10);     % make sure the limits are round numbers
    tlmin_one_mode = tlmax_one_mode-50;               % min for colorbar
    
    % Plot
    figure(fig_all_mode)
    imagesc(r,z,tl);
    hold on
    plot(r([1 end]),[19 19],'w--','linewidth',2)
    xlabel('Range (m)')
    ylabel('Depth (m)')
    colormap(flipud(parula))
    caxisrev([tlmin tlmax])
    title(sprintf('Modes included: %d, fc=%dHz, z0=%3.1fm',mode_incl_all(iM),fc,z0));

%     figure(fig_one_mode)
%     imagesc(r,z,tl_one_mode);
%     hold on
%     plot(r([1 end]),[19 19],'w--','linewidth',2)
%     xlabel('Range (m)')
%     ylabel('Depth (m)')
%     colormap(flipud(parula))
%     caxisrev([tlmin_one_mode tlmax_one_mode])
%     title(sprintf('Mode=%d, fc=%d Hz, z0=%d m',mode_incl_all(iM),fc,z0));
    
    % Save figures
    save_all_mode = sprintf('%s_incl_all_mode_%02d.png',script_name,mode_incl_all(iM));
    saveSameSize_100(fig_all_mode,'file',fullfile(save_path,save_all_mode),...
        'format','png','renderer','painters');
%     save_one_mode = sprintf('%s_mode_%02d.png',script_name,mode_incl_all(iM));
%     saveSameSize(fig_one_mode,'file',fullfile(save_path,save_one_mode),...
%         'format','png','renderer','painters');

end
close(fig_all_mode)
% close(fig_one_mode)
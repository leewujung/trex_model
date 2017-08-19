% 2016 05 16  Compare the excitation amplitude at different source depth
%             for different modes. This was done to investigate the
%             "null"-like pattern in the depth-dependent time series across
%             different modes

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
        addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    end
    base_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
else
    if ~kraken_on_path   % add path if Kraken is not already on path
        addpath(genpath(['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken']));
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    end
    base_path = ['C:\Users\',usrn,'\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests'];
end


% Save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Set params
freq = 3000;  % frequency [Hz]
mode_path = 'multifreq_mode_calc_20160427';
mode_path = fullfile(base_path,mode_path);
mode_file = sprintf('mfenv_%06.1f.mod',freq);

% Read modes
modes = read_modes(fullfile(mode_path,mode_file));
phi = get_component(modes,'N');

title_txt = sprintf('freq=%06.1f Hz',freq);
save_fname_pre = sprintf('%s_f%06.1fHz',script_name,freq);

% Plot all modes
fig_k = figure;
plot(real(modes.k),imag(modes.k),'.');
hold on
text(double(real(modes.k)),double(imag(modes.k)),num2str([1:length(modes.k)]'),'fontsize',8);
ylim([-0.3 0])
xlabel('Real(k)')
ylabel('Imag(k)')
title(title_txt);
saveas(fig_k,fullfile(save_path,[save_fname_pre,'_k.fig']),'fig');
saveSameSize_150(fig_k,'file',fullfile(save_path,[save_fname_pre,'_k.png']),...
    'format','png','renderer','painters');

% Get excitation amplitude across modes
src_depth = 5;
[~,z_idx] = min(abs(src_depth-modes.z));
exc_amp = real(modes.phi(z_idx,:));
title_txt = sprintf('freq=%06.1f Hz, src depth=%03.1f m',freq,src_depth);
save_fname_pre = sprintf('%s_f%06.1fHz_src%03.1fm',script_name,freq,src_depth);

fig_exc = figure;
plot(1:length(exc_amp),exc_amp,'o-');
hold on
plot(1:length(exc_amp),abs(exc_amp),'o-');
xlabel('Mode number');
ylabel('Mode amplitude (real part)');
title(title_txt);
legend('mode amp','abs(mode amp)');
saveas(fig_exc,fullfile(save_path,[save_fname_pre,'_exc_amp.fig']),'fig');
saveSameSize_150(fig_exc,'file',fullfile(save_path,[save_fname_pre,'_exc_amp.png']),...
    'format','png','renderer','painters');


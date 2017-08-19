% 2016 05 16  Compare the mode location due to variation of NMESH

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
mode_path = 'multifreq_mode_calc_20160427';
mode_path = fullfile(base_path,mode_path);
mode_num_wanted = 25;
freq_step = 8;
freq_all = freq_step:freq_step:4000;  % [Hz]

% Transmit signal
tau = 1e-3;  % [sec]
fc = 3000;   % [Hz]
S = exp(-pi^2*(freq_all(:)-fc).^2*tau^2);
% S_cmplx = [0;S;flipud(conj(S))];

% Get all mode info
modes_k = zeros(length(freq_all),mode_num_wanted);
fig = figure;
for iF=1:length(freq_all)
    % freq = 3000;  % frequency [Hz]
    freq = freq_all(iF);
    mode_file = sprintf('mfenv_%06.1f.mod',freq);
    
    % Read modes
    modes = read_modes(fullfile(mode_path,mode_file));

    modes_uplim = min([mode_num_wanted length(modes.k)]);
    modes_k(iF,1:modes_uplim) = modes.k(1:modes_uplim);
    
    figure(fig)
    plot(real(modes.k),imag(modes.k),'.');
    hold on
end


S_modify = repmat(S,1,mode_num_wanted).*modes_k;
fig_mode = figure;
fig_mode_norm = figure;
colorset = jet(mode_num_wanted);
for iM=1:mode_num_wanted
    figure(fig_mode)
    plot(freq_all,abs(S_modify(:,iM)),'color',colorset(iM,:))
    hold on
    
    figure(fig_mode_norm)
    plot(freq_all,abs(S_modify(:,iM))./max(abs(S_modify(:,iM))),'color',colorset(iM,:))
    hold on
end




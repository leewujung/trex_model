% 2016 03 31  Multi-freq kraken runs for time domain prediction

clear 

env_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160331_time_ir';
if ~exist(env_path,'dir')
    mkdir(env_path);
end

freq_all = 1:1:4000;  % [Hz]
% freq_all = 1000:1000:4000;
root_finder = 1;
param_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160331_time_ir';
param_file = 'env_param_trex13.mat';
env_file_pre = 'multifreq_envfile';
env_title = 'Multi-freq TREX13 env';

ori_path = pwd;
cd(env_path);

runkraken = which( 'krakenc.exe' );

env_param = load(fullfile(param_path,param_file));
z0 = 10;
r = 0:10:500;  % distance [m]
delta_z = 0.2;
% z = 0:delta_z:env_param.SSP{end}.z(end);
z = 0:delta_z:25;

tic
pp_freq = zeros(length(freq_all),length(z),length(r));
parfor iF=1:length(freq_all)
    freq = freq_all(iF);
    disp(['freq = ',num2str(freq),' Hz']);
    
    % generate ENV file
%     rnd_str = [num2str(round(rand(1)*1e5))];
    env_file = [env_file_pre,num2str(freq)];
    gen_env_file(freq,env_param,fullfile(env_path,env_file),env_title,root_finder);
    
    % run Kraken
    cmdstr = sprintf('%s %s',runkraken,env_file);
    system(cmdstr);

    % load modes
    mod_file = [env_file,'.mod'];
    modes = read_modes(fullfile(env_path,mod_file));
%     imag_kr_diff = diff(imag(modes.k));
%     idx = find(imag_kr_diff(2:end)./imag_kr_diff(1:end-1)>10,1);
%     mode_incl = idx+1;
    mode_incl = length(modes.k);
    
    % get field from modes
    pp_freq(iF,:,:) = get_pfield(modes,mode_incl,r,z0,z);
end
toc

% Transmit signal spectrum
% ftau = 100;
% S = exp(-(freq_all-mean(freq_all)).^2/ftau^2)';

% Assemble complex field for each field point
P = reshape(pp_freq,length(freq_all),[]);
% PS = P.*repmat(S,1,size(P,2));
PS = P;
PS_cmplx = [flipud(conj(PS));zeros(1,size(PS,2));PS];
p = ifft(PS_cmplx);
p2d = reshape(p,[],length(z),length(r));

% % Plot
% fig = figure;
% for iF=1:size(p2d,1)
%     figure(fig)
%     imagesc(10*log10(abs(squeeze(p2d(iF,:,:)))));
%     colormap(jet)
%     colorbar
%     pause(0.1)
% end







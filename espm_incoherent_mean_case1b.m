% 2016 06 06  Incoherent mean over water column, EPSM/TREX case Ib

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

if ~kraken_on_path   % add path if Kraken is not already on path
    addpath(genpath('F:\Dropbox\0_APL_normal_mode\kraken'));
end

addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
base_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';

mode_path = sprintf('%s/%s',base_path,'20160606_ESPM');

r = 0:10:10e3;  % distance [m]
z = 0.1:0.1:19;   % receiver depth [m]
z0 = 17.8;       % source depth [m]

mode_file1 = 'espm_case1b.mod';
modes1 = read_modes(fullfile(mode_path,mode_file1));
% mode_file2 = 'espm_case2b.mod';
% modes2 = read_modes(fullfile(mode_path,mode_file2));
% mode_file3 = 'espm_case3b.mod';
% modes3 = read_modes(fullfile(mode_path,mode_file3));

p1 = get_pfield(modes1,length(modes1.k),r,z0,z);
p1_traponly = get_pfield(modes1,32,r,z0,z);
p1_15modes = get_pfield(modes1,15,r,z0,z);
% p2 = get_pfield(modes2,length(modes2.k),r,z0,z);
% p3 = get_pfield(modes3,length(modes3.k),r,z0,z);

mean_incoh1 = mean(abs(p1).^2,1);
mean_incoh1_traponly = mean(abs(p1_traponly).^2,1);
mean_incoh1_15modes = mean(abs(p1_15modes).^2,1);

% mean_incoh2 = mean(abs(p2).^2,1);
% mean_incoh3 = mean(abs(p3).^2,1);

figure
plot(r/1e3,10*log10(mean_incoh1),'linewidth',2);
hold on
plot(r/1e3,10*log10(mean_incoh1_traponly),'linewidth',2);
plot(r/1e3,10*log10(mean_incoh1_15modes),'linewidth',2);
% plot(r/1e3,10*log10(mean_incoh2),'linewidth',2);
% plot(r/1e3,10*log10(mean_incoh3),'linewidth',2);
set(gca,'fontsize',14)
grid
legend('Case Ib, all modes','Case Ib, trap modes only','Case Ib, 15 modes')
% legend('Case Ib','Case IIb','Case IIIb')
title('Krakenc results');
xlabel('Range (km)');
ylabel('<|p|^2>, dB')



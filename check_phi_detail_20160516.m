% 2016 03 24  Compare the mode location due to variation of NMESH

mode_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\\multifreq_mode_calc_20160427';
% mode_file = 'trex13_2layers_5D_11401z.mod';
mode_file = 'mfenv_3000.0.mod';

clear read_modes_bin % to force rewind to beginning of mode file

modes = read_modes(fullfile(mode_path,mode_file));
phi = get_component(modes,'N');

phi_diff = diff(phi,1,1);

% Plot all modes
figure;
plot(real(modes.k),imag(modes.k),'.');
hold on
text(double(real(modes.k)),double(imag(modes.k)),num2str([1:length(modes.k)]'));
xlabel('Real(k)')
ylabel('Imag(k)')

% Get excitation amplitude across modes
src_depth = 10;
[~,z_idx] = min(abs(src_depth-modes.z));
exc_amp = real(modes.phi(z_idx,:));

figure
plot(1:length(exc_amp),exc_amp,'.-');
hold on
plot(1:length(exc_amp),abs(exc_amp),'.-');
xlabel('Mode number');
ylabel('Mode amplitude (real part)');

% modes_wanted = [10:20];
% % modes_wanted = 1:5:50;
% 
% figure
% corder = get(gca,'colororder');
% for iM=1:length(modes_wanted)
%     mnum = modes_wanted(iM);
%     subplot(1,length(modes_wanted),iM)
%     plot(imag(phi(:,mnum)),modes.z,'-','color',corder(2,:));
%     hold on
%     plot(real(phi(:,mnum)),modes.z,'-','color',corder(1,:));
%     axis ij
%     title(['Mode ',num2str(mnum)])
%     grid
%     if iM==1
%         L = legend('Imag','Real','location','south');
%         set(L,'fontsize',12)
%     end
% end


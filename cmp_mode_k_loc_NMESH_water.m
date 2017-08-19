% 2016 03 24  Compare the mode location due to variation of NMESH

folder = 'F:\Dropbox\0_APL_normal_mode\oalib.hlsresearch.com\Modes\AcousticsToolbox\atWin\at_Wintel\tests\wjlee_tests';

fpre = 'trex13_2layers';

sediment_H = 5;
nmesh = [500,800,1000,1500,1901,2500];

fig = figure;
corder = get(gca,'colororder');
for iN=1:length(nmesh)
    
    clear read_modes_bin % to force rewind to beginning of mode file
    
    fname = sprintf('%s_%dD_allModes_%dw.mod',fpre,sediment_H,nmesh(iN));
    modes = read_modes(fullfile(folder,fname));
    phi = get_component(modes,'N');

    fprintf('File: %s\n',fname);
    fprintf('Total number of modes: %d\n',length(modes.k));
    
    figure(fig)
    plot(real(modes.k),imag(modes.k)-0.005*(iN-1),'.');
    hold on
    
    k = modes.k;
    z = modes.z;
    save_fname = sprintf('%s_%dD_allModes_%dw_modeinfo.mat',fpre,sediment_H,nmesh(iN));
    save(fullfile(folder,save_fname),'k','phi','z');
end

figure(fig)
delta_mesh = 19./(nmesh-1);
L = legend(num2str(delta_mesh'));
set(L,'location','southeast','fontsize',12);
grid on
xlabel('Real(k)')
ylabel('Imag(k)')
ylim([min(ylim) 0.01])


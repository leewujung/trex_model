% 2016 03 18  Save mode shape for comparison

folder = 'F:\Dropbox\0_APL_normal_mode\oalib.hlsresearch.com\Modes\AcousticsToolbox\atWin\at_Wintel\tests\wjlee_tests';

fpre = 'trex13_2layers_';

sediment_H = [5];

% fig2 = figure;
% corder = get(gca,'colororder');
for iS=1:length(sediment_H)
    
    clear read_modes_bin % to force rewind to beginning of mode file
    
    fname = [fpre,num2str(sediment_H(iS)),'D_c.mod'];
    str = strsplit(fname(1:end-4),'trex13_2layers_');
    ss{iS} = str{2}(1:end-2);
    
    modes = read_modes(fullfile(folder,fname));
    phi = get_component(modes,'N');

    fprintf('File: %s\n',fname);
    fprintf('Total number of modes: %d\n',length(modes.k));
    
%     figure(fig2)
%     plot(real(modes.k),imag(modes.k)-0.005*(iS-1),'.','markersize',10);
%     hold on
    
    k = modes.k;
    z = modes.z;
    save_fname = [fpre,num2str(sediment_H(iS)),'D_modeshape.mat'];
    save(fullfile(folder,save_fname),'k','phi','z');
end

% figure(fig2)
% L = legend(ss);
% set(L,'location','northwest','fontsize',12);
% grid on
% xlabel('Real(k)')
% ylabel('Imag(k)')
% ylim([min(ylim) 0.01])
% 

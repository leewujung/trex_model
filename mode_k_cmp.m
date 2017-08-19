% 2016 03 17  Compare modes k

folder = 'F:\Dropbox\0_APL_normal_mode\oalib.hlsresearch.com\Modes\AcousticsToolbox\atWin\at_Wintel\tests\wjlee_tests';

fpre = 'trex13_2layers_';

sediment_H = [2,3,3.5,4,5];

% fig = figure;
fig2 = figure;
corder = get(gca,'colororder');
% symbol_list = {'.','o','x','^','s'};
for iS=1:length(sediment_H)
    
    clear read_modes_bin % to force rewind to beginning of mode file
    
    fname = [fpre,num2str(sediment_H(iS)),'D_c.mod'];
    str = strsplit(fname(1:end-4),'trex13_2layers_');
    ss{iS} = str{2}(1:end-2);
    
    modes = read_modes(fullfile(folder,fname));
    
    fprintf('File: %s\n',fname);
    fprintf('Total number of modes: %d\n',length(modes.k));
    
%     figure(fig)
%     subplot(length(sediment_H),2,(iS-1)*2+1);  % all roots
%     plot(real(modes.k),imag(modes.k),'.','color',corder(2,:));
%     L = legend(sprintf('sediment thickness = %sD',num2str(sediment_H(iS))));
%     set(L,'fontsize',10,'location','northwest');
%     axis([6 15 -0.3 0.05])
%     grid
%     
%     subplot(length(sediment_H),2,(iS-1)*2+2);  % blow up propagating modes
%     plot(real(modes.k),imag(modes.k),'.','color',corder(2,:));
%     L = legend(sprintf('sediment thickness = %sD',num2str(sediment_H(iS))));
%     set(L,'fontsize',10,'location','northwest');
%     axis([13.5 14.3 -0.01 0.01])
%     grid
    
    figure(fig2)
    plot(real(modes.k),imag(modes.k)-0.005*(iS-1),'.','markersize',10);
    hold on
    
    k = modes.k;
    save_fname = [fpre,num2str(sediment_H(iS)),'D_c.mat'];
    save(fullfile(folder,save_fname),'k');
end

figure(fig2)
L = legend(ss);
set(L,'location','northwest','fontsize',12);
grid on
xlabel('Real(k)')
ylabel('Imag(k)')
ylim([min(ylim) 0.01])


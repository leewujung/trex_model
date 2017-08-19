function [modes1,modes2] = mode_cmp( filename1,filename2, modenum )
% Compare two sets of modes and mode functions computed by KRAKEN


clear read_modes_bin % to force rewind to beginning of mode file

modes1 = read_modes( filename1 );
fprintf('file1: total %d modes',length(modes1.k));

clear read_modes_bin % to force rewind to beginning of mode file
fprintf('\n')
modes2 = read_modes( filename2 );
fprintf('file2: total %d modes',length(modes2.k));

% extract the specfied component from the stress-displacement vector
phi1 = get_component( modes1, 'N' );
phi2 = get_component( modes2, 'N' );

fname1 = filename1(1:end-4);
fname1(fname1=='_') = ' ';
fname2 = filename2(1:end-4);
fname2(fname2=='_') = ' ';

% plots of the wavenumbers in the complex plane
figure
plot(real(modes1.k), imag(modes1.k),'o');
hold on
plot(real(modes2.k), imag(modes2.k),'.');
xlabel('real(k)');
ylabel('imag(k)');
L = legend(fname1,fname2,'location','northwest');
set(L,'fontsize',12);


% Plots for the modes
figure
corder = get(gca,'colororder');
subplot_num = length(modenum);
for iplot=1:length(modenum)
    subplot(1,subplot_num,iplot);
    plot(real(phi1(:,modenum(iplot))), modes1.z,'color',corder(1,:),'linewidth',2);
    hold on
    plot(imag(phi1(:,modenum(iplot))), modes1.z,'color',corder(1,:));
    plot(real(phi2(:,modenum(iplot))), modes2.z,'color',corder(2,:),'linestyle','--','linewidth',2);
    plot(imag(phi2(:,modenum(iplot))), modes2.z,'color',corder(2,:),'linestyle','--');
    axis ij
end



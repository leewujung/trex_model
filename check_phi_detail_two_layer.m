% 2016 03 24  Compare the mode location due to variation of NMESH

mode_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
% mode_file = 'trex13_2layers_5D_11401z.mod';
mode_file = 'two_layers_50m.mod';

clear read_modes_bin % to force rewind to beginning of mode file

modes = read_modes(fullfile(mode_path,mode_file));
phi = get_component(modes,'N');

phi_diff = diff(phi,1,1);

% Plot all modes
figure;
plot(real(modes.k),imag(modes.k),'.');
hold on
text(double(real(modes.k)),double(-imag(modes.k)),num2str([1:length(modes.k)]'));
xlabel('Real(k)')
ylabel('Imag(k)')

figure;
plot(real(modes.k),imag(modes.k),'.');
xlabel('Real(k)')
ylabel('Imag(k)')


modes_wanted = [1:5];
% modes_wanted = 1:5:50;

fig = figure;
corder = get(gca,'colororder');
for iM=1:length(modes_wanted)
    mnum = modes_wanted(iM);
    figure(fig)
    subplot(1,length(modes_wanted),iM)
    plot(imag(phi(:,mnum)),modes.z,'-','color',corder(2,:));
    hold on
    plot(real(phi(:,mnum)),modes.z,'-','color',corder(1,:));
    axis ij
    title(['Mode ',num2str(mnum)])
    grid
    if iM==1
        L = legend('Imag','Real','location','south');
        set(L,'fontsize',12)
    end
    
    figure
    subplot(121)
    plot(real(phi(:,mnum)),modes.z,'-','color',corder(1,:));
    axis ij
    ylabel('Depth (m)')
    xlabel('Real part')
    title(['Mode ',num2str(mnum)])
    grid
    subplot(122)
    plot(imag(phi(:,mnum)),modes.z,'-','color',corder(1,:));
    axis ij
    ylabel('Depth (m)')
    xlabel('Imag part')
    title(['Mode ',num2str(mnum)])
    grid
end


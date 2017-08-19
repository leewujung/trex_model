% 2016 03 25  Calculate pressure field

clear

mode_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
mode_file = 'trex13_2layers_5D_allModes_1901w.mod';

clear read_modes_bin % to force rewind to beginning of mode file

modes = read_modes(fullfile(mode_path,mode_file));
phi = get_component(modes,'N');

cw = 1525;
kr = modes.k;

r = 0:1:500;  % distance [m]
z_calc = modes.z;     % depth [m]
delta_z = 0.1;  % depth spacing [m]
z = z_calc(1):delta_z:z_calc(end);
z0 = 10;       % source depth [m]

fig = figure;
mode_incl_all = 1:length(modes.k);
for iM=30%length(mode_incl_all)
    mode_incl = mode_incl_all(iM);
    phi_incl = phi(:,1:mode_incl);
    
    % interpolation on mode function
    phi_z = interp1(z_calc,phi_incl,z);
    if size(phi_z,1)~=length(z)
        phi_z = phi_z';
    end
    phi_z0 = interp1(z_calc,phi_incl,z0);
    
    kr_r = kr*r;
    H = besselh(0,2,kr_r);
    
    p = 1i*pi*repmat(phi_z0,size(phi_z,1),1).*phi_z *H(1:mode_incl,:);
    
    p_abs = abs(p);
    p_abs(p_abs<1e-7) = 1e-7;
    p_abs(isnan(p_abs)) = 1e-6;
    p_abs(isinf(p_abs)) = 1e-6;
    
    tl = -20*log10(p_abs);
    
    % Plot
    figure(fig)
    imagesc(r,z,tl);
    xlabel('Range (m)')
    ylabel('Depth (m)')
    colormap(fliplr(jet))
    caxisrev([30 80])
    title(['Modes included: ',num2str(mode_incl_all(iM))])
%     if iM<30
%         pause(1)
%     elseif iM==30
%         pause
%     else
%         pause(0.1)
%     end
    
end
pfcn = get_pfield(modes,mode_incl,r,z0,z);

% 2016 03 25  Calculate pressure field

clear

mode_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\multifreq_mode_calc_20160427';
fc = 3000;  % cetner frequency for the transmit signal [Hz]

mode_file = sprintf('mfenv_%06.1f.mod',fc);
modes = read_modes(fullfile(mode_path,mode_file));
% mode_file = 'trex13_2layers_5D_allModes_1901w.mod';

modes = read_modes(fullfile(mode_path,mode_file));
phi = get_component(modes,'N');

freq = modes.freq;
cw = 1525;
kw = 2*pi*freq/cw;
kr = modes.k;

r = 0:1:500;  % distance [m]
z_calc = modes.z;     % depth [m]
delta_z = 0.1;  % depth spacing [m]
z = z_calc(1):delta_z:z_calc(end);
z0_d = 0.2;  % element spacing [m]
z0_all = 8:z0_d:12;
theta = -20;  % steering angle relative to horizontal [deg]

strtime = datetime('now','format','yyyyMMdd');
save_path = sprintf('%s_pfield_steer_theta_%.1fdeg_zs_%.1f-%.1f-%.1fm',strtime,theta,z0_all(1),z0_d,z0_all(end));
if ~exist(fullfile(mode_path,save_path),'dir')
    mkdir(fullfile(mode_path,save_path));
end

fig = figure;
mode_incl_all = 1:length(modes.k);
for iM=1:30%length(mode_incl_all)
    % modes to be included
    mode_incl = mode_incl_all(iM);
    phi_incl = phi(:,1:mode_incl);
    
    % interpolation on mode function
    phi_z = interp1(z_calc,phi_incl,z);
    if size(phi_z,1)~=length(z)
        phi_z = phi_z';
    end
    
    for iN=1:length(z0_all)
        % source contribution
        z0 = z0_all(iN);
        steer_phase = exp(1i*kw*z0_d*(iN-1)*sin(theta/180*pi));
        phi_z0 = interp1(z_calc,phi_incl,z0)*steer_phase;
        
        kr_r = kr*r;
        H = besselh(0,2,kr_r);
        
        p_N(iN,:,:) = 1i*pi*repmat(phi_z0,size(phi_z,1),1).*phi_z *H(1:mode_incl,:);
    end
    p = sum(p_N,1);
    %     p_abs = abs(p);
    p_abs = abs(squeeze(p));
    icount = find(p_abs>1e-7);         % for stats, only these values count <-- taken from Kraken
    p_abs(p_abs<1e-7) = 1e-7;
    p_abs(isnan(p_abs)) = 1e-6;
    p_abs(isinf(p_abs)) = 1e-6;
    
    tl = -20*log10(p_abs);
    
    % compute some statistics to automatically set the color bar
    tlmed = median(tl(icount));	% median value
    tlstd = std(tl(icount));     % standard deviation
    tlmax = tlmed+0.75*tlstd;       % max for colorbar
    tlmax = 10*round(tlmax/10);     % make sure the limits are round numbers
    tlmin = tlmax-50;               % min for colorbar
    if iM==1
        tlmin_const = tlmin;
        tlmax_const = tlmax;
    else
        if tlmax<tlmax_const
            tlmax_const = tlmax;
            tlmin_const = tlmin;
        end
    end
    
    % Plot
    figure(fig)
    imagesc(r,z,tl);
    xlabel('Range (m)')
    ylabel('Depth (m)')
    colormap(fliplr(jet))
    caxisrev([tlmin_const tlmax_const])
    ttext = sprintf('\\theta=%.1fdeg, \\deltaz=%.2fm, modes included: %d',theta,z0_d,iM);
    title(ttext)
    ylim([0 25]);
%     pause(0.5)
    %     if iM<30
    %         pause(0.1)
    %     elseif iM==30
    %         pause
    %     else
    %         pause(0.1)
    %     end
    
    save_file = sprintf('pfield_mode_incl_%02d.png',iM);
    saveas(gcf,fullfile(mode_path,save_path,save_file));
    
end


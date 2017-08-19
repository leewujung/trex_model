function p = get_pfield(modes,mode_incl,r,z0,z)
% r             % distance [m]
% delta_z       % depth spacing [m]

kr = modes.k;

% interpolation on mode function
phi_incl = modes.phi(:,1:mode_incl);
phi_z = interp1(modes.z,phi_incl,z);
if size(phi_z,1)~=length(z)
    phi_z = phi_z.';
end

% source contribution
phi_z0 = interp1(modes.z,phi_incl,z0);

kr_r = kr*r;
H = besselh(0,2,kr_r);

p = 1i*pi*repmat(phi_z0,size(phi_z,1),1).*phi_z *H(1:mode_incl,:);

% p_abs = abs(p);
% p_abs(p_abs<1e-7) = 1e-7;
% p_abs(isnan(p_abs)) = 1e-6;
% p_abs(isinf(p_abs)) = 1e-6;
% 
% tl = -20*log10(p_abs);


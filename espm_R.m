% 2016 06 06  Plane wave reflection coefficients

% Reflection coeff
c1 = 1525;
rho1 = 1;
c2 = 1630;
rho2 = 2;
freq = 3500;

delta = 0.01;
k1 = 2*pi*freq/c1;
k2 = 2*pi*freq/c2*(1+1i*delta);
theta1 = 0:pi/360:pi/2;

k1z = k1*sin(theta1);
k2z = sqrt(k2.^2-k1.^2*cos(theta1).^2);

R = (rho2*k1z-rho1*k2z)./(rho2*k1z+rho1*k2z);

figure
plot(theta1/pi*180,abs(R));
grid

% % No loss case
% theta2 = acos(k1*cos(theta1)/k2);
% Z1 = rho1/c1./sin(theta1);
% Z2 = rho2/c2./sin(theta2);
% R = (Z2-Z1)./(Z2+Z1);

% Load mode function
base_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests';
mode_path = sprintf('%s/%s',base_path,'20160606_ESPM');
mode_file = 'espm_case1b.mod';
modes = read_modes(fullfile(mode_path,mode_file));

k_real = real(modes.k);
k_imag = imag(modes.k);
theta_k_real = acos(k_real/k1);

figure
plot(theta1/pi*180,abs(R));
hold on
plot(theta_k_real/pi*180,exp(+k_imag),'.');


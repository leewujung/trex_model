% 2016 03 31  TREX13 environmental parameters
% function env_param_trex13(out_path,out_fname)

clear

out_path = 'F:\Dropbox\0_APL_normal_mode\kraken\tests\wjlee_tests\20160401_time_ir';
out_fname = 'env_param_trex13.mat';
if ~exist(out_path,'dir')
    mkdir(out_path);
end

SSP{1}.z = [0, 19.0];
SSP{1}.cp = [1525, 1525];
SSP{1}.rho = 1;
SSP{1}.sigma = 0;
SSP{1}.attn_p = 0;

SSP{2}.z = [19.0, 57.0];
SSP{2}.cp = [1630, 1630];
SSP{2}.rho = 2;
SSP{2}.sigma = 0;
SSP{2}.attn_p = 0.545750541536736;   % attenuation [dB/lambda]

nmesh_delta = 0.01;  % nmesh resolution
top_opt = 'NVW ';
bot_opt = 'V';
bot_sigma = 0;
clim = [0, 3000];  % search limit of phase speed [m/s]
rmax = 0;
src_z = 17.8;
rcv_z_range = [0 SSP{end}.z(end)];
rcv_z_delta = 0.1;
rcv_z_num = rcv_z_range(2)/rcv_z_delta+1;

save(fullfile(out_path,out_fname));
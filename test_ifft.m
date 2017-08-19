% 2016 04 01  Test ifft

freq = 500;
fs = 10000;
t = 0:1/fs:0.005;
s = cos(2*pi*freq*t);
S = fft(s);

S_shift = fftshift(S);

s_direct_ifft = ifft(S);

S_pos_half = S(2:26);

S_cmplx = [0,S_pos_half,fliplr(conj(S_pos_half))];
s_cmplx_ifft = ifft(S_cmplx);

S_freq = linspace(0,fs,length(S));
S_freq = S_freq(1:round(length(S_freq)/2));
S_freq = [S_freq,-fliplr(S_freq(2:end))];

time_delta = 0.001;
phase_shift = exp(1j*2*pi*S_freq*time_delta);

s_direct_ifft_delay = ifft(S.*phase_shift);

%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 1000;
fs = 10000;
tt_delta = 1/fs;
tau = 1e-3;  % [sec]
s_len = 500;
s_t = (-s_len:s_len)*tt_delta;
s = exp(-s_t.^2/tau^2).*(exp(1i*2*pi*fc*s_t)+exp(-1i*2*pi*fc*s_t))*1/2;

S = fft(s);
s_direct_ifft = ifft(S);

S_freq = linspace(0,fs,length(S));
S_freq = S_freq(1:round(length(S_freq)/2));
S_freq = [S_freq,-fliplr(S_freq(2:end))];

time_delta = 0.01;
phase_shift = exp(-1j*2*pi*S_freq*time_delta);

s_direct_ifft_delay = ifft(S.*phase_shift);

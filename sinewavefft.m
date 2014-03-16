% Sine Wave Windowed FFT
clc;
clear all;
close all;

data = load('FFTOut.csv');
sinewave0 = load('sinewave0.csv');
snfft = data(1:256,:);

num_secs = 10;
a = 8;
b = 128;
t = 0:0.001:num_secs;
n = 0:num_secs*b;
total = num_secs*b;
size = 640;
fft_window = 256;
fft_inc = 16;

sn = sinewave0(2:end,3:end);

offset = (total/2)-(size/2);
k=1;
newsignal = zeros(24,fft_window);
for i=offset+1:fft_inc:offset+size-fft_window-1
    n=1;
    for j=i+1:i+fft_window-1
        newsignal(k,n)=sn(j);
        n=n+1;
    end
    k=k+1;
end
rotate = newsignal';
y = fft(rotate); 
snfft0 = snfft(:,1:2:end)+1j*(snfft(:,2:2:end));
z=y';
plot(1:fft_window, abs(z(1,:)));
figure();
w = snfft0';
plot(1:fft_window, abs(w(1,:)));

pxx = abs(w(1,:)).^2 * (1/128)^2 / 2;
figure();
plot(1:length(w(1,:)),pxx);

raw = load('raw.csv');
figure();
plot(1:fft_window,raw);


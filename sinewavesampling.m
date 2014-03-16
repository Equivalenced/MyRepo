% Sine Wave Sampling
clc;
clear all;
close all;

num_secs = 10;
a = 8;
b = 128;
t = 0:0.001:num_secs;
n = 0:num_secs*b;
signal = sin(a*pi*t);
sn = zeros(1,num_secs*b);
for i=0:num_secs*b
    sn(i+1) = sin((a/b)*pi*(i));
end
plot(t, signal);
figure();
stem(n, sn);
figure();
plot(n, abs(fft(sn)));
sn_fft = fft(sn);
format short;

dlmwrite('sn_fft.csv',fft(sn),'delimiter',',','precision','%.6f');
% dlmwrite('sinewave.csv',sn,'delimiter',',','precision','%.6f');
% DLMWRITE('file.txt',M,'delimiter','\t','precision','%.6f') writes M
%     to file file.txt with elements delimited by the tab character, using a
%     precision of 6 decimal places. 

sn = [123 sn];
emptyVector = zeros(1,length(sn));
numbers = 1:1:length(sn);
sn0 = [numbers' emptyVector' sn' sn' sn' sn' sn' sn' sn' sn' sn' sn' sn' sn' sn' sn'];
% Emotiv type formatting for .csv file
dlmwrite('sinewave0.csv',sn0,'delimiter',',','precision','%.6f');


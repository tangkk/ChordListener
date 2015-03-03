function [x, fs] = myInput(path)

[x,fs] = audioread(path);
DSR = fs / 11025;
x = (x(:,1)+x(:,2))/2;
x = resample(x, 1, DSR);
fs = fs / DSR;
songMax = max(abs(x));
x = x / songMax;



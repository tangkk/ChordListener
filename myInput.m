function [x, fs, songinfo, DSR] = myInput(path)

[x,fs] = audioread(path);
songinfo = audioinfo(path);
DSR = fs / 11025;
x = (x(:,1)+x(:,2))/2;
x = resample(x, 1, DSR);
fs = fs / DSR;
songMax = max(abs(x));
x = x / songMax;



% this implements the machine listening method for catching every chord
% and chord boundary within a piece of music
% note that the most part of this algorithm reimplements the work by
% M. Mauch - Mauch, M. (2010). Automatic chord transcription from audio using
% computational models of musical context (Doctoral dissertation,
% School of Electronic Engineering and Computer Science Queen Mary, University of London).

root = '../AudioSamples/';
audio = 'tuihou-1.mp3';
path = [root audio];

[x,fs] = audioread(path);
songinfo = audioinfo(path);
DSR = fs / 11025;
x = (x(:,1)+x(:,2))/2;
x = resample(x, 1, DSR);
fs = fs / DSR;
songMax = max(abs(x));
x = x/songMax;
len = length(x);

% player = audioplayer(song,fs);
% play(player);

% implement a size 4096 hamming windowed, hop size 512 STFT spectrogram
wl = 4096;
w = hamming(wl);
hopsize = 512;
X = zeros(wl/2, ceil((length(x)-wl)/hopsize));
idx = 1;
for i = wl/2+1:hopsize:len - wl/2
    raws = i-wl/2;
    rawe = i+wl/2-1;
    raw = x(raws:rawe).*w;
    fftraw = fft(raw);
    X(:,idx) = 2*abs(fftraw(1:wl/2));
    X(:,idx) = X(:,idx) / max(X(:,idx));
    idx = idx + 1;
end
f = fs/2*linspace(0,1,wl/2);
k = (1/fs)*(1:len);
image(k,f,X);
set(gca,'YDir','normal');


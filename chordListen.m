% this implements the machine listening method for catching every chord
% and chord boundary within a piece of music
% note that the first part of this algorithm reimplements the work by
% M. Mauch - Mauch, M. (2010). Automatic chord transcription from audio using
% computational models of musical context (Doctoral dissertation,
% School of Electronic Engineering and Computer Science Queen Mary, University of London).

root = '../AudioSamples/';
audio = 'tuihou-1.mp3';
path = [root audio];
close all;

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
    
%     % normalization
    X(:,idx) = X(:,idx) / max(X(:,idx));
    idx = idx + 1;
end
f = fs/2*linspace(0,1,wl/2);
k = (1/fs)*(1:len);
% image(k,f,X);
% set(gca,'YDir','normal');

% build the tone profiles for calculating note salience matrix
% each sinusoidal at frequency 'ftone' is generated via sin(2*pi*n*f/fs)
fmin = 27.5; % MIDI note 21
fmax = 1661; % MIDI note 92
fratio = 2^(1/36);
numtones = 72*3;
Ms = zeros(numtones, wl/2); % simple tone profiles
Mc = zeros(numtones, wl/2); % complex tone profiles
% the true frequency of the tone is supposed to lie on bin midinum*3-1
for toneidx = 1:1:numtones
    ftone = fmin*(fratio^(toneidx-2));
    stone = sin(2*pi*(1:wl)*ftone/fs).*w';
    ctone = (0.9*sin(2*pi*(1:wl)*ftone/fs) + 0.9^2*sin(2*pi*(1:wl)*2*ftone/fs) + ...
        0.9^3*sin(2*pi*(1:wl)*3*ftone/fs) + 0.9^4*sin(2*pi*(1:wl)*4*ftone/fs)).*w';
    
    ffttone = abs(fft(stone));
    ffttone = ffttone(1:wl/2);
    ffttone = ffttone / norm(ffttone,2);
    
    fftctone = abs(fft(ctone));
    fftctone = fftctone(1:wl/2);
    fftctone = fftctone / norm(fftctone,2);
    Ms(toneidx,:) = ffttone;
    Mc(toneidx,:) = fftctone;
end

% calculate note salience matrix of the stft spectrogram (cosine
% similarity)
Ss = Ms*X;
Sc = Mc*X;
sizeM = size(Ms);
sizeX = size(X);
% p = 1:sizeM(1);
% k = 1:sizeX(2);
% figure;
% image(k,p,Ss);
% set(gca,'YDir','normal');
% figure;
% image(k,p,Sc);
% set(gca,'YDir','normal');

% calculate running mean and runnind std matrix for every column
sizeSpre = size(Sc);
rmeanS = zeros(sizeSpre);
rstdS = zeros(sizeSpre);
rmeanC = zeros(sizeSpre);
rstdC = zeros(sizeSpre);
for j = 1:1:sizeSpre(2)
    % do a running mean and std of the 216 bins within a sliding window of <18 bin
    colS = Ss(:,j);
    colC = Sc(:,j);
    for i = 1:1:sizeSpre(1)
        wmean = max(i-8,1):min(i+9,sizeSpre(1));
        rmeanS(i,j) = mean(colS(wmean));
        rmeanC(i,j) = mean(colC(wmean));
%         rstdS(i,j) = std(colS(wmean));
%         rstdC(i,j) = std(colC(wmean));
    end
end

% compute preliminary salience matrix
Spre = Ss.*Sc;

for i = 1:1:sizeSpre(1)
    for j = 1:1:sizeSpre(2)
        if Ss(i,j) < rmeanS(i,j) || Sc(i,j) < rmeanC(i,j)
            Spre(i,j) = 0;
        end
    end
end
% p = 1:sizeSpre(1);
% k = 1:sizeSpre(2);
% figure;
% image(k,p,Spre);
% set(gca,'YDir','normal');
% title('preliminary salience matrix');

% % tuning algorithm (this algorithm seems to be unreliable)
% Sbar = (sum(Spre,2))/sizeSpre(2); % first sum over all time frames
% fftSbar = fft(Sbar, 2^nextpow2(sizeSpre(1)));
% lenfftSbar = length(fftSbar);
% tuningpoint = lenfftSbar/6;
% phi = angle(fftSbar(floor(tuningpoint)))*(tuningpoint - floor(tuningpoint))...
%     + angle(fftSbar(ceil(tuningpoint)))*(ceil(tuningpoint) - tuningpoint); % take the angle at pi/3
% d = - phi - 2*pi/3; % wrap d in [-pi,pi)
% if d < -pi
%     d = d + 2*pi;
% end
% if d > pi
%     d = d - 2*pi;
% end
% d = d / (2*pi);
% tau = 440*2^(d/12);
% if d < 0
%     for j = 1:1:sizeSpre(2)
%         col = Spre(:,j);
%         for i = 2:1:sizeSpre(1)
%             % do interpolation here
%             Spre(i,j) = -d*col(i-1) + (1+d)*col(i);
%         end
%     end
% else
%     for j = 1:1:sizeSpre(2)
%         col = Spre(:,j);
%         for i = 1:1:sizeSpre(1) - 1
%             % do interpolation here
%             Spre(i,j) = (1-d)*col(i) + d*col(i+1);
%         end
%     end
% end
sfactor = 1000;
S = zeros(sizeSpre(1)/3, sizeSpre(2));
for i = 1:3:sizeSpre(1)
    for j = 1:1:sizeSpre(2)
        S((i+2)/3,j) = sfactor*(Spre(i,j) + Spre(i+1,j) + Spre(i+2,j));
    end
end
sizeS = size(S);
p = 1:sizeS(1);
k = 1:sizeS(2);
figure;
image(k,p,S);
set(gca,'YDir','normal');
title('note salience matrix');

% compute bass and treble profiles (once for all time)
gb = zeros(1,sizeS(1));
gt = zeros(1,sizeS(1));
gw = ones(1,sizeS(1));
for i = 1:1:sizeS(1)
    midiidx = i + 21 - 1;
    if midiidx >= 21 && midiidx < 33
        % y = a + bx with two points (21,0) and (33,1)
        gb(i) = (1 - 33/12) + midiidx/12;
    end
    if midiidx >= 33 && midiidx <= 44
        gb(i) = 1;
    end
    if midiidx > 44 && midiidx <= 55
        % y = a + bx with two points (44,1) and (50,0.5)
        gb(i) = (1+44/12) - midiidx/12;
    end
    if midiidx > 55
        gb(i) = 0;
    end
end
for i = 1:1:sizeS(1)
    midiidx = i + 21 - 1;
    if midiidx < 45
        gt(i) = 0;
    end
    if midiidx >= 45 && midiidx < 56
        % y = a + bx with two points (50,0.5) and (56,1)
        gt(i) = (1/2 - 50/12) + midiidx/12;
    end
    if midiidx >= 56 && midiidx <= 68
        gt(i) = 1;
    end
    if midiidx > 68 && midiidx <= 92
        % y = a + bx with two points (68,1) and (92,0)
        gt(i) = 92/24 - midiidx/24;
    end
end

% compute the chromagram, bass chromagram and treble chromagram
chromagram = zeros(12, sizeS(2));
basschromagram = zeros(12, sizeS(2));
treblechromagram = zeros(12, sizeS(2));
for i = 1:1:12
    for k = 1:1:sizeS(1)/12
        chromagram(i,:) = chromagram(i,:) + S(i + 12*(k-1),:)*gw(i + 12*(k-1));
        basschromagram(i,:) = basschromagram(i,:) + S(i + 12*(k-1),:)*gb(i + 12*(k-1));
        treblechromagram(i,:) = treblechromagram(i,:) + S(i + 12*(k-1),:)*gt(i + 12*(k-1));
    end
end
notenames = {'A','A#','B','C','C#','D','D#','E','F','F#','G','G#'};
sizeCh = size(chromagram);
p = 1:sizeCh(1);
k = 1:sizeCh(2);
sfactor = 0.2;
figure;
image(k,p,sfactor*chromagram);
set(gca,'YDir','normal');
set(gca, 'YTick',1:12, 'YTickLabel', notenames);
title('chromagram');
figure;
image(k,p,sfactor*basschromagram);
set(gca,'YDir','normal');
set(gca, 'YTick',1:12, 'YTickLabel', notenames);
title('basschromagram');
figure;
image(k,p,sfactor*treblechromagram);
set(gca,'YDir','normal');
set(gca, 'YTick',1:12, 'YTickLabel', notenames);
title('treblechromagram');

% normalize the chromagrams
for i = 1:1:sizeS(2)
    chromagram(:,i) = chromagram(:,i) / max(chromagram(:,i));
    treblechromagram(:,i) = treblechromagram(:,i) / max(treblechromagram(:,i));
    basschromagram(:,i) = basschromagram(:,i) / max(basschromagram(:,i));
end

% try some DBN models for,
% chord progression inference,
% based on,
% chromagram observations.
%
% naming conventions:
% T = treble progression, B = bass progression,
% TC = treble chromagram, BC = bass chromagram (both Lmax normalized)
% Ch = chord progression(Ch(k) = T(k)/B(k), slash chord representation)
%
% the 1st model: 2 separate hidden markov chains
% treble structure:
% T(i-1) -> T(i)
%   |        |
%   v        v
% TC(i-1)   TC(i)
%
% bass structure:
% B(i-1) -> B(i)
%   |        |
%   v        v
% BC(i-1)   BC(i)
% 
% chord progression computation:
% Ch = T/B

% ******************************************* %
% treble model
ss = 2;
intra = zeros(ss);
intra(1,2) = 1;
inter = zeros(ss);
inter(1,1) = 1;
hnodes = 1;
onodes = 2;
dnodes = 1;
cnodes = 2;
eclass1 = [1 2];
eclass2 = [3 2];
eclass = [eclass1 eclass2];
Q = 25; % consider the maj min treble model, there are totally 24 chords  + 1 no-chord
O = 12; % 12 pitch-classes, each treble has a mean and cov for each pitch class
ns = [Q O];
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'observed', onodes, 'eclass1', eclass1, 'eclass2', eclass2);
% set CPDs: CPD{1} <- prior; CPD{2} <- emission CPD{3} <- transition
% let's order the trebles in such way:
% 1 2  3 4  5 6 7  8 9  10 11 12 13 14  15 16  17 18 19  20 21  22 23  24 25
% C C# D D# E F F# G G# A  A# B  Cm C#m Dm D#m Em Fm F#m Gm G#m Am A#m Bm N
% let's order the pitch classes in such way:
% 1 2  3 4  5 6 7  8 9  10 11 12
% C C# D D# E F F# G G# A  A# B
prior = normalise(ones(1,Q));
mu = zeros(O,Q);
sigma = zeros(O,O,Q);
for i = 1:1:12
    muimaj = zeros(1,12);
    muimin = zeros(1,12);
    muimaj(i) = 1;
    muimaj(mod(i+4-1,12)+1) = 1;
    muimaj(mod(i+7-1,12)+1) = 1;
    muimin(i) = 1;
    muimin(mod(i+3-1,12)+1) = 1;
    muimin(mod(i+7-1,12)+1) = 1;
    mu(:,i) = muimaj;
    mu(:,i+12) = muimin;
end
mu(:,Q) = ones(1,12);
for i = 1:1:Q
    sigma(:,:,i) = diag(ones(1,12)*0.2);
end
transmat = ones(Q,Q);
st = 10; % the self transition factor, with larger value yields stronger smoothy.
for i = 1:1:Q
    transmat(i,i) = transmat(i,i)*st;
end
transmat = mk_stochastic(transmat);

bnet.CPD{1} = tabular_CPD(bnet, 1, prior); % uniform distribution for prior
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', mu, 'cov', sigma); % gaussian emission probs
bnet.CPD{3} = tabular_CPD(bnet, 3, transmat); % a reasonable transition matrix

jengine = smoother_engine(jtree_2TBN_inf_engine(bnet));
T = sizeS(2); % the total number of time slices contained in the evidence
evidence = cell(ss,T);
% use real evidence from chromagram
for i = 1:1:T
    ev = treblechromagram(:,i);
    ev = [ev(4:end) ; ev(1:3)]; % rearrange so that C is at the bottom
    evidence(onodes,i) = num2cell(ev,1);
end
[jengine,llj] = enter_evidence(jengine, evidence);
mpe = find_mpe(jengine, evidence);
display(mpe(1,1:50));







% CPRSB - Chord Progression Recognition System B
% This implements a machine listening method for recognizing the chords
% and chord progression boundaries within a piece of music

% ********************************************************** %
% ********************* Input ****************************** %
% ********************************************************** %
close all;
clear;
clc;
root = '../AudioSamples/';
audio = 'heidongli-1.mp3';
path = [root audio];

% ********************************************************** %
% ********************* Front End ************************** %
% ********************************************************** %
display('frontend -- from input to harmonic salience matrix');
[x,fs] = audioread(path);
songinfo = audioinfo(path);
DSR = fs / 11025;
x = (x(:,1)+x(:,2))/2;
x = resample(x, 1, DSR);
fs = fs / DSR;
songMax = max(abs(x));
x = x/songMax;
len = length(x);

% implement a size 4096 hamming windowed, hop size 512 STFT spectrogram
% with sample rate 11025, each hop's duration is 0.0464s = 46.4399ms
wl = 4096;
w = hamming(wl);
hopsize = 512;
lenSlice = ceil((length(x)-wl)/hopsize);
X = zeros(wl/2, lenSlice);
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
kf = (1/fs)*(1:len);
sfactor = 100;
figure;
image(kf,f,sfactor*X);
set(gca,'YDir','normal');
title('spectrogram');

% build the tone profiles for calculating note salience matrix
% each sinusoidal at frequency 'ftone' is generated via sin(2*pi*n*f/fs)
fmin = 27.5; % MIDI note 21
fmax = 1661; % MIDI note 92
fratio = 2^(1/36);
numtones = 72*3;
Ms = zeros(numtones, wl/2); % simple tone profiles
Mc = zeros(numtones, wl/2); % complex tone profiles
% the true frequency of the tone is supposed to lie on bin notenum*3-1,
% e.g. A4 is bin 49*3-1 = 146, C4 is bin 40*3-1 = 119 (note that notenum is
% not midinum, note num is the index of the key on a piano with A0 = 1)
staticbassbound = 30;
statictreblebound = 60;
bassboot = 1.5;
trebleboot = 0.8;
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
    % bass and treble boot
    if toneidx < staticbassbound*3
        ffttone = ffttone * bassboot;
        fftctone = fftctone * bassboot;
    end
    if toneidx > statictreblebound*3
        ffttone = ffttone * trebleboot;
        fftctone = fftctone * trebleboot;
    end
    Ms(toneidx,:) = ffttone;
    Mc(toneidx,:) = fftctone;
end

% calculate note salience matrix of the stft spectrogram (cosine
% similarity) (note that this is an additive approach, as contrast to
% the nnls approach which is an deductive approach)
Ss = Ms*X;
Sc = Mc*X;
sizeM = size(Ms);
sizeX = size(X);
ps = 1:sizeM(1);
ks = 1:sizeX(2);
sfactor = 100;
figure;
image(ks,ps,sfactor*Ss);
set(gca,'YDir','normal');
title('simple tone salience matrix');
figure;
image(ks,ps,sfactor*Sc);
set(gca,'YDir','normal');
title('complex tone salience matrix');

Spre = Sc.*Sc;
figure;
image(ks,ps,sfactor*Spre);
set(gca,'YDir','normal');
title('preliminary salience matrix');

% calculate spectral centroid
Sec = Spre(:,1:100);
SX = sum(Sec,2);
sc = round(sum(SX.*(1:length(SX))') / sum(SX));
scw = sc/(length(SX));

% noise reduction process
sizeNS = size(Sc);
nt = 0.1;
Ssn = zeros(sizeNS(1),sizeNS(2));
Scn = zeros(sizeNS(1),sizeNS(2));
for j = 1:1:sizeNS(2)
    Ss(:,j) = Ss(:,j) / max(Ss(:,j));
    Sc(:,j) = Sc(:,j) / max(Sc(:,j));
    vs = Ss(:,j);
    vc = Sc(:,j);
    [pkss,locss] = findpeaks(vs,'MinPeakHeight',nt,'MinPeakProminence',nt);
    [pksc,locsc] = findpeaks(vc,'MinPeakHeight',nt,'MinPeakProminence',nt);
    vs = zeros(sizeNS(1),1);
    vs(locss) = pkss;
    Ssn(:,j) = vs;
    vc = zeros(sizeNS(1),1);
    vc(locsc) = pksc;
    Scn(:,j) = vc;
end
sfactor = 100;
figure;
image(ks,ps,sfactor*Ssn);
set(gca,'YDir','normal');
title('noised reduced simple tone salience matrix');
figure;
image(ks,ps,sfactor*Scn);
set(gca,'YDir','normal');
title('noised reduced complex tone salience matrix');

Spres = zeros(sizeNS(1)/3, sizeNS(2));
Sprec = zeros(sizeNS(1)/3, sizeNS(2));
for i = 1:3:sizeNS(1)
    for j = 1:1:sizeNS(2)
        Spres((i+2)/3,j) = (Ssn(i,j) + Ssn(i+1,j) + Ssn(i+2,j));
        Sprec((i+2)/3,j) = (Scn(i,j) + Scn(i+1,j) + Scn(i+2,j));
    end
end
S = Spres.*Sprec;
sizeS = size(S);
for j = 1:1:sizeNS(2)
    S(:,j) = S(:,j) / max(S(:,j)); 
end
sfactor = 100;
p = 1:sizeS(1);
k = 1:sizeS(2);
figure;
image(k,p,sfactor*S);
set(gca,'YDir','normal');
title('note salience matrix');

% calculate spectral centroid
Sec = S(:,1:100);
SX = sum(Sec,2);
sc = round(sum(SX.*(1:length(SX))') / sum(SX));
scw = sc/(length(SX));

% input from above, if a piece of salience is shorter than a gestalt window, ignore it
wg = 10;
Sgneg = zeros(sizeS(1), sizeS(2));
for i = 1:1:sizeS(1)
    trackidx = 1;
    islight = 0;
    for j = 1:1:sizeS(2)
        if S(i,j) == 0
            if islight == 1;
                lenLight = j - trackidx;
                if lenLight <= wg
                    Sgneg(i,trackidx:j-1) = zeros(1,lenLight);
                end
            end
            trackidx = j;
            islight = 0;
        else
            Sgneg(i,j) = S(i,j);
            islight = 1;
        end
    end
end
figure;
image(k,p,sfactor*Sgneg);
set(gca,'YDir','normal');
title('note gestalt salience matrix - 1');

% if within a gestalt window ahead there's a non-zero bin, compensate the
% blank in the middle
wg = 20;
Sgpos = zeros(sizeS(1), sizeS(2));
for i = 1:1:sizeS(1)
    trackidx = 1;
    isblank = 0;
    for j = 1:1:sizeS(2)
        if Sgneg(i,j) > 0
            % compensate the gestalt
            if isblank == 1
                lenBlank = j - trackidx;
                if lenBlank <= wg
                    Sgpos(i,trackidx:j-1) = Sgneg(i,trackidx)*ones(1,lenBlank);
                end
            end
            Sgpos(i,j) = Sgneg(i,j);
            isblank = 0;
            trackidx = j;
        else
            isblank = 1;
        end
    end
end
figure;
image(k,p,sfactor*Sgpos);
set(gca,'YDir','normal');
title('note gestalt salience matrix - 2');

Sg = Sgpos; % gestalt salience matrix

% onset filter (roughly detect the note onsets)
So = zeros(sizeS(1), sizeS(2)); % onset matrix
ot = 0.2;
for i = 1:1:sizeS(1)
    for j = 1:1:sizeS(2)
        if j == 1
            hi = Sg(i,j) - 0;
        else
            hi = Sg(i,j) - Sg(i,j-1);
        end
        if hi > ot
            So(i,j) = hi;
        end
    end
end
figure;
image(k,p,sfactor*So);
set(gca,'YDir','normal');
title('onset matrix');

% bassline filter (roughly set the dynamic bass bounds)
Sb = zeros(1, sizeS(2)); % bassline vector
bt = 0.3;
for j = 1:1:sizeS(2)
    for i = 1:1:sizeS(1)
        if Sg(i,j) >= bt
            Sb(j) = i;
            break;
        end
    end
end
figure;
plot(1:length(Sb),Sb,'*');
title('rough bassline');

% harmonic change filter (detect harmonic change boundaries)
Sh = zeros(sizeS(1),sizeS(2)); % harmonic bounded salience matrix (one slice per col)
Shv = zeros(sizeS(1),sizeS(2)); % harmonic change matrix (one chord per col)
Shc = zeros(1,sizeS(2)); % harmonic change moments
ht = ot;
whs = 0;
whe = 0;
shidx = 1;
firsttime = 1;
for j = 1:1:sizeS(2)
    for i = 1:1:min(Sb(j)+5,sizeS(1))
        if (So(i,j) > ht && (j - whs > 10 || firsttime == 1)) || j == sizeS(2)
            if firsttime == 1
                firsttime = 0;
                whs = j;
                break;
            end
            % take the mean over the harmonic window in terms of row
            if j == sizeS(2)
                whe = j;
                wh = whs:whe;
            else
                whe = j-1;
                wh = whs:whe;
            end
            for ii = 1:1:sizeS(1)
                gesiiwh = mean(Sg(ii,wh));
                if gesiiwh > 0.10
                    Sh(ii,wh) = ones(1,length(wh))*gesiiwh;
                end
            end
            % normalize the content within harmonic window in terms of col
            for jj = whs:1:whe
                tmp = Sh(:,jj);
                if max(tmp) ~= 0
                    tmp = tmp / max(tmp);
                end
                tmp(tmp < 0.2) = 0;
                Sh(:,jj) = tmp;
            end
            % fill the harmonic change vector
            Shv(:,shidx) = Sh(:,whs);
            Shc(shidx) = whe;
            shidx = shidx + 1;
            whs = j;
            break;
        end
    end
end
nchords = shidx - 1;
Shc(shidx) = sizeS(2);
Shv = Shv(:,(1:nchords));
Shc = Shc(:,(1:nchords));
figure;
image(k,p,sfactor*Sh);
set(gca,'YDir','normal');
title('harmonic bounded salience matrix');
figure;
image(k(1:nchords),p,sfactor*Shv);
set(gca,'YDir','normal');
title('harmonic change matrix');
figure;
plot(1:length(Shc),Shc,'o');
title('haromonic change moments');

% ********************************************************** %
% ********************* Mid End - B************************* %
% ********************************************************** %
display('midend-B -- treblechromagram and basschromagram');
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
SS = Sh; % select the salience matrix to be used
chromagram = zeros(12, sizeS(2));
basschromagram = zeros(12, sizeS(2));
treblechromagram = zeros(12, sizeS(2));
for i = 1:1:12
    for kk = 1:1:sizeS(1)/12
        chromagram(i,:) = chromagram(i,:) + SS(i + 12*(kk-1),:)*gw(i + 12*(kk-1));
        basschromagram(i,:) = basschromagram(i,:) + SS(i + 12*(kk-1),:)*gb(i + 12*(kk-1));
        treblechromagram(i,:) = treblechromagram(i,:) + SS(i + 12*(kk-1),:)*gt(i + 12*(kk-1));
    end
end

% normalize the chromagrams and transform them to C base
for i = 1:1:sizeS(2)
    if max(chromagram(:,i)) ~= 0
        chromagram(:,i) = chromagram(:,i) / max(chromagram(:,i));
        tmp = chromagram(:,i);
        chromagram(:,i) = [tmp(4:end) ; tmp(1:3)];
    end
    if max(treblechromagram(:,i)) ~= 0
        treblechromagram(:,i) = treblechromagram(:,i) / max(treblechromagram(:,i));
        tmp = treblechromagram(:,i);
        treblechromagram(:,i) = [tmp(4:end) ; tmp(1:3)];
    end
    if max(basschromagram(:,i)) ~= 0
        basschromagram(:,i) = basschromagram(:,i) / max(basschromagram(:,i));
        tmp = basschromagram(:,i);
        basschromagram(:,i) = [tmp(4:end) ; tmp(1:3)];
    end
end

% low-cut various chromagrams
lowT = 0.1;
chromagram(chromagram < lowT) = 0;
treblechromagram(treblechromagram < lowT) = 0;
basschromagram(basschromagram < lowT) = 0;

% plot varoius chromagrams
notenames = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
sizeCh = size(chromagram);
pc = 1:sizeCh(1);
kc = 1:sizeCh(2);
T = sizeS(2); % the total number of time slices contained in the evidence
t = ((hopsize/fs)*(1:T));
plotsize = 1:T;
sfactor = 100;
figure;
image(t(plotsize),pc,sfactor*chromagram(:,plotsize));
set(gca,'YDir','normal');
set(gca, 'YTick',1:12, 'YTickLabel', notenames);
title('chromagram');
xlabel('time(s)');
ylabel('pitch class');
figure;
image(t(plotsize),pc,sfactor*basschromagram(:,plotsize));
set(gca,'YDir','normal');
set(gca, 'YTick',1:12, 'YTickLabel', notenames);
title('basschromagram');
xlabel('time(s)');
ylabel('pitch class');
figure;
image(t(plotsize),pc,sfactor*treblechromagram(:,plotsize));
set(gca,'YDir','normal');
set(gca, 'YTick',1:12, 'YTickLabel', notenames);
title('treblechromagram');
xlabel('time(s)');
ylabel('pitch class');

% ********************************************************** %
% ********************* Back End - B************************ %
% ********************************************************** %
display('backend-B -- dbn');
% try DBN models for,
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
% consider the maj min treble model, as well as sus2 and power treble model
Qt = 37;
% 12 pitch-classes, each treble has a mean and cov for each pitch class
O = 12;
ns = [Qt O];
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'observed', onodes, 'eclass1', eclass1, 'eclass2', eclass2);
% set CPDs: CPD{1} <- prior; CPD{2} <- emission; CPD{3} <- transition
% let's order the trebles in such way:
% 1 2  3 4  5 6 7  8 9  10 11 12 13 14  15 16  17 18 19  20 21  22 23  24
% C C# D D# E F F# G G# A  A# B  Cm C#m Dm D#m Em Fm F#m Gm G#m Am A#m Bm
% 25 26  27 28  29 30 31  32 33  34  35  36 (37)
% C2 C#2 D2 D#2 E2 F2 F#2 G2 G#2 A2  A#2 B2 N
% 37 38  39 40  41 42 43  44 45  46  47  48 49
% C5 C#5 D5 D#5 E5 F5 F#5 G5 G#5 A5  A#5 B5 N
% let's order the pitch classes in such way:
% 1 2  3 4  5 6 7  8 9  10 11 12
% C C# D D# E F F# G G# A  A# B

% prior probabilities
prior = normalise(ones(1,Qt));
% emission probabilities
mu = zeros(O,Qt);
sigma = zeros(O,O,Qt);
for i = 1:1:12
    muimaj = zeros(1,12);
    muimin = zeros(1,12);
    mui2 = zeros(1,12);
    mui5 = zeros(1,12);
    
    % major treble
    muimaj(mod(i+0-1,12)+1) = 1;
    muimaj(mod(i+4-1,12)+1) = 1;
    muimaj(mod(i+7-1,12)+1) = 1;
    
    % minor treble
    muimin(mod(i+0-1,12)+1) = 1;
    muimin(mod(i+3-1,12)+1) = 1;
    muimin(mod(i+7-1,12)+1) = 1;
    
    % sus2 treble
    mui2(mod(i+0-1,12)+1) = 1;
    mui2(mod(i+2-1,12)+1) = 1;
    mui2(mod(i+7-1,12)+1) = 1;
    
    % power treble
    mui5(mod(i+0-1,12)+1) = 1;
    mui5(mod(i+7-1,12)+1) = 1;
    
    mu(:,i) = muimaj;
    mu(:,i+12) = muimin;
    mu(:,i+24) = mui2;
%     mu(:,i+36) = mui5;
end
mu(:,Qt) = ones(1,12); % the last one is no-treble
for i = 1:1:Qt
    sigma(:,:,i) = diag(ones(1,12))*0.2;
end

% transition probabilities
transmat = ones(Qt,Qt);
st = 500; % the self transition factor, with larger value yields stronger smoothy.
for i = 1:1:Qt
    transmat(i,i) = transmat(i,i)*st;
end
transmat = mk_stochastic(transmat);

bnet.CPD{1} = tabular_CPD(bnet, 1, prior); % uniform distribution for prior
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', mu, 'cov', sigma); % gaussian emission probs
bnet.CPD{3} = tabular_CPD(bnet, 3, transmat); % a reasonable transition matrix

jengine = smoother_engine(jtree_2TBN_inf_engine(bnet));
evidence = cell(ss,T);
% use real evidence from chromagram
for i = 1:1:T
    ev = treblechromagram(:,i);
    evidence(onodes,i) = num2cell(ev,1);
end
[jengine,llj] = enter_evidence(jengine, evidence);
mpe = find_mpe(jengine, evidence);

treblenames = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B',...
    'Cm','C#m','Dm','D#m','Em','Fm','F#m','Gm','G#m','Am','A#m','Bm',...
    'C2','C#2','D2','D#2','E2','F2','F#2','G2','G#2','A2','A#2','B2',...
    'N'};
%     'C5','C#5','D5','D#5','E5','F5','F#5','G5','G#5','A5','A#5','B5',...

figure;
plot(t(plotsize),cell2mat(mpe(1,plotsize)));
title('treble progression');
xlabel('time(s)');
ylabel('treble');
set(gca, 'YTick',1:length(treblenames), 'YTickLabel', treblenames);
display('treble mpe');
display(mpe(1,plotsize));

% ******************************************* %
% bass model
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
% consider the bass C, C#, ... B as well as N
Qb = 13;
% 12 pitch-classes, each treble has a mean and cov for each pitch class
O = 12;
ns = [Qb O];
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'observed', onodes, 'eclass1', eclass1, 'eclass2', eclass2);
% set CPDs: CPD{1} <- prior; CPD{2} <- emission; CPD{3} <- transition
% let's order the basses in such way:
% 1 2  3 4  5 6 7  8 9  10 11 12 13
% C C# D D# E F F# G G# A  A# B  N
% prior probabilities
prior = normalise(ones(1,Qb));
% emission probabilities
mu = zeros(O,Qb);
sigma = zeros(O,O,Qb);
for i = 1:1:12
    mubass = zeros(1,12);
    mubass(i) = 1;
    mu(:,i) = mubass;
end
mu(:,Qb) = ones(1,12); % the last one is no-bass

for i = 1:1:Qb
    sigma(:,:,i) = diag(ones(1,12))*0.2;
end

% transition probabilities
transmat = ones(Qb,Qb);
st = 500; % the self transition factor, with larger value yields stronger smoothy.
for i = 1:1:Qb
    transmat(i,i) = transmat(i,i)*st;
end
transmat = mk_stochastic(transmat);

bnet.CPD{1} = tabular_CPD(bnet, 1, prior); % uniform distribution for prior
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', mu, 'cov', sigma); % gaussian emission probs
bnet.CPD{3} = tabular_CPD(bnet, 3, transmat); % a reasonable transition matrix

jengine = smoother_engine(jtree_2TBN_inf_engine(bnet));
evidence = cell(ss,T);
% use real evidence from chromagram
for i = 1:1:T
    ev = basschromagram(:,i);
    evidence(onodes,i) = num2cell(ev,1);
end
[jengine,llj] = enter_evidence(jengine, evidence);
mpeb = find_mpe(jengine, evidence);

bassnames = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B',...
    'N'};

figure;
plot(t(plotsize),cell2mat(mpeb(1,plotsize)));
title('bass progression');
xlabel('time(s)');
ylabel('bass');
set(gca, 'YTick',1:length(bassnames), 'YTickLabel', bassnames);
display('bass mpe');
display(mpeb(1,plotsize));

% ******************************************* %
% chord model
lenChord = length(mpe);
chordprogression = cell(1,lenChord);
treblebassprogression = cell(2,lenChord);
for i = 1:1:lenChord
    % simply print slash chords except for maj, min
    tnum = mpe{1,i};
    bnum = mpeb{1,i};
    treblebassprogression{1,i} = num2treble(tnum);
    treblebassprogression{2,i} = num2bass(bnum);
    % map sus2 to maj
    if tnum >= 25 && tnum <= 36
        tnum = tnum - 24;
    end
    tname = num2treble(tnum);
    bname = num2bass(bnum);
    if tnum == Qt || bnum == Qb
        chordprogression{1,i} = 'N';
    elseif mod(tnum - 1,12) + 1 == bnum
        chordprogression{1,i} = tname;
    else
        chordprogression{1,i} = [tname '/' bname];
    end
end
display('chord progression');
display(chordprogression(1,plotsize));

% smooth the chord progression
lenChordPrint = 100;
chordprint = cell(2,lenChordPrint);
oldchord = chordprogression{1,1};
chordcount = 1;
chordidx = 1;
st = 10;
for i = 2:1:lenChord
    newchord = chordprogression{1,i};
    if strcmp(newchord,oldchord) == 1
        chordcount = chordcount + 1;
        if i == lenChord
            chordprint(:,chordidx) = {oldchord, chordcount};
        end
    else
        if chordcount <= st && chordidx > 1
            tmp = chordprint(:,chordidx - 1);
            tmp{2} = tmp{2} + chordcount;
            chordprint(:,chordidx - 1) = tmp;
            chordcount = 1;
        else
            chordprint(:,chordidx) = {oldchord, chordcount};
            chordidx = chordidx + 1;
            chordcount = 1;
        end
    end
    oldchord = newchord;
end
display('chordprint');
display(chordprint);

% write results
fw = fopen([audio(1:end-4) '.dbn.txt'],'w');
formatSpec1 = '%s';
formatSpec2 = '%s\n';
sumhop = 0;
for i = 1:1:lenChordPrint
    if ~isempty(chordprint{2,i})
        if sumhop == 0
            s = [chordprint{1,i} '===>' num2str(0)];
        else
            s = [chordprint{1,i} '===>' num2str(t(sumhop))];
        end
        fprintf(fw, formatSpec1, s);
        sumhop = sumhop + chordprint{2,i};
        s = ['-' num2str(t(sumhop))];
        fprintf(fw, formatSpec2, s);
    else
        break;
    end
end
fclose(fw);

% ********************* End of System B ******************** %
display('end of system B...');
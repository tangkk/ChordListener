% this implements the machine listening method for catching every chord
% and chord boundary within a piece of music
% note that part of this algorithm reimplements the work by
% M. Mauch - Mauch, M. (2010). Automatic chord transcription from audio using
% computational models of musical context (Doctoral dissertation,
% School of Electronic Engineering and Computer Science Queen Mary, University of London).

% ********************************************************** %
% ********************* Input ****************************** %
% ********************************************************** %
root = '../AudioSamples/';
audio = 'tuihou-1.mp3';
path = [root audio];
close all;

% ********************************************************** %
% ********************* Front End *************************** %
% ********************************************************** %
display('frontend');
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
bassbound = 28;
treblebound = 60;
bassboot = 2;
trebleboot = 0.5;
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
    % bass boot
    if toneidx < bassbound*3
        ffttone = ffttone * bassboot;
        fftctone = fftctone * bassboot;
    end
    if toneidx > treblebound*3
        ffttone = ffttone * trebleboot;
        fftctone = fftctone * trebleboot;
    end
    Ms(toneidx,:) = ffttone;
    Mc(toneidx,:) = fftctone;
end

% calculate note salience matrix of the stft spectrogram (cosine
% similarity)
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

% if within a gestalt window ahead there's a non-zero bin, compensate the
% blank in the middle
wg = 20;
Sg = zeros(sizeS(1), sizeS(2));
for i = 1:1:sizeS(1)
    trackidx = 1;
    isblank = 0;
    for j = 1:1:sizeS(2)
        if S(i,j) > 0
            % compensate the gestalt
            if isblank == 1
                lenBlank = j - trackidx;
                if lenBlank <= wg
                    Sg(i,trackidx:j-1) = S(i,trackidx)*ones(1,lenBlank);
                end
            end
            Sg(i,j) = S(i,j);
            isblank = 0;
            trackidx = j;
        else
            isblank = 1;
        end
    end
end
figure;
image(k,p,sfactor*Sg);
set(gca,'YDir','normal');
title('note gestalt salience matrix - 1');
% input from above, if a piece of salience is shorter than a gestalt window, ignore it
wg = 10;
for i = 1:1:sizeS(1)
    trackidx = 1;
    islight = 0;
    for j = 1:1:sizeS(2)
        if Sg(i,j) == 0
            if islight == 1;
                lenLight = j - trackidx;
                if lenLight <= wg
                    Sg(i,trackidx:j-1) = zeros(1,lenLight);
                end
            end
            trackidx = j;
            islight = 0;
        else
            islight = 1;
        end
    end
end
figure;
image(k,p,sfactor*Sg);
set(gca,'YDir','normal');
title('note gestalt salience matrix - 2');

% onset filter (roughly detect the note onsets)
So = zeros(sizeS(1), sizeS(2));
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

% harmonic change filter (detect harmonic change boundaries)
Sh = zeros(sizeS(1),sizeS(2));
Shv = zeros(sizeS(1),sizeS(2)); % harmonic change matrix (one chord per col)
Shc = zeros(1,sizeS(2)); % harmonic change moments
bassbound = 30;
ht = ot;
whs = 0;
whe = 0;
shidx = 1;
firsttime = 1;
for j = 1:1:sizeS(2)
    for i = 1:1:bassbound
        if (So(i,j) > ht && (j - whs > 10 || firsttime == 1)) || j == sizeS(2)
            if firsttime == 1
                firsttime = 0;
                whs = j;
                display(whs);
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
                tmp = tmp / max(tmp);
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

% ********************************************************** %
% ********************* Mid End *************************** %
% ********************************************************** %
display('midend -- various 12-bin grams or 1-bin grams');
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

% compute basegram and uppergram (based on harmonic change matrix)
sizeShv = size(Shv);
basegram = zeros(2,sizeShv(2));
uppergram = zeros(12,sizeShv(2));
for j = 1:1:sizeShv(2)
    for i = 1:1:sizeShv(1)
        % find out the lowest salience
        if Shv(i,j) > 0
            basegram(1,j) = mod((i+9)-1,12) + 1; % turn into a normal format
            basegram(2,j) = Shv(i,j);
            break;
        end
    end
end
for i = 1:1:12
    for kk = 1:1:sizeShv(1)/12
        uppergram(i,:) = uppergram(i,:) + Shv(i + 12*(kk-1),:);
    end
end
% normalize uppergram
for j = 1:1:sizeShv(2)
    if max(uppergram(:,j)) ~= 0
        uppergram(:,j) = uppergram(:,j) / max(uppergram(:,j));
        tmp = uppergram(:,j);
        uppergram(:,j) = [tmp(4:end) ; tmp(1:3)];
    end
end

notenames = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
ph = 1:sizeShv(2);
kh = 1:12;
sfactor = 100;
figure;
image(ph,kh,sfactor*uppergram);
set(gca,'YDir','normal');
set(gca, 'YTick',1:12, 'YTickLabel', notenames);
title('uppergram');
figure;
plot(ph,basegram(1,:),'o');
title('basegram');
set(gca, 'YTick',1:12, 'YTickLabel', notenames);

% ********************************************************** %
% ********************* Back End *************************** %
% ********************************************************** %

% try chord recognition tree method (or chord dictionary method)
% the chord tree is built this way with priority from top to down
% --music--            --digital--     --digital difference--
% 1->3->5 maj          1,5,8            4,7
% 1->3->5# aug         1,5,9            4,8
% 1->3->7b dom7        1,5,11           4,10
% 1->3->7 maj7         1,5,12           4,11
% 1->3b->5 min         1,4,8            3,7
% 1->3b->5b dim        1,4,7            3,6
% 1->3b->6b maj/3       1,4,9           3,8
% 1->3b->7b min7       1,4,11           3,10
% 1->3b->7 minmaj7     1,4,12           3,11
% 1->4->6 maj/5        1,6,10           5,9
% 1->4->5 sus4         1,6,8            5,7
% 1->2->5 sus2         1,3,8            2,7
% 1->3b min            1,4              3
% 1->3 maj             1,5              4
% 1->5 5               1,8              7
% 1->2 2               1,3              2
% default N            x,x,x            x,x
% run the bass through this tree, if hit, then that's the chord rooted
% on the bass
% if miss, then run every other pitch through this tree, pick a hit with
% root closest to the bass to form a slash chord
display('backend -- chordtree');
nchordtype = 16;
chordtree = cell(2,nchordtype);
chordtree{1,1} = [4,7];
chordtree{2,1} = 'maj';
chordtree{1,2} = [4,8];
chordtree{2,2} = 'aug';
chordtree{1,3} = [4,10];
chordtree{2,3} = 'dom7';
chordtree{1,4} = [4,11];
chordtree{2,4} = 'maj7';
chordtree{1,5} = [3,7];
chordtree{2,5} = 'min';
chordtree{1,6} = [3,6];
chordtree{2,6} = 'dim';
chordtree{1,7} = [3,8];
chordtree{2,7} = 'maj/3';
chordtree{1,8} = [3,10];
chordtree{2,8} = 'min7';
chordtree{1,9} = [3,11];
chordtree{2,9} = 'minmaj7';
chordtree{1,10} = [5,9];
chordtree{2,10} = 'maj/5';
chordtree{1,11} = [5,7];
chordtree{2,11} = 'sus4';
chordtree{1,12} = [2,7];
chordtree{2,12} = 'sus2';
chordtree{1,13} = 3;
chordtree{2,13} = 'min';
chordtree{1,14} = 4;
chordtree{2,14} = 'maj';
chordtree{1,15} = 7;
chordtree{2,15} = '5';
chordtree{1,16} = 2;
chordtree{2,16} = '2';
% walk the tree using basegram and uppergram
chordogram = cell(1,sizeShv(2));
for j = 1:1:sizeShv(2)
    bass = basegram(1,j);
    upper = uppergram(:,j);
    ismatchout = 0;
    for i = 1:1:nchordtype
        ismatchout = 0;
        ismatchin = 1;
        chordentry = chordtree{1,i};
        lenentry = length(chordentry);
        for jj = 1:1:lenentry
            matchpos = mod(bass+chordentry(jj)-1,12) + 1;
            if upper(matchpos) == 0
                ismatchin = 0;
            end
        end
        if ismatchin == 1
            ismatchout = 1;
            % convert maj/3 and maj/5 chord to correct bass
            if i == 7
                bass = mod(bass-4-1,12)+1;
            end
            if i == 10
                bass = mod(bass-7-1,12)+1;
            end
            chordogram{j} = [num2bass(bass) chordtree{2,i}];
            break;
        end
    end
    if ismatchout == 0
        chordogram{j} = [num2bass(bass) 'n'];
    end
end

% write results
fw = fopen([audio(1:end-4) '.ct.txt'],'w');
formatSpec1 = '%s';
formatSpec2 = '%s\n';
sumhop = 0;
lenchordogram = length(chordogram);
for i = 1:1:lenchordogram
    if i < lenchordogram
        if strcmp(chordogram{i}, chordogram{i+1}) == 1
            continue;
        end
    end
    if sumhop == 0
        s = [chordogram{i} '===>' num2str(0)];
    else
        s = [chordogram{i} '===>' num2str(t(sumhop))];
    end
    fprintf(fw, formatSpec1, s);
    sumhop = Shc(i);
    s = ['-' num2str(t(sumhop))];
    fprintf(fw, formatSpec2, s);
end
fclose(fw);

display('backend -- dbn');
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

% ******************************************* %
% play sound
% sound(x(1:hopsize*length(plotsize)),fs);

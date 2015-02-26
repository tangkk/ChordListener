% CPRSA - Chord Progression Recognition System A
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
x = x / songMax;
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
    if max(S(:,j)) > 0
        S(:,j) = S(:,j) / max(S(:,j)); 
    end
end
sfactor = 100;
p = 1:sizeS(1);
k = 1:sizeS(2);
figure;
image(k,p,sfactor*S);
set(gca,'YDir','normal');
title('note salience matrix');

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
% ********************* Mid End - A************************* %
% ********************************************************** %
display('midend-A -- uppergram and basegram');

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
        % if the whole col of Shv are zero, then basegram(:,j) is [0;0]
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

treblenotenames = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
ph = 1:sizeShv(2);
kh = 1:12;
sfactor = 100;
figure;
image(ph,kh,sfactor*uppergram);
set(gca,'YDir','normal');
set(gca, 'YTick',1:12, 'YTickLabel', treblenotenames);
title('uppergram');
bassnotenames = {'N','C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
figure;
plot(ph,basegram(1,:),'o');
title('basegram');
set(gca, 'YTick',0:12, 'YTickLabel', bassnotenames);

% ********************************************************** %
% ********************* Back End - A************************ %
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
display('backend-A -- chordtree');
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
chordogram = cell(3,sizeShv(2));
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
            chordogram{1,j} = bass;
            chordogram{2,j} = chordtree{2,i};
            chordogram{3,j} = num2bass(bass);
            break;
        end
    end
    if ismatchout == 0
        chordogram{1,j} = bass;
        chordogram{2,j} = 'n';
        chordogram{3,j} = num2bass(bass);
    end
end

% gestalt chordogram (merge chord types as human would do)
% only keep maj, min, aug, dim, dom, minmaj, /3 and /5
% if 2, 5, n type, merge to the nearest maj or min type with same bass
%   if there is no nearest maj or min, merge n to 2 or 5, if both 2 and 5
%   exist in the neigbourhood, merge them to be sus2
% if /3 and /5, merge from nearest whatever type with same bass
% if sus2, sus4, merge to maj with same bass
chordogram = [{0;'N';'N'} chordogram {0;'N';'N'}];
yes = 1;
while yes
    yes = 0;
    for i = 2:1:length(chordogram)-1
        cb = chordogram{1,i};
        ct = chordogram{2,i};
        pcb = chordogram{1,i-1};
        pct = chordogram{2,i-1};
        ncb = chordogram{1,i+1};
        nct = chordogram{2,i+1};
        if strcmp(ct,'n')
            if pcb == cb && ~(strcmp(pct,'n'))&& ~(strcmp(pct,'maj/3'))&& ~(strcmp(pct,'maj/5'))
                ct = pct;
                chordogram{2,i} = ct;
                yes = 1;
                continue;
            end
            if ncb == cb && ~(strcmp(nct,'n'))&& ~(strcmp(nct,'maj/3'))&& ~(strcmp(nct,'maj/5'))
                ct = nct;
                chordogram{2,i} = ct;
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'2')
            if pcb == cb && strcmp(pct,'5')
                ct = 'sus2';
                chordogram{2,i} = ct;
                yes = 1;
                continue;
            end
            if ncb == cb && strcmp(nct,'5')
                ct = 'sus2';
                chordogram{2,i} = ct;
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'5')
            if pcb == cb && strcmp(pct,'2')
                ct = 'sus2';
                chordogram{2,i} = ct;
                yes = 1;
                continue;
            end
            if ncb == cb && strcmp(nct,'2')
                ct = 'sus2';
                chordogram{2,i} = ct;
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'5') || strcmp(ct,'2')|| strcmp(ct,'n')
            if pcb == cb && ~(strcmp(pct,'5') || strcmp(pct,'2')|| strcmp(pct,'n')||strcmp(pct,'maj/3')||strcmp(pct,'maj/5'))
                ct = pct;
                chordogram{2,i} = ct;
                yes = 1;
                continue;
            end
            if ncb == cb && ~(strcmp(nct,'5') || strcmp(nct,'2')|| strcmp(nct,'n')||strcmp(nct,'maj/3')||strcmp(nct,'maj/5'))
                ct = nct;
                chordogram{2,i} = ct;
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'maj/3')
            if mod(pcb-4-1,12)+1 == cb && ~strcmp(pct,'maj/3') && ~strcmp(pct,'maj/5')
                chordogram{1,i-1} = cb;
                chordogram{2,i-1} = 'maj/3';
                chordogram{3,i-1} = num2bass(cb);
                yes = 1;
                continue;
            end
            if mod(ncb-4-1,12)+1 == cb && ~strcmp(nct,'maj/3') && ~strcmp(nct,'maj/5')
                chordogram{1,i+1} = cb;
                chordogram{2,i+1} = 'maj/3';
                chordogram{3,i+1} = num2bass(cb);
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'maj/5')
            if mod(pcb-7-1,12)+1 == cb && ~strcmp(pct,'maj/3') && ~strcmp(pct,'maj/5')
                chordogram{1,i-1} = cb;
                chordogram{2,i-1} = 'maj/5';
                chordogram{3,i-1} = num2bass(cb);
                yes = 1;
                continue;
            end
            if mod(ncb-7-1,12)+1 == cb && ~strcmp(nct,'maj/3') && ~strcmp(nct,'maj/5')
                chordogram{1,i+1} = cb;
                chordogram{2,i+1} = 'maj/5';
                chordogram{3,i+1} = num2bass(cb);
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'sus2') || strcmp(ct,'sus4')
            chordogram{2,i} = 'maj';
            yes = 1;
            continue;
        end
    end
end
chordogram = chordogram(:,2:end-1);

% write results
maxOutNum = 200;
outputstrings = cell(1,maxOutNum);

fw = fopen([audio(1:end-4) '.ct.txt'],'w');
formatSpec1 = '%s';
formatSpec2 = '%s\n';
sumhop = 0;
T = sizeS(2); % the total number of time slices contained in the evidence
t = ((hopsize/fs)*(1:T));
lenchordogram = length(chordogram);
chordboundaries = zeros(2,lenchordogram+1);
chordboundaries(2,:) = ones(1,lenchordogram+1);
outidx = 1;
for i = 1:1:lenchordogram
    outstring = '';
    if i < lenchordogram
        if strcmp(chordogram{3,i}, chordogram{3,i+1}) == 1  && strcmp(chordogram{2,i}, chordogram{2,i+1}) == 1
            continue;
        end
    end
    if sumhop == 0
        s = [chordogram{3,i} chordogram{2,i} '===>' num2str(0)];
    else
        s = [chordogram{3,i} chordogram{2,i} '===>' num2str(t(sumhop))];
    end
    fprintf(fw, formatSpec1, s);
    outstring = s;
    sumhop = Shc(i);
    s = ['-' num2str(t(sumhop))];
    if i < lenchordogram
        fprintf(fw, formatSpec2, s);
    else
        fprintf(fw, formatSpec1, s);
    end
    outstring = strcat(outstring,s);
    outputstrings{outidx} = outstring;
    outidx = outidx + 1;
    chordboundaries(2,outidx) = sumhop;
end
outputstrings = outputstrings(1:outidx-1);
chordboundaries = chordboundaries(:,1:outidx);
fclose(fw);

% visualize chord progression
lenOut = length(outputstrings);
chordprogression = cell(4,lenOut);

for i = 1:1:lenOut
    tline = outputstrings{i};
    s1 = strsplit(tline, '===>');
    chordname = s1{1};
    setime = s1{2};
    s2 = strsplit(setime, '-');
    starttime = s2{1};
    endtime = s2{2};
    if strcmp(chordname(2),'#')
        bassname = chordname(1:2);
        bassnum = bass2num(bassname);
        treblename = chordname(3:end);
    else
        bassname = chordname(1);
        bassnum = bass2num(bassname);
        treblename = chordname(2:end);
    end
    % transform /3 and /5 back to original bassname and bassnum
    if length(treblename) > 1
        if strcmp(treblename(end-1:end),'/3')
            bassnum = mod(bassnum+4-1,12) + 1;
            bassname = num2bass(bassnum);
        end
        if strcmp(treblename(end-1:end),'/5')
            bassnum = mod(bassnum+7-1,12) + 1;
            bassname = num2bass(bassnum);
        end
    end
    chordprogression{1,i} = chordname;
    chordprogression{2,i} = bassname;
    chordprogression{3,i} = treblename;
    chordprogression{4,i} = bassnum;
    chordboundaries(1,i) = str2double(starttime);
    chordboundaries(1,i+1) = str2double(endtime);
end
figure;
hold on;
Y = -10:0.1:10;
for i = 1:1:lenOut+1
    X = chordboundaries(1,i)*ones(size(Y));
    plot(X,Y);
end
hold off;
div = (max(Y) - min(Y) - 1) / lenOut;
for i = 1:1:lenOut
%     x = (chordboundaries(i)+ chordboundaries(i+1)) / 2;
    x = chordboundaries(1,i);
    text(x,10 - i*div,chordprogression{1,i});
end
xlabel('time');
ylabel('chord');
title('chordprogression vs. time');

figure;
hold on;
Y = -10:0.1:10;
for i = 1:1:lenOut+1
    X = chordboundaries(2,i)*ones(size(Y));
    plot(X,Y);
end
hold off;
div = (max(Y) - min(Y) - 1) / lenOut;
for i = 1:1:lenOut
%     x = (chordboundaries(i)+ chordboundaries(i+1)) / 2;
    x = chordboundaries(2,i);
    text(x,10 - i*div,chordprogression{1,i});
end
xlabel('slice');
ylabel('chord');
title('chordprogression vs. slices');

% ********************************************************** %
% ********************* Feedback Once - A******************* %
% ********************************************************** %

% feedback the chordboundaries information to note salience matrix
% and update the uppergram
display('feedback-A -- use chord boundaries information to do it again');
newbasegram = zeros(1,lenOut);
for i = 1:1:lenOut
    newbasegram(i) = chordprogression{4,i};
end

newuppergram = zeros(12,lenOut);
upg = zeros(12,1);
ut = 1;
for i = 1:1:lenOut
    % update note salience matrix in terms of boundaries window
    wb = chordboundaries(2,i):chordboundaries(2,i+1);
    sm = sum(S(:,wb),2);
    sm(sm < ut) = 0;
    if max(sm) > 0
        sm = sm ./ max(sm);
    end
    upg = sum(reshape(sm,12,6),2);
    upg = [upg(4:end);upg(1:3)];
    newuppergram(:,i) = upg;
end
treblenotenames = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
ph = 1:lenOut;
kh = 1:12;
sfactor = 100;
figure;
image(ph,kh,sfactor*newuppergram);
set(gca,'YDir','normal');
set(gca, 'YTick',1:12, 'YTickLabel', treblenotenames);
title('newuppergram');
bassnotenames = {'N','C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
figure;
plot(ph,newbasegram,'o');
title('newbasegram');
set(gca, 'YTick',0:12, 'YTickLabel', bassnotenames);

% ****** feedback back end ****** %
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
newchordogram = cell(3,lenOut);
for j = 1:1:lenOut
    bass = newbasegram(1,j);
    upper = newuppergram(:,j);
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
            newchordogram{1,j} = bass;
            newchordogram{2,j} = chordtree{2,i};
            newchordogram{3,j} = num2bass(bass);
            break;
        end
    end
    if ismatchout == 0
        newchordogram{1,j} = bass;
        newchordogram{2,j} = 'n';
        newchordogram{3,j} = num2bass(bass);
    end
end

% gestalt newchordogram (merge chord types as human would do)
% only keep maj, min, aug, dim, dom, minmaj, /3 and /5
% if 2, 5, n type, merge to the nearest maj or min type with same bass
%   if there is no nearest maj or min, merge n to 2 or 5, if both 2 and 5
%   exist in the neigbourhood, merge them to be sus2
% if /3 and /5, merge from nearest whatever type with same bass
% if sus2, sus4, merge to maj with same bass
newchordogram = [{0;'N';'N'} newchordogram {0;'N';'N'}];
yes = 1;
while yes
    yes = 0;
    for i = 2:1:length(newchordogram)-1
        cb = newchordogram{1,i};
        ct = newchordogram{2,i};
        pcb = newchordogram{1,i-1};
        pct = newchordogram{2,i-1};
        ncb = newchordogram{1,i+1};
        nct = newchordogram{2,i+1};
        if strcmp(ct,'n')
            if pcb == cb && ~(strcmp(pct,'n'))&& ~(strcmp(pct,'maj/3'))&& ~(strcmp(pct,'maj/5'))
                ct = pct;
                newchordogram{2,i} = ct;
                yes = 1;
                continue;
            end
            if ncb == cb && ~(strcmp(nct,'n'))&& ~(strcmp(nct,'maj/3'))&& ~(strcmp(nct,'maj/5'))
                ct = nct;
                newchordogram{2,i} = ct;
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'2')
            if pcb == cb && strcmp(pct,'5')
                ct = 'sus2';
                newchordogram{2,i} = ct;
                yes = 1;
                continue;
            end
            if ncb == cb && strcmp(nct,'5')
                ct = 'sus2';
                newchordogram{2,i} = ct;
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'5')
            if pcb == cb && strcmp(pct,'2')
                ct = 'sus2';
                newchordogram{2,i} = ct;
                yes = 1;
                continue;
            end
            if ncb == cb && strcmp(nct,'2')
                ct = 'sus2';
                newchordogram{2,i} = ct;
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'5') || strcmp(ct,'2')|| strcmp(ct,'n')
            if pcb == cb && ~(strcmp(pct,'5') || strcmp(pct,'2')|| strcmp(pct,'n')||strcmp(pct,'maj/3')||strcmp(pct,'maj/5'))
                ct = pct;
                newchordogram{2,i} = ct;
                yes = 1;
                continue;
            end
            if ncb == cb && ~(strcmp(nct,'5') || strcmp(nct,'2')|| strcmp(nct,'n')||strcmp(nct,'maj/3')||strcmp(nct,'maj/5'))
                ct = nct;
                newchordogram{2,i} = ct;
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'maj/3')
            if mod(pcb-4-1,12)+1 == cb && ~strcmp(pct,'maj/3') && ~strcmp(pct,'maj/5')
                newchordogram{1,i-1} = cb;
                newchordogram{2,i-1} = 'maj/3';
                newchordogram{3,i-1} = num2bass(cb);
                yes = 1;
                continue;
            end
            if mod(ncb-4-1,12)+1 == cb && ~strcmp(nct,'maj/3') && ~strcmp(nct,'maj/5')
                newchordogram{1,i+1} = cb;
                newchordogram{2,i+1} = 'maj/3';
                newchordogram{3,i+1} = num2bass(cb);
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'maj/5')
            if mod(pcb-7-1,12)+1 == cb && ~strcmp(pct,'maj/3') && ~strcmp(pct,'maj/5')
                newchordogram{1,i-1} = cb;
                newchordogram{2,i-1} = 'maj/5';
                newchordogram{3,i-1} = num2bass(cb);
                yes = 1;
                continue;
            end
            if mod(ncb-7-1,12)+1 == cb && ~strcmp(nct,'maj/3') && ~strcmp(nct,'maj/5')
                newchordogram{1,i+1} = cb;
                newchordogram{2,i+1} = 'maj/5';
                newchordogram{3,i+1} = num2bass(cb);
                yes = 1;
                continue;
            end
        end
        if strcmp(ct,'sus2') || strcmp(ct,'sus4')
            newchordogram{2,i} = 'maj';
            yes = 1;
            continue;
        end
    end
end
newchordogram = newchordogram(:,2:end-1);

figure;
hold on;
Y = -10:0.1:10;
for i = 1:1:lenOut+1
    X = chordboundaries(1,i)*ones(size(Y));
    plot(X,Y);
end
hold off;
div = (max(Y) - min(Y) - 1) / lenOut;
for i = 1:1:lenOut
%     x = (chordboundaries(i)+ chordboundaries(i+1)) / 2;
    x = chordboundaries(1,i);
    text(x,10 - i*div,strcat(newchordogram{3,i},newchordogram{2,i}));
end
xlabel('time');
ylabel('chord');
title('updated chordprogression vs. time');

% ********************* End of System A ******************** %
display('end of system A...');
% CPRSA - Chord Progression Recognition System A
% This implements a machine listening method for recognizing the chords
% and chord progression boundaries within a piece of music

% ********************************************************** %
% ********************* Input ****************************** %
% ********************************************************** %
close all;
clear;
clc;

display('input stage -- read audio from path');
% input stage
root = '../AudioSamples/';
audio = 'tuihou/tuihou-002.mp3';
path = [root audio];
[x, fs, songinfo, DSR] = myInput(path);

% ********************************************************** %
% ********************* Front End ************************** %
% ********************************************************** %
display('frontend -- from input to harmonic salience matrix');

% transform time domain into frequency domain
wl = 4096;
hopsize = 512;
X = mySpectrogram(x, wl, hopsize);
tt = (1/fs)*(1:length(x));
kk = (1:length(x));
ff = fs/2*linspace(0,1,wl/2);
myImagePlot(X, kk, ff, 'slice', 'Hz', 'spectrogram');

fmin = 27.5; % MIDI note 21
fmax = 1661; % MIDI note 92
fratio = 2^(1/36);
numtones = 72*3;
[Ms,Mc] = toneProfileGen(wl, numtones, fmin, fmax, fratio, fs);
myImagePlot(Ms, wl/2, numtones, 'time', '1/3 semitone', 'simple tone profile');
myImagePlot(Mc, wl/2, numtones, 'time', '1/3 semitone', 'complex tone profile');

% calculate note salience matrix of the stft spectrogram (cosine
% similarity) (note that this is an additive approach, as contrast to
% the nnls approach which is an deductive approach)
Ss = Ms*X;
Sc = Mc*X;
Spre = Sc.*Sc;
sizeM = size(Ms);
sizeX = size(X);
ps = 1:sizeM(1);
myImagePlot(Ss, kk, ps, 'slice', '1/3 semitone', 'simple tone salience matrix');
myImagePlot(Sc, kk, ps, 'slice', '1/3 semitone', 'complex tone salience matrix');
myImagePlot(Ss, kk, ps, 'slice', '1/3 semitone', 'preliminary salience matrix');

% noise reduction process
sizeNS = size(Spre);
nt = 0.1;
Ssn = noteSalienceNoiseReduce(Ss, nt);
Scn = noteSalienceNoiseReduce(Sc, nt);
myImagePlot(Ssn, kk, ps, 'slice', '1/3 semitone', 'noised reduced simple tone salience matrix');
myImagePlot(Scn, kk, ps, 'slice', '1/3 semitone', 'noised reduced complex tone salience matrix');

% computer note salience matrix by combining 1/3 semitones into semitones
Spres = Ssn(1:3:end,:) + Ssn(2:3:end,:) + Ssn(3:3:end,:);
Sprec = Scn(1:3:end,:) + Scn(2:3:end,:) + Scn(3:3:end,:);
S = Spres.*Sprec;
sizeS = size(S);
ntones = sizeS(1);
nslices = sizeS(2);
for j = 1:1:nslices
    if max(S(:,j)) > 0
        S(:,j) = S(:,j) / max(S(:,j)); 
    end
end
p = 1:ntones;
myImagePlot(S, kk, p, 'slice', 'semitone', 'note salience matrix');

% gestaltize the note salience matrix
wgmax = 20;
wpg = 0;
wng = 10;
Sg = gestaltNoteSalience(S,wgmax, wpg, wng);
myImagePlot(Sg, kk, p, 'slice', 'semitone', 'gestalt note salience matrix');

% onset filter (roughly detect the note onsets)
ot = 0.0;
wo = 10;
So = noteOnsetFilter(Sg, ot, wo);
myImagePlot(So, kk, p, 'slice', 'semitone', 'note onset matrix');

% bassline filter (roughly set the dynamic bass bounds)
bt = 0.0;
wb = 10;
cb = 1;
Sb = bassLineFilter(Sg, bt, wb, cb);
myLinePlot(1:length(Sb), Sb, 'slice', 'semitone', nslices, ntones, '*', 'rough bassline');

% harmonic change filter (detect harmonic change boundaries)ht = 0.1;
ht = 0.1;
[Sh, Shv, Shc, nchords] = harmonicChangeFilter(Sg, Sb, So, ht);
myImagePlot(Sh, kk, p, 'slice', 'semitone', 'harmonic bounded salience matrix');
myImagePlot(Shv, kk(1:nchords), p, 'chord progression order', 'semitone', 'harmonic change matrix');
myLinePlot(1:length(Shc), Shc, 'chord progression order', 'slice', nchords, nslices, 'o', 'haromonic change moments');

% ********************************************************** %
% ********************* Mid End - A************************* %
% ********************************************************** %
display('midend-A -- uppergram and basegram');

% compute basegram and uppergram (based on harmonic change matrix)
basegram = computeBasegram(Shv);
uppergram = computeUppergram(Shv);

treblenotenames = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
bassnotenames = {'N','C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
kh = 1:nchords;
ph = 1:12;
myImagePlot(uppergram, kh, ph, 'chord progression order', 'semitone', 'uppergram', ph,treblenotenames);
myLinePlot(kh, basegram(1,:), 'chord progression order', 'semitone', nchords, 12, 'o', 'basegram', 0:12, bassnotenames);

% ********************************************************** %
% ********************* Back End - A************************ %
% ********************************************************** %
display('backend-A -- chordtree');
chordtree = buildChordTree;
chordogram = computeChordogram(basegram, uppergram, chordtree);

chordogram = gestaltizeChordogram(chordogram, chordtree);

[outchordogram, outboundaries] = combineSameChords(chordogram, Shc);

visualizeChordProgression(outchordogram, outboundaries);

writeChordProgression(audio(1:end-4), nslices, hopsize, fs, outchordogram, outboundaries);

chordprogression = fullInfoChordProgression(outchordogram);

% ********************************************************** %
% ********************* Feedback Once - A******************* %
% ********************************************************** %
display('press enter to continue with feedback stage...');
pause;

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
    wb = chordboundaries(i):chordboundaries(i+1);
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
nchordtype = 28;
chordtree = cell(2,nchordtype);

chordtree{1,1} = [3,7];
chordtree{2,1} = 'min';

chordtree{1,2} = [4,7];
chordtree{2,2} = 'maj';

chordtree{1,3} = [4,8];
chordtree{2,3} = 'aug';

chordtree{1,4} = [3,6];
chordtree{2,4} = 'dim';

chordtree{1,5} = [3,10];
chordtree{2,5} = 'min7';

chordtree{1,6} = [4,11];
chordtree{2,6} = 'maj7';

chordtree{1,7} = [3,11];
chordtree{2,7} = 'minmaj7';

chordtree{1,8} = [4,10];
chordtree{2,8} = 'dom7';

chordtree{1,9} = [3,9];
chordtree{2,9} = 'min6';

chordtree{1,10} = [4,9];
chordtree{2,10} = 'maj6';

chordtree{1,11} = [5,7];
chordtree{2,11} = 'sus4';

chordtree{1,12} = [2,7];
chordtree{2,12} = 'sus2';

chordtree{1,13} = [3,8];
chordtree{2,13} = 'maj/3';

chordtree{1,14} = [5,9];
chordtree{2,14} = 'maj/5';

chordtree{1,15} = [2,5];
chordtree{2,15} = 'min/7';

chordtree{1,16} = [2,6];
chordtree{2,16} = 'maj/7';

chordtree{1,17} = [1,5];
chordtree{2,17} = 'maj/7+';

chordtree{1,18} = [2,10];
chordtree{2,18} = 'maj/2';

chordtree{1,19} = 3;
chordtree{2,19} = 'min';

chordtree{1,20} = 4;
chordtree{2,20} = 'maj';

chordtree{1,21} = 7;
chordtree{2,21} = '5';

chordtree{1,22} = 2;
chordtree{2,22} = '2';

chordtree{1,23} = 5;
chordtree{2,23} = '4';

chordtree{1,24} = 6;
chordtree{2,24} = '4#';

chordtree{1,25} = 8;
chordtree{2,25} = '6b';

chordtree{1,26} = 9;
chordtree{2,26} = '6';

chordtree{1,27} = 10;
chordtree{2,27} = '7b';

chordtree{1,28} = 11;
chordtree{2,28} = '7';

newchordogram = cell(4,lenOut);
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
            % convert slash chords to correct treble name, keep bass num
            upperbass = bass;
            if i == 13 % maj/3
                upperbass = mod(bass-4-1,12)+1;
            end
            if i == 14 % maj/5
                upperbass = mod(bass-7-1,12)+1;
            end
            if i == 15 % min/7
                upperbass = mod(bass-10-1,12)+1;
            end
            if i == 16 % maj/7
                upperbass = mod(bass-10-1,12)+1;
            end
            if i == 17 % maj/7+
                upperbass = mod(bass-11-1,12)+1;
            end
            if i == 18 % maj/2
                upperbass = mod(bass-2-1,12)+1;
            end
            newchordogram{1,j} = bass;
            newchordogram{2,j} = chordtree{2,i};
            newchordogram{3,j} = num2bass(upperbass);
            newchordogram{4,j} = i;
            break;
        end
    end
    if ismatchout == 0
        newchordogram{1,j} = bass;
        newchordogram{2,j} = '0';
        newchordogram{3,j} = num2bass(bass);
        newchordogram{4,j} = 0;
    end
end

% gestaltize newchordogram ---
% only keep triad, seventh and slash chords
% if diad type, merge to the nearest triad or seventh type with same bass
% if there is no nearest triad or seventh types, merge diads to try to form one
newchordogram = [{0;'N';'N';0} newchordogram {0;'N';'N';0}];
yes = 1;
while yes
    yes = 0;
    for i = 2:1:length(newchordogram)-1
        cb = newchordogram{1,i};
        ct = newchordogram{2,i};
        ctn = newchordogram{3,i};
        cte = newchordogram{4,i};
        pcb = newchordogram{1,i-1};
        pct = newchordogram{2,i-1};
        pctn = newchordogram{3,i-1};
        pcte = newchordogram{4,i-1};
        ncb = newchordogram{1,i+1};
        nct = newchordogram{2,i+1};
        nctn = newchordogram{3,i+1};
        ncte = newchordogram{4,i+1};
        % merge diad with nearest triad or seventh with same bass
        if (str2double(ct(1))) >= 0  % this indicates diad
            if pcb == cb && strcmp(ct,'0')
                newchordogram{2,i} = pct;
                newchordogram{3,i} = pctn;
                newchordogram{4,i} = pcte;
                yes = 1;
                continue;
            end
            
            if ncb == cb && strcmp(ct,'0')
                newchordogram{2,i} = nct;
                newchordogram{3,i} = nctn;
                newchordogram{4,i} = ncte;
                yes = 1;
                continue;
            end
            
            if pcb == cb && ~(str2double(pct(1)) >=0)
                newchordogram{2,i} = pct;
                newchordogram{3,i} = pctn;
                newchordogram{4,i} = pcte;
                yes = 1;
                continue;
            end
            
            if ncb == cb && ~(str2double(nct(1)) >=0)
                newchordogram{2,i} = nct;
                newchordogram{3,i} = nctn;
                newchordogram{4,i} = ncte;
                yes = 1;
                continue;
            end
            
            % if other diads are near by, try to merge them
            if pcb == cb && (str2double(pct(1)) >=0) && ~strcmp(pct,ct) && length(pct) <= 2 && length(ct) <= 2
                tri = sort([pcte cte]);
                ismatchout = 0;
                for ci = 1:1:length(chordtree)
                    ctri = chordtree{1,ci};
                    if length(ctri) > 1
                        if ctri(1) == tri(1) && ctri(2) == tri(2)
                            % there's a triad match
                            upperbass = bass;
                            if ci == 7 % maj/3
                                upperbass = mod(cb-4-1,12)+1;
                            end
                            if ci == 10 % maj/5
                                upperbass = mod(cb-7-1,12)+1;
                            end
                            if ci == 13 % min/7
                                upperbass = mod(cb-10-1,12)+1;
                            end
                            if ci == 14 % maj/7
                                upperbass = mod(cb-10-1,12)+1;
                            end
                            if ci == 15 % maj/7+
                                upperbass = mod(cb-11-1,12)+1;
                            end
                            if ci == 16 % maj/2
                                upperbass = mod(cb-2-1,12)+1;
                            end
                            newchordogram{2,i} = chordtree{2,ci};
                            newchordogram{3,i} = num2bass(upperbass);
                            newchordogram{4,i} = ci;
                            ismatchout = 1;
                            break;
                        end
                    end
                end
                if ismatchout == 0
                    % if there isn't a triad match, create a new label
                    newchordogram{2,i} = strcat(ct,pct);
                end
                yes = 1;
                continue;
            end
            
            % if other diads are near by, try to merge them
            if ncb == cb && (str2double(nct(1)) >=0)  && ~strcmp(nct,ct) && length(nct) <= 2 && length(ct) <= 2
                tri = sort([ncte cte]);
                ismatchout = 0;
                for ci = 1:1:length(chordtree)
                    ctri = chordtree{1,ci};
                    if length(ctri) > 1
                        if ctri(1) == tri(1) && ctri(2) == tri(2)
                            % there's a triad match
                            upperbass = bass;
                            if ci == 7 % maj/3
                                upperbass = mod(cb-4-1,12)+1;
                            end
                            if ci == 10 % maj/5
                                upperbass = mod(cb-7-1,12)+1;
                            end
                            if ci == 13 % min/7
                                upperbass = mod(cb-10-1,12)+1;
                            end
                            if ci == 14 % maj/7
                                upperbass = mod(cb-10-1,12)+1;
                            end
                            if ci == 15 % maj/7+
                                upperbass = mod(cb-11-1,12)+1;
                            end
                            if ci == 16 % maj/2
                                upperbass = mod(cb-2-1,12)+1;
                            end
                            newchordogram{2,i} = chordtree{2,ci};
                            newchordogram{3,i} = num2bass(upperbass);
                            newchordogram{4,i} = ci;
                            ismatchout = 1;
                            break;
                        end
                    end
                end
                if ismatchout == 0
                    % if there isn't a triad match, create a new label
                    newchordogram{2,i} = strcat(ct,nct);
                end
                yes = 1;
                continue;
            end
        end
    end
end
newchordogram = newchordogram(:,2:end-1);

% combine the same items
prevchord = strcat(newchordogram{3,1},newchordogram{2,1});
outchordogram = cell(1,lenOut);
outchordogram{1} = prevchord;
outboundaries = zeros(1,lenOut+1);
outboundaries(1) = chordboundaries(1);
outidx = 2;
for i = 2:1:lenOut
    curchord = strcat(newchordogram{3,i},newchordogram{2,i});
    if strcmp(curchord,prevchord)
        continue;
    else
        outchordogram{outidx} = curchord;
        outboundaries(outidx) = chordboundaries(i);
        outidx = outidx + 1;
    end
    prevchord = curchord;
end
outchordogram = outchordogram(1:outidx-1);
outboundaries(outidx) = chordboundaries(end);
outboundaries = outboundaries(1:outidx);

figure;
hold on;
Y = -10:0.1:10;
for i = 1:1:length(outboundaries)
    X = outboundaries(1,i)*ones(size(Y));
    plot(X,Y);
end
hold off;
div = (max(Y) - min(Y) - 1) / length(outchordogram);
for i = 1:1:length(outchordogram)
    x = outboundaries(i);
    text(x,10 - i*div,outchordogram{i});
end
xlabel('time');
ylabel('chord');
title('updated chordprogression vs. slice');

audiopath = audio(1:end-4);
if ~isempty(strfind(audio,'/'))
    tmp = strsplit(audiopath,'/');
    filename = tmp{2};
else
    filename = audiopath;
end
fw = fopen([filename '.ct.txt'],'w');
formatSpec1 = '%s';
formatSpec2 = '%s\n';
T = nslices; % the total number of slices contained in the evidence
tw = ((hopsize/fs)*(1:T));
lenoutchordogram = length(outchordogram);
for i = 1:1:lenoutchordogram
    if i == 1
        s = [outchordogram{i} '===>' num2str(0)];
    else
        s = [outchordogram{i} '===>' num2str(tw(outboundaries(i)))];
    end
    fprintf(fw, formatSpec1, s);
    s = ['-' num2str(tw(outboundaries(i+1)))];
    if i < lenoutchordogram
        fprintf(fw, formatSpec2, s);
    else
        fprintf(fw, formatSpec1, s);
    end
end
fclose(fw);

% ********************* End of System A ******************** %
display(strcat('end of system A recognizing...',filename));
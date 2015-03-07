% CPRSA - Chord Progression Recognition System A
% This implements a machine listening method for recognizing the chords
% and chord progression boundaries within a piece of music

% ********************************************************** %
% ********************* Input ****************************** %
% ********************************************************** %
close all;
clear;
clc;

feedbackpause = 0;

display('input stage -- read audio from path');
% input stage
root = '../AudioSamples/';
audio = 'haoting/haoting.02.mp3';
path = [root audio];
[x, fs] = myInput(path);


% ********************************************************** %
% ********************* Front End ************************** %
% ********************************************************** %
display('frontend -- from input to harmonic salience matrix');

% transform time domain into frequency domain
wl = 4096;
hopsize = 512;
X = mySpectrogram(x, wl, hopsize);
sizeX = size(X);
nslices = sizeX(2);
tk = (1/fs)*(1:length(x));
kk = (1:nslices);
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
S = normalizeGram(S);
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
myLinePlot(1:length(Sb), Sb, 'slice', 'semitone',...
    nslices, ntones, '*', 'rough bassline');

% harmonic change filter (detect harmonic change boundaries)ht = 0.1;
ht = 0.1;
bc = 6;
[Sh, Shv, Shc, nchords] = harmonicChangeFilter(Sg, Sb, So, ht, bc);
myImagePlot(Sh, kk, p, 'slice', 'semitone', 'harmonic bounded salience matrix');
myImagePlot(Shv, kk(1:nchords), p, 'chord progression order',...
    'semitone', 'harmonic change matrix');
myLinePlot(1:length(Shc), Shc, 'chord progression order', 'slice',...
    nchords, nslices, 'o', 'haromonic change moments');

% ********************************************************** %
% ********************* Mid End - A************************* %
% ********************************************************** %
display('midend-A -- uppergram and basegram');

% compute basegram and uppergram (based on harmonic change matrix)
basegram = computeBasegram(Shv);
uppergram = computeUppergram(Shv);
bassnotenames = {'N','C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
treblenotenames = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
kh = 1:nchords;
ph = 1:12;
myLinePlot(kh, basegram(1,:), 'chord progression order', 'semitone',...
    nchords, 12, 'o', 'basegram', 0:12, bassnotenames);
myImagePlot(uppergram, kh, ph, 'chord progression order', 'semitone',...
    'uppergram', ph,treblenotenames);

% ********************************************************** %
% ********************* Back End - A************************ %
% ********************************************************** %
display('backend-A -- chordmode');

chordmode = buildChordMode;

chordogram = computeChordogram(basegram, uppergram, chordmode);

[outchordogram, outbassgram, outboundaries] = combineSameChords(chordogram, Shc);

[outchordogram, outbassgram, outboundaries] = mergeDifferentChords(outchordogram, outbassgram, outboundaries); 

myLinePlot(1:length(outbassgram), outbassgram, 'chord progression order', 'semitone',...
    length(outbassgram), 12, 'o', 'outbassgram', 0:12, bassnotenames);
myLinePlot(1:length(outboundaries), outboundaries, 'chord progression order', 'slice',...
    length(outboundaries), nslices, 'o', 'outboundaries');

visualizeChordProgression(outchordogram, outboundaries);

writeChordProgression(path, nslices, hopsize, fs, outchordogram, outboundaries);

% ********************************************************** %
% ********************* Feedback Once - A******************* %
% ********************************************************** %
if feedbackpause
    display('press enter to continue with the feedback stage...');
    pause;
end
% ****** feedback mid end ****** %
display('feedback-A -- use chord boundaries information to do it again');

ut = 1;
[newbasegram, newuppergram] = updateBaseUpperGram(outbassgram, outboundaries, S, ut);
kh = 1:length(newbasegram);
ph = 1:12;
myLinePlot(kh, newbasegram(1,:), 'chord progression order', 'semitone', nchords, 12, 'o', 'basegram', 0:12, bassnotenames);
myImagePlot(newuppergram, kh, ph, 'chord progression order', 'semitone', 'uppergram', ph,treblenotenames);

% ****** feedback back end ****** %
newchordogram = computeChordogram(newbasegram, newuppergram, chordmode);

[newoutchordogram, newoutbassgram, newoutboundaries] = combineSameChords(newchordogram, outboundaries);

[newoutchordogram, newoutbassgram, newoutboundaries] = mergeDifferentChords(newoutchordogram, newoutbassgram, newoutboundaries); 

myLinePlot(1:length(newoutbassgram), newoutbassgram, 'chord progression order', 'semitone',...
    length(newoutbassgram), 12, 'o', 'newoutbassgram', 0:12, bassnotenames);
myLinePlot(1:length(newoutboundaries), newoutboundaries, 'chord progression order', 'slice',...
    length(newoutboundaries), nslices, 'o', 'newoutboundaries');

visualizeChordProgression(newoutchordogram, newoutboundaries);

writeChordProgression(path, nslices, hopsize, fs, newoutchordogram, newoutboundaries);

% ********************* End of System A ******************** %
display(strcat('end of system A recognizing:',path));
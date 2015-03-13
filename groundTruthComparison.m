% compare the output of this algorithm with ground truth.
% the comparison goes with the following settings:
% bass relative correct overlap (brco) = 
%   # bass matching frames / total # of frames
% treble relative correct overlap (trco) = 
%   # treble matching frames / total # of frames
% chord relative correct overlap (rco) =
%   # chord matching frames / total # of frames
% segmentation quality (sq) =
%   (# of correct boundaries - # of incorrect boundaries) / total # of boundaries
%   (or alternatively, use "Hamming divergence")
% chord progression quality (cpq) = 
%   sum of # of frames of 3 longest matching sequence / total # of frames
%
% definition of correct bass matching:
% definition of correct treble matching:
% definition of correct chord matching:
%
% definition of significant difference between algorithms:
%   Friedman rank
function rco = groundTruthComparison(target)
gtpath = strcat('gt/',target,'.gt.lrc'); % groundtruth output path
cdpath = strcat('cd/',target,'.cd.lrc'); % chordino output path
cppath = strcat('cp/',target,'.cp.lrc'); % cprs output path

bassnotenames = {'N','C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
% treble type 1 = major, type 2 = minor, type 3 = universal

% read groundtruth output into bassgram, treblegram and boundaries
gtbassgram = [];
gttreblegram = [];
gtboundaries = [];
fr = fopen(gtpath,'r');
tline = fgets(fr);
while ischar(tline) && length(tline) > 1
    disp(tline);
    [st,ch] = strtok(tline,']');
    st = st(2:end); % start time
    ch = ch(2:end); % chord
    % parse bass and treble
    bass = 0;
    treble = 'N';
    if (length(ch) == 1 && ~strcmp(ch,'N'))
        bass = bass2num(ch);
        treble = '';
    end
    if (length(ch) > 1)
        if ch(2) ~= '#' || ch(2) ~= 'b'
            bass = bass2num(ch(1));
            treble = ch(2:end);
        else
            bass = bass2num(ch(1:2));
            treble = ch(3:end);
        end
    end
    gtbassgram = [gtbassgram bass];
    gttreblegram = [gttreblegram treble];
    
    % parse seconds
    [mm,ss] = strtok(st,':');
    ss = ss(2:end);
    sec = str2double(mm)*60 + str2double(ss);
    tline = fgets(fr);
end
fclose(fr);

% read chordino output into bassgram, treblegram and boundaries
fr = fopen(cdpath,'r');
tline = fgets(fr);
while ischar(tline) && length(tline) > 1
    disp(tline);
    tline = fgets(fr);
end
fclose(fr);

% read cprs output into bassgram, treblegram and boundaries
fr = fopen(cppath,'r');
tline = fgets(fr);
while ischar(tline) && length(tline) > 1
    disp(tline);
    tline = fgets(fr);
end
fclose(fr);
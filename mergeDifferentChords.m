% this is another gestalize process that merge chords with very short
% duration into the chords nearby with reasonable duration
function [chordogram, bassgram, treblegram, boundaries] = mergeDifferentChords(chordogram, bassgram, treblegram, boundaries, grain)

% chord beat duration calculation
nchords = length(chordogram);
difbdrys = boundaries(2:end) - boundaries(1:end-1);
chordbeats = round(difbdrys ./ (median(difbdrys) / 4));

% merge candidates selection:
% 1. findout short chord indexes
shortchords = chordbeats;
shortchords(shortchords <= grain) = 0;
shortchordidx = find(shortchords == 0);

% 2. findout non chord indexes that's with a relatively short duration
nonchords = find(treblegram == 0);
nonchordidx = [];
for i = 1:1:length(nonchords)
    nci = nonchords(i);
    if chordbeats(nci) < 2
        nonchordidx = [nonchordidx nci];
    end
end

mergeecdgidxes = [shortchordidx nonchordidx];
if ~isempty(mergeecdgidxes)
    if mergeecdgidxes(end) == nchords
        % the last chord cannot be merge
        mergeecdgidxes(end) = [];
    end
end
if ~isempty(mergeecdgidxes)
    if mergeecdgidxes(1) == 1
        % the first chord cannot be merge
        mergeecdgidxes(1) = [];
    end
end

% merge process
chordogram(mergeecdgidxes) = [];
bassgram(mergeecdgidxes) = [];
treblegram(mergeecdgidxes) = [];
boundaries(mergeecdgidxes) = [];


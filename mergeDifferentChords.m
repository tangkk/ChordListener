% this is another gestalize process that merge chords with very short
% duration into the chords nearby with reasonable duration
function [chordogram, bassgram, boundaries] = mergeDifferentChords(chordogram, bassgram, boundaries)

% chord beat duration calculation
nchords = length(chordogram);
difbdrys = boundaries(2:end) - boundaries(1:end-1);
chordbeats = round(difbdrys ./ (median(difbdrys) / 4));

% merge candidates selection
mergeeCandidates = chordbeats;
mergeeCandidates(mergeeCandidates ~= 0) = -1;
mergeecdgidxes = find(mergeeCandidates == 0);
if ~isempty(mergeecdgidxes)
    if mergeecdgidxes(end) == nchords
        % the last chord cannot be merge
        mergeecdgidxes(end) = [];
    end
    if mergeecdgidxes(1) == 1
        % the first chord cannot be merge
        mergeecdgidxes(1) = [];
    end
end

mergeebdryidxes = mergeecdgidxes + 1;

% merge process
chordogram(mergeecdgidxes) = [];
bassgram(mergeecdgidxes) = [];
boundaries(mergeebdryidxes) = [];


% this is another gestalize process that merge chords with very short
% duration into the chords nearby with reasonable duration
function [outchordogram, outbassgram, outboundaries] = mergeDifferentChords(outchordogram, outbassgram, outboundaries)

% chord beat duration calculation
nchords = length(outchordogram);
difbdrys = outboundaries(2:end) - outboundaries(1:end-1);
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
outchordogram(mergeecdgidxes) = [];
outbassgram(mergeecdgidxes) = [];
outboundaries(mergeebdryidxes) = [];


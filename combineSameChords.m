% combine the same chords
function [outchordogram, outbassgram, chordboundaries] = combineSameChords(chordogram, Shc)

prevchord = strcat(chordogram{3,1},chordogram{2,1});
outchordogram = cell(1,length(chordogram));
outbassgram = zeros(1,length(chordogram));
outchordogram{1} = prevchord;
outbassgram(1) = chordogram{1,1};
chordboundaries = zeros(1,length(chordogram)+1);
chordboundaries(1) = Shc(1);
outidx = 2;
for i = 2:1:length(chordogram)
    curchord = strcat(chordogram{3,i},chordogram{2,i});
    curbass = chordogram{1,i};
    if strcmp(curchord,prevchord)
        continue;
    else
        outchordogram{outidx} = curchord;
        outbassgram(outidx) = curbass;
        chordboundaries(outidx) = Shc(i);
        outidx = outidx + 1;
    end
    prevchord = curchord;
end
outchordogram = outchordogram(1:outidx-1);
outbassgram = outbassgram(1:outidx-1);
chordboundaries(outidx) = Shc(end);
chordboundaries = chordboundaries(1:outidx);
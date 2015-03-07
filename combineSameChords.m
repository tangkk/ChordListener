% combine the same chords
function [outchordogram, outbassgram, outtreblegram, chordboundaries] = combineSameChords(inchordogram, Shc)

prevchord = strcat(inchordogram{3,1},inchordogram{2,1});
outchordogram = cell(1,length(inchordogram));
outbassgram = zeros(1,length(inchordogram));
outtreblegram = zeros(1,length(inchordogram));
outchordogram{1} = prevchord;
outbassgram(1) = inchordogram{1,1};
outtreblegram(1) = inchordogram{4,1};
chordboundaries = zeros(1,length(inchordogram)+1);
chordboundaries(1) = Shc(1);
outidx = 2;
for i = 2:1:length(inchordogram)
    curchord = strcat(inchordogram{3,i},inchordogram{2,i});
    curbass = inchordogram{1,i};
    curtreble = inchordogram{4,i};
    if strcmp(curchord,prevchord)
        continue;
    else
        outchordogram{outidx} = curchord;
        outbassgram(outidx) = curbass;
        outtreblegram(outidx) = curtreble;
        chordboundaries(outidx) = Shc(i);
        outidx = outidx + 1;
    end
    prevchord = curchord;
end
outchordogram = outchordogram(1:outidx-1);
outbassgram = outbassgram(1:outidx-1);
outtreblegram = outtreblegram(1:outidx-1);
chordboundaries(outidx) = Shc(end);
chordboundaries = chordboundaries(1:outidx);
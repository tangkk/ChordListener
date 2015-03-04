function [treble, ctidx] = chordTreeEntryMatching(tri, chordtree)
treble = strcat(chordtree{2,tri(1)}, chordtree{2,tri(2)});
ctidx = 0;
nchordtype = length(chordtree);
for ci = 1:1:nchordtype
    ctri = chordtree{1,ci};
    if length(ctri) > 1
        if ctri(1) == tri(1) && ctri(2) == tri(2)
            % there's a triad match
            treble = chordtree{2,ci};
            ctidx = ci;
            break;
        end
    end
end
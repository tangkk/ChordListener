% try chord recognition tree method (or chord dictionary method)
% run the bass through this tree, if hit, then that's the chord rooted
% on the bass
% if miss, then run every other pitch through this tree, pick a hit with
% root closest to the bass to form a slash chord

% walk the tree using basegram and uppergram
% chordogram: 1st row: bass num; 2st row: treble type;
% 3st row: treble root; 4th row: chord entry index
function chordogram = computeChordogram(basegram, uppergram, chordtree)

nchordtype = length(chordtree);
nchords = length(basegram);
chordogram = cell(4,nchords);
for j = 1:1:nchords
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
            chordogram{1,j} = bass;
            chordogram{2,j} = chordtree{2,i};
            chordogram{3,j} = num2bass(upperbass);
            chordogram{4,j} = i;
            break;
        end
    end
    if ismatchout == 0
        chordogram{1,j} = bass;
        chordogram{2,j} = '0';
        chordogram{3,j} = num2bass(bass);
        chordogram{4,j} = 0;
    end
end
function [treble, ctidx] = trebleMatching(bass,upper,template, templatetype)

if sum(upper) >= 5 % this is the 'no-chord' condition
    treble = '0';
    ctidx = 0;
    return;
end

% match by chord tree priority
if strcmp(templatetype,'chordtree')
    chordtree = template;
    treble = '0';
    ctidx = 0;
    nchordtype = length(chordtree);
    for i = 1:1:nchordtype
        ismatchin = 1;
        chordentry = chordtree{1,i};
        lenentry = length(chordentry);
        for jj = 1:1:lenentry
            matchpos = pitchTranspose(bass, chordentry(jj));
            if upper(matchpos) == 0
                ismatchin = 0;
            end
        end
        if ismatchin == 1
            treble = chordtree{2,i};
            ctidx = i;
            break;
        end
    end
end

% match by final score in each chord modes, if ties, selects the one
% with lower index
if strcmp(templatetype, 'chordmode')
    chordmode = template;
    nchordtype = length(chordmode);
    treblescore = zeros(1,nchordtype);
    for i = 1:1:nchordtype
        uppermode = zeros(12,1);
        uppermode(pitchTranspose(bass, chordmode{1,i})) = chordmode{3,i};
        treblescore(i) = upper'*uppermode;
    end
    [maxscore,ctidx] = max(treblescore);
    if maxscore > 0
        treble = chordmode{2,ctidx};
    else
        treble = '0';
        ctidx = 0;
    end
end

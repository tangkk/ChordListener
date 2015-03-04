function [treble, ctidx] = trebleMatching(bass,upper,template, templatetype)

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

if strcmp(templatetype, 'chordmode')
end
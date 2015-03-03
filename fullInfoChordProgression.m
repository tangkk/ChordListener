function chordprogression = fullInfoChordProgression(outchordogram)

lenOut = length(outchordogram);
chordprogression = cell(4,lenOut);

for i = 1:1:lenOut
    chordname = outchordogram{i};
    if strcmp(chordname(2),'#')
        bassname = chordname(1:2);
        bassnum = bass2num(bassname);
        treblename = chordname(3:end);
    else
        bassname = chordname(1);
        bassnum = bass2num(bassname);
        treblename = chordname(2:end);
    end
    % transform slash chord back to original bassname and bassnum
    if length(treblename) > 1
        if strcmp(treblename(end-1:end),'/3')
            bassnum = mod(bassnum+4-1,12) + 1;
            bassname = num2bass(bassnum);
        end
        if strcmp(treblename(end-1:end),'/5')
            bassnum = mod(bassnum+7-1,12) + 1;
            bassname = num2bass(bassnum);
        end
        if strcmp(treblename(end-1:end),'/7')
            bassnum = mod(bassnum+10-1,12) + 1;
            bassname = num2bass(bassnum);
        end
        if strcmp(treblename(end-1:end),'/7+')
            bassnum = mod(bassnum+11-1,12) + 1;
            bassname = num2bass(bassnum);
        end
        if strcmp(treblename(end-1:end),'/2')
            bassnum = mod(bassnum+2-1,12) + 1;
            bassname = num2bass(bassnum);
        end
    end
    chordprogression{1,i} = chordname;
    chordprogression{2,i} = bassname;
    chordprogression{3,i} = treblename;
    chordprogression{4,i} = bassnum;
end
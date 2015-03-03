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
    bassnum = upperbass2bass(bassnum, treblename);
    bassname = num2bass(bassnum);
    
    chordprogression{1,i} = chordname;
    chordprogression{2,i} = bassname;
    chordprogression{3,i} = treblename;
    chordprogression{4,i} = bassnum;
end
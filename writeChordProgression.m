% write out the chord progression as lrc file and if possible, combine
% the chord lrc file with the word lrc file
function writeChordProgression(audiopath, nslices, hopsize, fs, outchordogram, outbassgram, outboundaries)

chordlrc = [audiopath(1:end-4) '.cp.lrc'];
fw = fopen(chordlrc,'w');
formatSpec1 = '%s';
formatSpec2 = '%s\n';
tw = ((hopsize/fs)*(1:nslices));
lenoutchordogram = length(outchordogram);
for i = 1:1:lenoutchordogram
    if i == 1
        sec = 0;
    else
        sec = tw(outboundaries(i));
    end
    timestr = strcat(num2str(floor(sec/60)),':',num2str(mod(sec,60)));
    if isempty(strfind(outchordogram{i}, '/'))
        chordstr = outchordogram{i};
    else
        originchordstr = strsplit(outchordogram{i},'/');
        treble = originchordstr{1};
        bassstr = num2bass(outbassgram(i));
        chordstr = strcat(treble, '/', bassstr);
    end
    s = strcat('[',timestr,']',chordstr);
    fprintf(fw, formatSpec2, s);
end
sec = tw(outboundaries(end));
timestr = strcat(num2str(floor(sec/60)),':',num2str(mod(sec,60)));
s = strcat('[',timestr,']');
fprintf(fw, formatSpec1, s);
fclose(fw);

wordlrc = [audiopath(1:end-4) '.word.lrc'];
if exist(wordlrc, 'file')
    fword = fopen(wordlrc,'r');
    fsheet = fopen(chordlrc,'r');
    binword = fread(fword);
    binchord = fread(fsheet);
    binsheet = [binword ; binchord];
    fclose(fsheet);
    fclose(fword);
    sheetlrc = [audiopath(1:end-4) '.sheet.lrc'];
    fw = fopen(sheetlrc,'w');
    fwrite(fw, binsheet);
    fclose(fw);
end
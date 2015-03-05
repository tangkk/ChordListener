% write out the chord progression as lrc file and if possible, combine
% the chord lrc file with the word lrc file
function writeChordProgression(audiopath, nslices, hopsize, fs, outchordogram, outboundaries)

chordlrc = [audiopath(1:end-4) '.chord.lrc'];
fw = fopen(chordlrc,'w');
formatSpec1 = '%s';
formatSpec2 = '%s\n';
tw = ((hopsize/fs)*(1:nslices));
lenoutchordogram = length(outchordogram);
for i = 1:1:lenoutchordogram
    sec = tw(outboundaries(i));
    timestr = strcat(num2str(floor(sec/60)),':',num2str(mod(sec,60)));
    s = strcat('[',timestr,']',outchordogram{i});
    fprintf(fw, formatSpec2, s);
end
sec = tw(outboundaries(end));
timestr = strcat(num2str(floor(sec/60)),':',num2str(mod(sec,60)));
s = strcat('[',timestr,']');
fprintf(fw, formatSpec1, s);
fclose(fw);

% wordlrc = [audiopath(1:end-4) '.word.lrc'];
% if exist(wordlrc, 'file')
%     textword = fileread(wordlrc);
%     textchord = fileread(chordlrc);
%     textsheet = [textword textchord];
%     sheetlrc = [audiopath(1:end-4) '.sheet.lrc'];
%     fw = fopen(sheetlrc,'w');
%     fprintf(fw, formatSpec1,textsheet);
%     fclose(fw);
% end
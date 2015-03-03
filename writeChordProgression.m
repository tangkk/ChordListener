function writeChordProgression(audiopath, nslices, hopsize, fs, outchordogram, outboundaries)

if ~isempty(strfind(audiopath,'/'))
    tmp = strsplit(audiopath,'/');
    filename = tmp{2};
else
    filename = audiopath;
end
fw = fopen([filename '.ct.txt'],'w');
formatSpec1 = '%s';
formatSpec2 = '%s\n';
T = nslices; % the total number of time slices contained in the evidence
tw = ((hopsize/fs)*(1:T));
lenoutchordogram = length(outchordogram);
for i = 1:1:lenoutchordogram
    if i == 1
        s = [outchordogram{i} '===>' num2str(0)];
    else
        s = [outchordogram{i} '===>' num2str(tw(outboundaries(i)))];
    end
    fprintf(fw, formatSpec1, s);
    s = ['-' num2str(tw(outboundaries(i+1)))];
    if i < lenoutchordogram
        fprintf(fw, formatSpec2, s);
    else
        fprintf(fw, formatSpec1, s);
    end
end
fclose(fw);
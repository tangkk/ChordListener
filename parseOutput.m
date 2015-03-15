function [bassgram, treblegram, treblenames, boundaries] = parseOutput(path)

bassgram = [];
treblegram = [];
treblenames = {};
boundaries = [];
fr = fopen(path,'r');
tline = fgets(fr);
while ischar(tline) && length(tline) > 1
    [st,ch] = strtok(tline,']');
    st = st(2:end); % start time
    ch = ch(2:end-1); % chord (delete the newline label)
    
    % parse bass and treble
    if ~isempty(ch)
        bass = 0;
        treble = 'N';
        if ~isempty(strfind(ch,'/'))
            [bass,treble] = slash2BassTreble(ch);
        elseif (length(ch) == 1 && ~strcmp(ch,'N'))
            bass = bass2num(ch);
            treble = '';
        elseif (length(ch) > 1)
            if ch(2) ~= '#' || ch(2) ~= 'b'
                bass = bass2num(ch(1));
                treble = ch(2:end);
            else
                bass = bass2num(ch(1:2));
                treble = ch(3:end);
            end
        end
        
        if strcmp(treble,'')
            treble = 'M';
        end
        % transform treble to built-in way
        treblenames = [treblenames treble];
        trebletype = trebleTypeMapping(treble);
        bassgram = [bassgram bass];
        treblegram = [treblegram trebletype];
    end
    
    % parse seconds
    [mm,ss] = strtok(st,':');
    ss = ss(2:end);
    sec = str2double(mm)*60 + str2double(ss);
    boundaries = [boundaries sec];
    tline = fgets(fr);
end
fclose(fr);

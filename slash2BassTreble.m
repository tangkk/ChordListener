function [bass,treble] = slash2BassTreble(ch)

if ~isempty(strfind(ch,'m/'))
    [ustr,bstr] = strtok(ch,'m/');
    bstr = bstr(3:end);
    upperbass = bass2num(ustr);
    bass = bass2num(bstr);

    trebleList = {'m/3','m/5','m/7','m/7+','m/2'};
    treble = 'N';
    for i = 1:1:length(trebleList)
        if upperbass == bass2upperbass(bass, trebleList{i})
            treble = trebleList{i};
            break;
        end
    end
else
    [ustr,bstr] = strtok(ch,'/');
    bstr = bstr(2:end);
    upperbass = bass2num(ustr);
    bass = bass2num(bstr);

    trebleList = {'/3','/5','/7','/7+','/2'};
    treble = 'N';
    for i = 1:1:length(trebleList)
        if upperbass == bass2upperbass(bass, trebleList{i})
            treble = trebleList{i};
            break;
        end
    end
end
% assume the files are all with '.mp3' suffix, this function sort the
% files into a correct chord progression order which is divided manually
% this is incase the number of files exceeeds 100
function files = sortFiles(path)
files = dir(path);
files = files(3:end);
fileidx = zeros(1, length(files));
for i = 1:1:length(files)
    name = files(i).name;
    name = name(1:end - 4);
    [~, idxstr] = strtok(name, '.');
    idxstr = idxstr(2:end);
    idx = str2double(idxstr);
    fileidx(i) = idx;
end
[fileidx, idx] = sort(fileidx);
files = files(idx);

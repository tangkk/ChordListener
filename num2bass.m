% 1 2  3 4  5 6 7  8 9  10 11 12 13
% C C# D D# E F F# G G# A  A# B  N
function bass = num2bass(num)

switch num
    case 0
        bass = 'N';
    case 1
        bass = 'C';
    case 2
        bass = 'C#';
    case 3
        bass = 'D';
    case 4
        bass = 'D#';
    case 5
        bass = 'E';
    case 6
        bass = 'F';
    case 7
        bass = 'F#';
    case 8
        bass = 'G';
    case 9
        bass = 'G#';
    case 10
        bass = 'A';
    case 11
        bass = 'A#';
    case 12
        bass = 'B';
end
% 1 2  3 4  5 6 7  8 9  10 11 12 13
% C C# D D# E F F# G G# A  A# B  N
function treble = num2bass(num)

switch num
    case 1
        treble = 'C';
    case 2
        treble = 'C#';
    case 3
        treble = 'D';
    case 4
        treble = 'D#';
    case 5
        treble = 'E';
    case 6
        treble = 'F';
    case 7
        treble = 'F#';
    case 8
        treble = 'G';
    case 9
        treble = 'G#';
    case 10
        treble = 'A';
    case 11
        treble = 'A#';
    case 12
        treble = 'B';
    case 13
        treble = 'N';
end
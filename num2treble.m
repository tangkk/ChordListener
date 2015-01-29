% 1 2  3 4  5 6 7  8 9  10 11 12 13 14  15 16  17 18 19  20 21  22 23  24
% C C# D D# E F F# G G# A  A# B  Cm C#m Dm D#m Em Fm F#m Gm G#m Am A#m Bm
% 25 26  27 28  29 30 31  32 33  34  35  36 37
% C5 C#5 D5 D#5 E5 F5 F#5 G5 G#5 A5  A#5 B5 N
function treble = num2treble(num)

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
        treble = 'Cm';
    case 14
        treble = 'C#m';
    case 15
        treble = 'Dm';
    case 16
        treble = 'D#m';
    case 17
        treble = 'Em';
    case 18
        treble = 'Fm';
    case 19
        treble = 'F#m';
    case 20
        treble = 'Gm';
    case 21
        treble = 'G#m';
    case 22
        treble = 'Am';
    case 23
        treble = 'A#m';
    case 24
        treble = 'Bm';
    case 25
        treble = 'C5';
    case 26
        treble = 'C#5';
    case 27
        treble = 'D5';
    case 28
        treble = 'D#5';
    case 29
        treble = 'E5';
    case 30
        treble = 'F5';
    case 31
        treble = 'F#5';
    case 32
        treble = 'G5';
    case 33
        treble = 'G#5';
    case 34
        treble = 'A5';
    case 35
        treble = 'A#5';
    case 36
        treble = 'B5';
    case 37
        treble = 'N';
end
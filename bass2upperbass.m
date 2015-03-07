function upperbass = bass2upperbass(bass, treble)

switch treble
    case '/3'
        transval = -4;
    case '/5'
        transval = -7;
    case 'm/7'
        transval = -10;
    case '/7'
        transval = -10;
    case '/7+'
        transval = -11;
    case '/2'
        transval = -2;
    otherwise
        transval = 0;
end

upperbass = pitchTranspose(bass, transval);
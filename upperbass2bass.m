function bass = upperbass2bass(upperbass, treble)

switch treble
    case '/3'
        transval = 4;
    case 'm/3'
        transval = 4;
    case '/5'
        transval = 7;
    case 'm/5'
        transval = 7;
    case '/7'
        transval = 10;
    case 'm/7'
        transval = 10;
    case '/7+'
        transval = 11;
    case 'm/7+'
        transval = 11;
    case '/2'
        transval = 2;
    case 'm/2'
        transval = 2;
    otherwise
        transval = 0;
end

bass = pitchTranspose(upperbass, transval);
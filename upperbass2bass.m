function bass = upperbass2bass(upperbass, treble)

switch treble
    case 'maj/3'
        transval = 4;
    case 'maj/5'
        transval = 7;
    case 'min/7'
        transval = 10;
    case 'maj/7'
        transval = 10;
    case 'maj/7+'
        transval = 11;
    case 'maj/2'
        transval = 2;
    otherwise
        transval = 0;
end

bass = pitchTranspose(upperbass, transval);
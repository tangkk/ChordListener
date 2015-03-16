% this maps treble to a treble class
% more trebles are to be added
function tt = trebleTypeMapping(treble)

switch treble
    % type 1 cases:
    case '0'
        tt = 1; % 0 type of treble also 1
    case ''
        tt = 1;
    case 'M'
        tt = 1;
    case 'maj'
        tt = 1;
    case 'maj6'
        tt = 1;
    case 'M6'
        tt = 1;
    case 'maj7'
        tt = 1;
    case 'M7'
        tt = 1;
    case '7'
        tt = 1;
    case 'dom7'
        tt = 1;
    case 'aug'
        tt = 1;
    case '+'
        tt = 1;
    case '+7'
        tt = 1;
        
    % type 2 cases
    case 'm'
        tt = 2;
    case 'min'
        tt = 2;
    case 'm6'
        tt = 2;
    case 'min6'
        tt = 2;
    case 'min7'
        tt = 2;
    case 'm7'
        tt = 2;
    case 'minmaj7'
        tt = 2;
    case 'mM7'
        tt = 2;
    case 'm7-5'
        tt = 2;
    case 'm7b5'
        tt = 2;
    case '-'
        tt = 2;
    case 'dim'
        tt = 2;
    case 'dim7'
        tt = 2;
    case '-7'
        tt = 2;
    case '/3'
        tt = 2;
        
    % type 3 cases
    case 'b2'
        tt = 3;
    case '2'
        tt = 3;
    case '4'
        tt = 3;
    case '4#'
        tt = 3;
    case '5'
        tt = 3;
    case '6b'
        tt = 3;
    case '6'
        tt = 3;
    case '7b'
        tt = 3;
    case '7+'
        tt = 3;
    case 'sus4'
        tt = 3;
    case 'sus2'
        tt = 3;
    case 'm/5'
        tt = 3;
    case 'm/7'
        tt = 3;
    case 'm/7+'
        tt = 3;
    case 'm/2'
        tt = 3;
    case '/5'
        tt = 3;
    case '/7'
        tt = 3;
    case '/7+'
        tt = 3;
    case '/2'
        tt = 3;

    otherwise
        tt = 0;
end
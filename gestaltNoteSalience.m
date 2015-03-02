% gestalize process
% if within a gestalt window ahead there's a non-zero bin, compensate the
% blank in between, the length of the gestalt window vary according to
% the accumulated non-blank length, but with maximum value of 20 slices
function Sg = gestaltNoteSalience(S,wgmax, wpg, wng)

sizeS = size(S);
Sg = zeros(sizeS(1), sizeS(2));
for i = 1:1:sizeS(1)
    trackidx = 1;
    isblank = 0;
    for j = 1:1:sizeS(2)
        if S(i,j) > 0
            % compensate the gestalt
            if isblank == 1
                lenBlank = j - trackidx;
                if lenBlank <= min(wgmax, wpg)
                    Sg(i,trackidx:j-1) = mean(S(i,max(trackidx-wpg,1):trackidx))*ones(1,lenBlank);
                    wpg = wpg + lenBlank;
                else
                    wpg = 0;
                end
            end
            Sg(i,j) = S(i,j);
            isblank = 0;
            wpg = wpg + 1;
            trackidx = j;
        else
            isblank = 1;
        end
    end
end

% input from above, if a piece of salience is shorter than a gestalt window, ignore it
for i = 1:1:sizeS(1)
    trackidx = 1;
    islight = 0;
    for j = 1:1:sizeS(2)
        if Sg(i,j) == 0
            if islight == 1;
                lenLight = j - trackidx;
                if lenLight <= wng
                    Sg(i,trackidx:j-1) = zeros(1,lenLight);
                end
            end
            trackidx = j;
            islight = 0;
        else
            islight = 1;
        end
    end
end

% gestaltize chordogram ---
% only keep triad, seventh and slash chords
% if diad type, merge to the nearest triad or seventh type with same bass
% if there is no nearest triad or seventh types, merge diads to try to form one
function chordogramout = gestaltizeChordogram(chordogramin, chordtree)

chordogramin = [{0;'N';'N';0} chordogramin {0;'N';'N';0}];
yes = 1;
while yes
    yes = 0;
    for i = 2:1:length(chordogramin)-1
        cb = chordogramin{1,i};
        ct = chordogramin{2,i};
        ctn = chordogramin{3,i};
        cte = chordogramin{4,i};
        pcb = chordogramin{1,i-1};
        pct = chordogramin{2,i-1};
        pctn = chordogramin{3,i-1};
        pcte = chordogramin{4,i-1};
        ncb = chordogramin{1,i+1};
        nct = chordogramin{2,i+1};
        nctn = chordogramin{3,i+1};
        ncte = chordogramin{4,i+1};
        % merge diad with nearest triad or seventh with same bass
        if (str2double(ct(1))) >= 0  % this indicates diad
            if pcb == cb && strcmp(ct,'0')
                chordogramin{2,i} = pct;
                chordogramin{3,i} = pctn;
                chordogramin{4,i} = pcte;
                yes = 1;
                continue;
            end
            
            if ncb == cb && strcmp(ct,'0')
                chordogramin{2,i} = nct;
                chordogramin{3,i} = nctn;
                chordogramin{4,i} = ncte;
                yes = 1;
                continue;
            end
            
            if pcb == cb && ~(str2double(pct(1)) >=0)
                chordogramin{2,i} = pct;
                chordogramin{3,i} = pctn;
                chordogramin{4,i} = pcte;
                yes = 1;
                continue;
            end
            
            if ncb == cb && ~(str2double(nct(1)) >=0)
                chordogramin{2,i} = nct;
                chordogramin{3,i} = nctn;
                chordogramin{4,i} = ncte;
                yes = 1;
                continue;
            end
            
            % if other diads are near by, try to merge them
            if pcb == cb && (str2double(pct(1)) >=0) && ~strcmp(pct,ct) && length(pct) <= 2 && length(ct) <= 2
                tri = sort([pcte cte]);
                ismatchout = 0;
                for ci = 1:1:length(chordtree)
                    ctri = chordtree{1,ci};
                    if length(ctri) > 1
                        if ctri(1) == tri(1) && ctri(2) == tri(2)
                            % there's a triad match
                            treble = chordtree{2,ci};
                            ucb = bass2upperbass(cb, treble);
                            chordogramin{2,i} = chordtree{2,ci};
                            chordogramin{3,i} = num2bass(ucb);
                            chordogramin{4,i} = ci;
                            ismatchout = 1;
                            break;
                        end
                    end
                end
                if ismatchout == 0
                    % if there isn't a triad match, create a new label
                    chordogramin{2,i} = strcat(ct,pct);
                end
                yes = 1;
                continue;
            end
            
            % if other diads are near by, try to merge them
            if ncb == cb && (str2double(nct(1)) >=0)  && ~strcmp(nct,ct) && length(nct) <= 2 && length(ct) <= 2
                tri = sort([ncte cte]);
                ismatchout = 0;
                for ci = 1:1:length(chordtree)
                    ctri = chordtree{1,ci};
                    if length(ctri) > 1
                        if ctri(1) == tri(1) && ctri(2) == tri(2)
                            % there's a triad match
                            treble = chordtree{2,ci};
                            ucb = bass2upperbass(cb, treble);
                            chordogramin{2,i} = chordtree{2,ci};
                            chordogramin{3,i} = num2bass(ucb);
                            chordogramin{4,i} = ci;
                            ismatchout = 1;
                            break;
                        end
                    end
                end
                if ismatchout == 0
                    % if there isn't a triad match, create a new label
                    chordogramin{2,i} = strcat(ct,nct);
                end
                yes = 1;
                continue;
            end
        end
    end
end
chordogramout = chordogramin(:,2:end-1);

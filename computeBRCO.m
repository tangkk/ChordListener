function brco = computeBRCO(gtbassgram, gtboundaries, bassgram, boundaries, TD)

BCD = 0; % bass correct duration
for i = 1:1:length(gtbassgram)
    st = gtboundaries(i);
    et = gtboundaries(i+1);
    boundst = locatebound(st,boundaries, 'st');
    bst = bassgram(boundst(1));
    boundet = locatebound(et,boundaries, 'et');
    bet = bassgram(boundet(1));
    boundset = [st boundaries(boundst(2):boundet(1)) et];
    bassset = [bst bassgram(boundst(2):boundet(1))];
    
    gtbass = gtbassgram(i);
    for j = 1:1:length(bassset)
        if bassset(j) == gtbass
            BCD = BCD + boundset(j+1) - boundset(j);
        end
    end
end
brco = BCD / TD;
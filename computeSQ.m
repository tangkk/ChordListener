function sq = computeSQ(gtboundaries, boundaries, BE)

SCN = 0;
TB = length(gtboundaries);
TDB = length(boundaries);
for i = 1:1:length(gtboundaries)
    st = gtboundaries(i);
    boundst = locatebound(st,boundaries, 'st');
    tbst = boundaries(boundst);
    if st - tbst(1) < BE || tbst(2) - st < BE % less than an error bound
        SCN = SCN + 1;
    end
end
sq = (SCN - (TDB - SCN)) / TB; % (correct bounds - incorrect bounds) / total bounds
% part of the feedback route, update the basegram and uppergram
function [newbasegram, newuppergram] = updateBaseUpperGram(chordprogression, chordboundaries, S, ut)

lenOut = length(chordprogression);
newbasegram = zeros(1,lenOut);
for i = 1:1:lenOut
    newbasegram(i) = chordprogression{4,i};
end
newuppergram = zeros(12,lenOut);

for i = 1:1:lenOut
    % update note salience matrix in terms of boundaries window
    wb = chordboundaries(i):chordboundaries(i+1);
    sm = sum(S(:,wb),2);
    sm(sm < ut) = 0;
    if max(sm) > 0
        sm = sm ./ max(sm);
    end
    upg = sum(reshape(sm,12,6),2);
    upg = [upg(4:end);upg(1:3)];
    newuppergram(:,i) = upg;
end
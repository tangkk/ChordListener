function uppergram = computeUppergram(Shv)

sizeShv = size(Shv);
ntones = sizeShv(1);
nchords = sizeShv(2);
uppergram = zeros(12,nchords);

for i = 1:1:12
    for k = 1:1:ntones/12
        uppergram(i,:) = uppergram(i,:) + Shv(i + 12*(k-1),:);
    end
end
% normalize uppergram
for j = 1:1:nchords
    if max(uppergram(:,j)) ~= 0
        uppergram(:,j) = uppergram(:,j) / max(uppergram(:,j));
        tmp = uppergram(:,j);
        % revise the order to be starting from C
        uppergram(:,j) = [tmp(4:end) ; tmp(1:3)];
    end
end
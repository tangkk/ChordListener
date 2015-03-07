function base = findBase(Sw, Swo)

ntones = size(Sw,1);
nslices = size(Sw,2);

basevec = zeros(1,ntones);
for j = 1:1:nslices
    for i = 1:1:ntones
        if Swo(i,j) > 0 % onset salience with priority
            basevec(i) = basevec(i) + Swo(i,j);
            break;
        end
        if Sw(i,j) > 0
            basevec(i) = basevec(i) + Sw(i,j);
            break;
        end
    end
end

[baseval,base] = max(basevec);
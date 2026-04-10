function idxs = getNonDominated(solutions)

% This function obtains the nondominates point is a set of solutions
%
% idxs = getNonDominates(solutions)
% 
% solutions is a m x n matrix, where m is the number of points (solutions)
% and n is the number of objectives.
%
% idxs contains the indexes for the non-dominated points. 

[ni,no] = size(solutions);
idxs    = zeros(size(solutions,1),1);
resTemp = zeros(size(solutions,1) - 1, 3);
counter = 0;

if (ni == 1)
    idxs = 1;
    return;
end

for j = 1:ni
    resTmp  = zeros(size(solutions));
    thisSolution = solutions(j,:);
    for k = 1:no
        resTmp(:,k) =  thisSolution(k) - solutions(:,k);
    end
    resTmp(j,:) = [];
    resTmp = sign(resTmp);
    resTemp(:,1) = sum(resTmp < 0,2);
    resTemp(:,2) = sum(resTmp == 0,2);
    resTemp(:,3) = sum(resTmp > 0,2);
    if (min(resTemp(:,1)) > 0)
        counter = counter + 1;
        idxs(counter) = j;
    else
        resTempz = resTemp(resTemp(:,1) == 0,:);
        idxz     = find(resTempz(:,2) == no);
        if (length(idxz) == size(resTempz,1))
            counter = counter + 1;
            idxs(counter) = j;
        end
    end
end
idxs(counter + 1:end) = [];
end
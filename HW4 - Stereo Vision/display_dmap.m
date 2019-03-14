function [dColor] = display_dmap(dMap)
    maxD = max(max(dMap));
    minD = min(min(dMap));
    dMap = dMap-minD;
    diff = maxD-minD;
    dMap = dMap/diff;
    dColor = repmat(dMap, 1, 1, 3);
    [rows, cols,~] = size(dColor);
    for r = 1:rows
        for c = 1:cols
            if isnan(dColor(r,c,:))
                dColor(r,c,1) = 1;dColor(r,c,2) = 0; dColor(r,c,3) = 0;
            end     
        end
    end

end
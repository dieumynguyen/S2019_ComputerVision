function [outlierMap] = check_consistency(LR, RL, TLR)
    [rows, cols] = size(LR);outlierMap = zeros(rows, cols);
    for row = 1:rows
        for col = 1:cols
            dlr = LR(row,col);x = col;
            if x + dlr < 1 || x + dlr > cols                
                outlierMap(row, col) = 1;
            else
                xdlr = round(x + dlr);drl =  RL(row, xdlr);
                if abs(dlr - drl) <= TLR
                    outlierMap(row,col) = 0;
                else
                    outlierMap(row,col) = 1;
                end
            end
        end
    end
end
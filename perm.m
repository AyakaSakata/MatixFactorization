function [permindx, parity, overlap] = perm( F0, F, N )
% return permutation index and parity
%
    permindx = zeros(N,1);
    parity = zeros(N,1);

    overlap = F'*F0;

    rowindx = 1:N;
    columnindx = 1:N;
    
    for n = 1: N
    
        tmp_matrix = abs(overlap(rowindx,columnindx));
    
        [MAX, I] = max(tmp_matrix);
        [MAX2, I2] = max(MAX);
    
        permindx(columnindx(I2)) = rowindx(I(I2));
        parity(columnindx(I2)) = sign(overlap(rowindx(I(I2)),columnindx(I2)));
        tmp_rowindx = rowindx;
        tmp_columnindx = columnindx;
    
        rowindx = zeros(1,N-n);
        columnindx = zeros(1,N-n);
    
        row_carryover = 0;
        column_carryover = 0;
    
        for k = 1: N-n+1
            if k ~= I2
                column_carryover = column_carryover+1;
                columnindx(column_carryover) = tmp_columnindx(k);
            end
            if k ~= I(I2)
                row_carryover = row_carryover+1;
                rowindx(row_carryover) = tmp_rowindx(k);
            end
        end

    end

end


function [ F ] = F_initial( Y , N )

    M = size(Y,1);
    P = size(Y,2);
    sqM = sqrt(M);
    
    ColumnNum = 0;
    DrawAttempts = 0;
    Ycond = cond(Y);
    
    while(ColumnNum <= N) && (DrawAttempts < 20)
            
        % Initialization
        atomOrder = randperm(P); % different random order each time
        F = zeros(M,N);
        F(:,1) = Y(:,atomOrder(1))/norm(Y(:,atomOrder(1)));
        ColumnNum = 2;
        counter = 2;
            
        % Try to draw the initial dictionary
        while ColumnNum <= N && counter <= P
                
            % Assign the new column
            F(:,ColumnNum) = Y(:,atomOrder(counter)) /...
                    norm(Y(:,atomOrder(counter)));
                
            % Check inner product and increment if different
            % also check cond to watch for rank deficient cases
            if (max(abs(F(:,1:(ColumnNum-1))'*F(:,ColumnNum))) < 0.9) && ...
                    (cond(F(:,1:ColumnNum)) < 10*Ycond)
                ColumnNum = ColumnNum + 1;
            end
                
            %Increment counter
            counter = counter + 1;
        end
            
        DrawAttempts = DrawAttempts + 1;
    end
    
    % Fill in with random if needed
    if ColumnNum <= N
        F(:,ColumnNum:end) = randn(size(F(:,ColumnNum:end)));
    end
        
    % Normalization
    F_norm = sqrt(diag(F'*F));
    F = F./(ones(M,1)*F_norm')*sqM;

end


function d = se3Distance(A1, A2)
    % Compute the geodesic distance between two transformations in SE(3)
    %
    % Inputs:
    %  
    % Outputs:
    %   d - the geodesic distance between (R1, t1) and (R2, t2)

    % Rotation part using the Frobenius norm of the logarithm of the
    % rotation difference 
    R1 = A1(1:3,1:3); 
    t1 = A1(1:3,4); 
    
    R2 = A2(1:3,1:3); 
    t2 = A2(1:3,4); 
    
    R = R1 * R2';
    angle = norm(logm(R), 'fro');
    
    % Translation part using Euclidean distance
    translation_dist = norm(t1 - t2);
    
    % Combine the distances
    d = sqrt(angle^2 + translation_dist^2);
end 
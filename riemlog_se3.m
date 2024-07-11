function L = riemlog_se3(A,B) % calcule le log Riemannien de A^{-1} * B 
    L = zeros(4,4); 
    L(1:3,1:3) = logm ( A(1:3,1:3)' * B(1:3,1:3) ); 
    L(1:3,4) = - A(1:3,4) + B(1:3,4); 
end 

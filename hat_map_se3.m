function A = hat_map_se3(v)
    A = zeros(4,4);
    A(1:3,1:3) = [0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0] ;
    A(1:3,4) = [v(4) ; v(5) ; v(6)] ;
end
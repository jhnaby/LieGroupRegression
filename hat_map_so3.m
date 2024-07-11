function A = hat_map_so3(v)
    A = zeros(3,3);
    A(1:3,1:3) = [0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0] ;
end
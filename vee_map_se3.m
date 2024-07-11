function t = vee_map_se3(A)
    
    
  %  if(double(sum(abs(diag(A)))) > 1e-2)
           % disp('Error in vee_map() - Nonzero diagonal') ;
           %  pause ;
  %  end
  
    
    t = [A(3,2) ;
         A(1,3) ;
         A(2,1) ;
         A(1,4) ;
         A(2,4) ;
         A(3,4)] ;
end
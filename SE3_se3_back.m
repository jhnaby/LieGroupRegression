function [ se3] = SE3_se3_back( SE3 )
    % The input is a 4-by-4 matrix, element of SE(3)
    % The output is a vetor of 6-by-1 coordinates,  the first 3
    % corresponding to the "angular velocity" omega, the last three being
    % the " translational velocity" v  
    R=SE3(1:3,1:3);
    theta=acos((trace(R)-1)/2);
    %lnR=(theta/(2*sin(theta)))*(R-R');
    %w=[-lnR(2,3) lnR(1,3) -lnR(1,2)];
    %wx=[0 -w(3) w(2);w(3) 0 -w(1);-w(2) w(1) 0];
    if(theta==0)
        lnR=zeros(3,3);
        Vin=eye(3);
        w=[0 0 0];
    else
        lnR=(theta/(2*sin(theta)))*(R-R');
        w=[-lnR(2,3) lnR(1,3) -lnR(1,2)];
        wx=[0 -w(3) w(2);w(3) 0 -w(1);-w(2) w(1) 0];
        A=sin(theta)/theta;
        B=(1-cos(theta))/(theta^2);
        Vin=eye(3)-(1/2)*wx+(1/(theta^2))*(1-(A/(2*B)))*(wx*wx);
        
    end
    u=Vin*SE3(1:3,4);
    se3=[ w u'];

end
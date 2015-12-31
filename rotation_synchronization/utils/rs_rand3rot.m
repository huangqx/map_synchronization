function rM = rs_rand3rot()
% rand3rot computes a random rotation matrix in 3D.
% Based on Graphics Gems 3 C code by James Arvo

if 0
    x1 = rand(1);
    x2 = rand(1);
    x3 = rand(1);

    R = [ cos(2*pi*x1)  sin(2*pi*x1) 0; -sin(2*pi*x1)  cos(2*pi*x1) 0; 0 0 1];
    v = [ sqrt(x3)*cos(2*pi*x2)  sqrt(x3)*sin(2*pi*x2) sqrt(1.0 - x3) ];
    H = eye(3) - 2 * v'*v;
    rM = -H*R;
else
    u = rand(1,3);
    u1 = u(1);
    u2 = u(2)*2*pi;
    u3 = u(3)*2*pi;
    t1 = sqrt(1-u1);
    t2 = sqrt(u1);
    qw = t1*sin(u2);
    qx = t1*cos(u2);
    qy = t2*sin(u3);
    qz = t2*cos(u3);  
    rM = [1- 2*(qy^2 + qz^2), 2*(qx*qy-qz*qw), 2*(qx*qz + qy*qw);
          2*(qx*qy+qz*qw), 1-2*(qx^2 + qz^2), 2*(qy*qz - qx*qw);
          2*(qx*qz-qy*qw), 2*(qy*qz + qx*qw), 1 - 2*(qx^2+qy^2)];
%     theta = 2*pi*rand(1);
%     rM = [cos(theta), -sin(theta);
%         sin(theta), cos(theta)];
%     
end
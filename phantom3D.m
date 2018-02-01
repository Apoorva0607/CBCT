function [P] = phantom3D(xoff,yoff,zoff,dx,dy,dz,Nx,Ny,Nz,E)

P = zeros(Nx,Ny,Nz);

x = -xoff*dx+(0:Nx-1)*dx;
y = -yoff*dy+(0:Ny-1)*dy;
z = -zoff*dz+(0:Nz-1)*dz;

[xcoord,ycoord,zcoord] = meshgrid(x,y,z);

for k=1:length(E(:,1))
    
    x0  = E(k,1);
    y0  = E(k,2);
    z0  = E(k,3);
    a   = E(k,4);
    b   = E(k,5);
    c   = E(k,6);
    ang1= E(k,7)*pi/180;
    ang2= E(k,8)*pi/180;
    ang3= E(k,9)*pi/180;
    f0  = E(k,10);
%     t   = E(k,11);
    
    t=0; %for elipsoid, t=1 for cylinder

    % ang1, ang2 and ang3 are three angles controls the direction of the
    % elipsoil (or cylinder). ang1 and ang2 are general spherical
    % coordinates of n3, and ang3 is the related angle between n1 and e1 (
    % or n2 and e2).
    
    %If the object is cylinder, we create it by first create a very long (
    %in n3 direction) ellipsoid and then cut it. c = 10*Nx*Ny*Nz makes it
    %big, and we read the clipping length into c0.
    
    if t==1 
        c0 = c;
        c  = 100*Nx*Ny*Nz;
    end
    
    e1 = [cos(ang2)*cos(ang1) sin(ang2)*cos(ang1) -sin(ang1)];
    e2 = [-sin(ang2) cos(ang2) 0];
    
    n1 = cos(ang3)*e1+sin(ang3)*e2;
    n2 = -sin(ang3)*e1+cos(ang3)*e2;
    n3 = [cos(ang2)*sin(ang1) sin(ang1)*sin(ang2) cos(ang1)];
    
    aa = 1/a^2;
    bb = 1/b^2;
    cc = 1/c^2;
    
    p1 = n1(1)^2*aa+n2(1)^2*bb+n3(1)^2*cc;
    p2 = n1(2)^2*aa+n2(2)^2*bb+n3(2)^2*cc;
    p3 = n1(3)^2*aa+n2(3)^2*bb+n3(3)^2*cc;
    p4 = n1(1)*n1(2)*aa+n2(1)*n2(2)*bb+n3(1)*n3(2)*cc;
    p5 = n1(1)*n1(3)*aa+n2(1)*n2(3)*bb+n3(1)*n3(3)*cc;
    p6 = n1(3)*n1(2)*aa+n2(3)*n2(2)*bb+n3(3)*n3(2)*cc;
    
    equation1 = p1*(xcoord-x0).^2+p2*(ycoord-y0).^2+p3*(zcoord-z0).^2;
    equation2 = p4*(xcoord-x0).*(ycoord-y0)+p5*(xcoord-x0).*(zcoord-z0)+p6*(ycoord-y0).*(zcoord-z0);
    equation = equation1+2*equation2;
    

    i      = find(equation<1.0);
    ell    = zeros(size(equation));
    ell(i) = 1;
 
    %If it is a cylinder, we cut it into half length c0.
    
    if t==1
        
        cylindereq = (xcoord-x0)*n3(1)+(ycoord-y0)*n3(2)+(zcoord-z0)*n3(3);
        cylindereq = abs(cylindereq)-c0;
        ii = find(cylindereq>=0);
        ell(ii)= 0;
        
    end
    
    P      = P+f0*ell;
    
end
        
    

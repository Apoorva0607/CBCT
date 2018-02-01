function G0=cbdataApoor(uoff, voff, du, dv, Nu, Nv, Np, R0, D, E)
G0  = zeros(Nv,Nu,Np);
u_1=((0:Nu-1) - uoff) * du;
v_1=((0:Nv-1) - voff) * dv;
[u , v]=meshgrid(u_1,v_1);
dlambda = 2*pi/(Np);
lambda=dlambda*(0:1:Np-1);
for i = 1:length(lambda)
    lam = lambda(i);
    
    for m = 1:size(E,1)
        
        
        x0  = E(m,1);
        y0  = E(m,2);
        z0  = E(m,3);
        a   = E(m,4);
        b   = E(m,5);
        c   = E(m,6);
        phi= E(m,7)*pi/180;
        psi= E(m,8)*pi/180;
        mu= E(m,9)*pi/180;
        f0  = E(m,10);
        
        
        e1 = [cos(psi)*cos(phi) sin(psi)*cos(phi) -sin(phi)]';
        e2 = [-sin(psi) cos(psi) 0]';
        alpha = zeros(Nv, Nu, 3);
        n1 = cos(mu)*e1+sin(mu)*e2;
        n2 = -sin(mu)*e1+cos(mu)*e2;
        n3 = [cos(psi)*sin(phi) sin(phi)*sin(psi) cos(phi)]';
        
        aa = 1/a.^2;
        bb = 1/b.^2;
        cc = 1/c.^2;
        
        A = [aa 0 0;
            0 bb 0;
            0 0 cc];
        
        Q = [n1';n2';n3'];
        
        E0= Q'*A*Q;
        
        ew = [cos(lam) sin(lam) 0]';
        eu = [-sin(lam) cos(lam) 0]';
        ev = [0 0 1]';
        for k = 1:3
            U=(u.*eu(k) + v*ev(k) - D*ew(k));
            alpha(:,:,k) = 1./sqrt(u.^2 + v.^2 + D.^2) .*U ;
        end
        alpha1 = alpha(:,:,1);
        alpha2 = alpha(:,:,2);
        alpha3 = alpha(:,:,3);
        ahat = [R0*cos(lam) R0*sin(lam) 0]';
        Xhat = ahat - [x0 y0 z0]';
        el0 = alpha1 .* ( alpha1.*E0(1,1) + alpha2.*E0(1,2) + alpha3.*E0(1,3) ) + ...
            alpha2 .* ( alpha1.*E0(2,1) + alpha2.*E0(2,2) + alpha3.*E0(2,3) ) + ...
            alpha3 .* ( alpha1.*E0(3,1) + alpha2.*E0(3,2) + alpha3.*E0(3,3) );
        
        el1 = alpha1 .* ( Xhat(1).*E0(1,1) + Xhat(2).*E0(1,2) + Xhat(3).*E0(1,3) ) + ...
            alpha2 .* ( Xhat(1).*E0(2,1) + Xhat(2).*E0(2,2) + Xhat(3).*E0(2,3) ) + ...
            alpha3 .* ( Xhat(1).*E0(3,1) + Xhat(2).*E0(3,2) + Xhat(3).*E0(3,3) );
        
        el2 = Xhat' * E0 * Xhat - 1;
        el00=1./el0;
        tp = (-el1 - real(sqrt(el1.^2-el0.*el2))).*el00;
        tq = (-el1 + real(sqrt(el1.^2-el0.*el2))).*el00;
        G0(:,:,i) = G0(:,:,i) + f0.*(tq-tp);
        
        
    end
end
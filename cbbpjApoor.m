function FilteredG = cbbpjApoor(G, uoff, voff, du, dv, R0, D, xoff, yoff, zoff, dx, dy, dz, Nx, Ny, Nz)
[Nv Nu Np]=size(G); %CORRECTED
u_1 = ((0:Nu-1) - uoff)*du;
v_1 = ((0:Nv-1) - voff)*dv;
[u, v] = meshgrid(u_1, v_1);

dLambda = (2*pi)/Np;
lambda=dLambda*(0:1:Np-1)';
%lambda = lambda(1:end-1); %CORRECTED
eu = [-sin(lambda) cos(lambda) zeros(length(lambda),1)];
ev = [zeros(length(lambda),1) zeros(length(lambda),1) ones(length(lambda),1)];
ew = [cos(lambda) sin(lambda) zeros(length(lambda),1)];

xc = ((0:Nx-1)-xoff) * dx;
yc = ((0:Ny-1)-yoff) * dy;
zc = ((0:Nz-1)-zoff) * dz;
[x, y] = meshgrid(xc, yc);

FilteredG = zeros(Ny, Nx, Nz); %CORRECTED

for i = 1:Nz
    z = zc(i);
    
    for j = 1:length(lambda)
        
        lam = lambda(j);
        alamx = R0 .* cos(lam) .* ones(size(x)); %CORRECTED
        alamy = R0 .* sin(lam) .* ones(size(y)); %CORRECTED
        ustar_num = (x-alamx).*eu(j,1) + (y-alamy).*eu(j,2) + z.*eu(j,3);
        ustar_den = (x-alamx).*ew(j,1) + (y-alamy).*ew(j,2) + z.*ew(j,3);
        ustar = -D.*(ustar_num./ustar_den);
        vstar_num = (x-alamx).*ev(j,1) + (y-alamy).*ev(j,2) + z.*ev(j,3);
        vstar_den = ustar_den;
        vstar = -D.*(vstar_num./vstar_den);
        Gmod = interp2(u, v, G(:,:,j), ustar, vstar, 'linear');
        Gmod(isnan(Gmod)) = 0;
        Dd = ( (alamx-x).*ew(j,1) + (alamy-y).*ew(j,2) + (0-z).*ew(j,3) ).^2;
        B0 = (R0.*D) ./ Dd;
        FilteredG(:,:,i) = FilteredG(:,:,i) + B0.*Gmod .* dLambda;
        
    end
    
end

FilteredG = FilteredG ./ 2;

end
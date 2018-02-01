function [F] = RampApoor(G,uoff,voff,du,dv,a,D)

[Nv,Nu,Nw]=size(G);

for deg=1:Nw
    for i=1:Nu
        for j=1:Nv
            
            u = du*(i-1-uoff);
            v = dv*(j-1-voff);
            
            
            for k=1:Nu
                
                u_1 = du*(k-1-uoff);
                dist  = D/sqrt(u_1^2+v^2+D^2)*du;
                
                T_1 = pi*(u-u_1)/du;
                T_2= T_1+pi;
                T_3= T_1-pi;
                
                if abs(T_1) < 1e-4
                    R1 = 0.5;
                else
                    R1 =sin(T_1)/T_1+(cos(T_1)-1)/(T_1^2);
                end
                
                if abs(T_2) < 1e-4
                    R2 = 0.5;
                else
                    
                  R2 = sin(T_2)/T_2+(cos(T_2)-1)/(T_2^2);
                end
                
                if abs(T_3) < 1e-4
                    R3 = 0.5;
                else
                    R3 = sin(T_3)/T_3+(cos(T_3)-1)/(T_3^2);
                end
                
                hram = a*R1+0.5*(1-a)*(R2+R3);
               
                hramp = 0.5/du^2*hram;
                
                g(k) = hramp*dist*G(j,k,deg);
                
            end
            
            F(j,i,deg) = sum(g);
            
        end
    end
end
end




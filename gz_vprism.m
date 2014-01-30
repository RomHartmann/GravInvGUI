function [gz] = gz_vprism(x1,y1,z1,x2,y2,z2,p,x0,y0,z0)

% gz_vprsim computes the vertical attraction due to a single vertical prism
% (x0,y0,z0) is the observation point
% prism extends from x1 to x2, y1 to y2, z1 to z2
% all directions are parallel to x, y, z axes
% p = density of rock in kg/m3
% distances are in km
% gz is in mGal


G=6.670e-11;

x=[0,0];
y=[0,0];
z=[0,0];

isign=[-1;1];

x(1)=x0-x1;
y(1)=y0-y1;
z(1)=z0-z1;
x(2)=x0-x2;
y(2)=y0-y2;
z(2)=z0-z2;


sum=0;
for i=1:2
    for j=1:2
        for k=1:2
            rijk = sqrt(x(i)^2 +y(j)^2 +z(k)^2);
            ijk = isign(i)*isign(j)*isign(k);
            arg1=atan((x(i)*y(j))/(z(k)*rijk));
            arg2=log(rijk+y(j));
            arg3=log(rijk+x(i));
            sum=sum+ijk*(z(k)*arg1-x(i)*arg2-y(j)*arg3);
        end
    end
end

gz=p*G*sum*1e5*1e3;



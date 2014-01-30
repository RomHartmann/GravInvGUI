function [calc] = Fmodel(xloc,yloc,D,p,xobs,yobs,Lx,Ly,dxsig,dysig,rh)

% Fmodel calculates the signal response for a single prism
% (x0,y0,z0) is the observation point
% prism extends from x1 to x2, y1 to y2, z1 to z2
% all directions are parallel to x, y, z axes
% distances are in km
% D is the depth of the prism
% p = density of rock in kg/m3
% xloc and yloc is the location of the prism on the x-y plane
% rh = reference height above the surface.  Minimum 1m

rh = rh/1000; %(convert to km)

x1=-Lx/2+xloc;
x2=Lx/2+xloc;
y1=-Ly/2+yloc;
y2=Ly/2+yloc;
z1=(D-D)+rh;
z2=D;

z0=0;


% preallocate calc
calc(1:xobs, 1:yobs) = 0;


% calculate the signal for an individual column    
for x0=0:dxsig:(dxsig*(xobs-1)); 
    xcount = round(1+(x0/dxsig));
    for y0=0:dysig:(dysig*(yobs-1));
        ycount = round(1+(y0/dysig));
        calc(xcount,ycount) = gz_vprism(x1,y1,z1,x2,y2,z2,p,x0,y0,z0);
    end
end




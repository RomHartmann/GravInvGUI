function mapplot(obsmap,calc0map,calcmap,xpos,ypos,Lx,Ly,depthmap,iter,x0,y0,x1,y1,dxsig,dysig,xobs,yobs,azi,elev) 
% plot maps of observed, inverted and depth distributions


% find max and min of obs and iteration for purpose of having same limits
% as prism plot
obsmax = max(max(obsmap));
obsmin = min(min(obsmap));
calc0max = max(max(calc0map));
calc0min = min(min(calc0map));
calcmax = max(max(calcmap));
calcmin = min(min(calcmap));

if obsmax>calc0max
    plt1max = obsmax;
else
    plt1max = calc0max;
end
if obsmax>calcmax
    plt2max = obsmax;
else
    plt2max = calcmax;
end

if obsmin<calc0min
    plt1min = obsmin;
else
    plt1min = calc0min;
end
if obsmin<calcmin
    plt2min = obsmin;
else
    plt2min = calcmin;
end


subplot(3,3,1); 
imagesc(0:dysig:dysig*yobs,0:dxsig:dxsig*xobs,calc0map,[plt1min,plt1max]); 
title('INITIAL ESTIMATE of the Gravity map')
xlabel('x-length (km)')
ylabel('y-length (km)')
zlabel('gz in mGal')
axis xy; axis equal; axis tight; 

subplot(3,3,4); 
imagesc(0:dysig:dysig*yobs,0:dxsig:dxsig*xobs,obsmap,[plt1min,plt1max]); 
title('OBSERVED Gravity map')
xlabel('x-length (km)')
ylabel('y-length (km)')
zlabel('gz in mGal')
axis xy; axis equal; axis tight; 
% colorbar('location','westoutside'):  Colorbar causes bug in plotting

subplot(3,3,7); 
imagesc(0:dysig:dysig*yobs,0:dxsig:dxsig*xobs,calcmap,[plt2min,plt2max]); 
title('INVERTED Gravity map')
xlabel('x-length (km)')
ylabel('y-length (km)')
zlabel('gz in mGal')
axis xy; axis equal; axis tight; 

% this mesh needs the y grid first, else dimensions do not agree. Strange.
subplot(3,3,[2,3]); grid on;

mesh((0:dysig:dysig*(yobs-1)),(0:dxsig:dxsig*(xobs-1)),(calcmap-obsmap))
% hold off; 
title(['DIFFERENCE between Observed and iterated signal at iteration number ', num2str(iter)])
xlabel('x-length (km)')
ylabel('y-length (km)')
zlabel('gz in mGal')
axis xy;
view(azi+90,elev+10);


subplot(3,3,[5,6,8,9]); grid on;
for xprismshow = 1:length(xpos)
    for yprismshow = 1:length(ypos)


        X1 = [xpos(xprismshow)-Lx/2,xpos(xprismshow)-Lx/2,xpos(xprismshow)-Lx/2,xpos(xprismshow)-Lx/2];
        Y1 = [ypos(yprismshow)+Ly/2,ypos(yprismshow)-Ly/2,ypos(yprismshow)-Ly/2,ypos(yprismshow)+Ly/2];
        Z1 = [0.001,0.001,-depthmap(xprismshow,yprismshow),-depthmap(xprismshow,yprismshow)];

        X2 = [xpos(xprismshow)-Lx/2,xpos(xprismshow)+Lx/2,xpos(xprismshow)+Lx/2,xpos(xprismshow)-Lx/2];
        Y2 = [ypos(yprismshow)+Ly/2,ypos(yprismshow)+Ly/2,ypos(yprismshow)+Ly/2,ypos(yprismshow)+Ly/2];
        Z2 = [0.001,0.001,-depthmap(xprismshow,yprismshow),-depthmap(xprismshow,yprismshow)];

        X3 = [xpos(xprismshow)+Lx/2,xpos(xprismshow)+Lx/2,xpos(xprismshow)+Lx/2,xpos(xprismshow)+Lx/2];
        Y3 = [ypos(yprismshow)+Ly/2,ypos(yprismshow)-Ly/2,ypos(yprismshow)-Ly/2,ypos(yprismshow)+Ly/2];
        Z3 = [0.001,0.001,-depthmap(xprismshow,yprismshow),-depthmap(xprismshow,yprismshow)];

        X4 = [xpos(xprismshow)-Lx/2,xpos(xprismshow)+Lx/2,xpos(xprismshow)+Lx/2,xpos(xprismshow)-Lx/2];
        Y4 = [ypos(yprismshow)-Ly/2,ypos(yprismshow)-Ly/2,ypos(yprismshow)-Ly/2,ypos(yprismshow)-Ly/2];
        Z4 = [0.001,0.001,-depthmap(xprismshow,yprismshow),-depthmap(xprismshow,yprismshow)];

        X5 = [xpos(xprismshow)-Lx/2,xpos(xprismshow)+Lx/2,xpos(xprismshow)+Lx/2,xpos(xprismshow)-Lx/2];
        Y5 = [ypos(yprismshow)+Ly/2,ypos(yprismshow)+Ly/2,ypos(yprismshow)-Ly/2,ypos(yprismshow)-Ly/2];
        Z5 = [-depthmap(xprismshow,yprismshow),-depthmap(xprismshow,yprismshow),-depthmap(xprismshow,yprismshow),-depthmap(xprismshow,yprismshow)];

        X6 = [xpos(xprismshow)-Lx/2,xpos(xprismshow)+Lx/2,xpos(xprismshow)+Lx/2,xpos(xprismshow)-Lx/2];
        Y6 = [ypos(yprismshow)+Ly/2,ypos(yprismshow)+Ly/2,ypos(yprismshow)-Ly/2,ypos(yprismshow)-Ly/2];
        Z6 = [0.001,0.001,0.001,0.001];

        patch(X1,Y1,-Z1,[0.6,0.3,0.7]);
        patch(X2,Y2,-Z2,[0.6,0.3,0.7]);
        patch(X3,Y3,-Z3,[0.6,0.3,0.7]);
        patch(X4,Y4,-Z4,[0.6,0.3,0.7]);
        patch(X5,Y5,-Z5,[0.6,0.3,1]);
        patch(X6,Y6,-Z6,[0.6,0.3,0.7]);

        hold on

    end
end

axis xy; set(gca,'YDir','reverse');
axis([min(xpos)-Lx,max(xpos)+Lx,min(ypos)-Ly,max(ypos)+Ly,0,1.1*max(max(depthmap))]); 
title(['Prisms depths depicted as heights at iteration number =  ', num2str(iter)])
ylabel('x-length (km)')
xlabel('y-length (km)')
zlabel('Depth of prisms (km)')
view(azi,elev);

hold off


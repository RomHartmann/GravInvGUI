function profileplots(xobs,yobs,obs,calc0,calc,xpos,ypos,Lx,Ly,depthmap,dxsig,dysig,iter,x0,y0,x1,y1) 
% plot profiles



% find max and min of obs and iteration for purpose of having same limits
% as prism plot
obsmax = max(obs);
obsmin = min(obs);
calc0max = max(calc0);
calc0min = min(calc0);
calcmax = max(calc);
calcmin = min(calc);

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


% find angle between (x0,y0) and(x1,y1) and the oblique intersection factor
t = atan((y1-y0)/(x1-x0));

% x and y coordinates of observed points
xobspos = 0:dxsig:(xobs-1)*dxsig;
yobspos = 0:dysig:(yobs-1)*dysig;
count = 0;
xyobs = zeros([length(xobspos)*length(yobspos) 2]);
for valx = xobspos
    for valy = yobspos
         count = count + 1;
         xyobs(count,1) = valx;
         xyobs(count,2) = valy;
    end
end


% x and y coordinates for interpolated points
if x1>=x0 && y1>=y0
    if floor(1e6*cos(t)) == 0   %Matlab rounding error
        xtheta = x0*ones(size(y0:dysig*sin(t):y1));  
    else
        xtheta = x0:dxsig*cos(t):x1;
    end
    if floor(1e6*sin(t)) == 0
        ytheta = y0*ones(size(x0:dxsig*cos(t):x1));
    else
        ytheta = y0:dysig*sin(t):y1;
    end
elseif x0>=x1 && y1>=y0
    if floor(1e6*cos(t+pi)) == 0   %Matlab rounding error
        xtheta = x0*ones(size(y0:dysig*sin(t+pi):y1));  
    else
        xtheta = x0:dxsig*cos(t+pi):x1;
    end
    if floor(1e6*sin(t+pi)) == 0
        ytheta = y0*ones(size(x0:dxsig*cos(t+pi):x1));
    else
        ytheta = y0:dysig*sin(t+pi):y1;
    end
elseif x1>=x0 && y0>=y1
    if floor(1e6*cos(t)) == 0   %Matlab rounding error
        xtheta = x0*ones(size(y0:dysig*sin(t):y1));  
    else
        xtheta = x0:dxsig*cos(t):x1;
    end
    if floor(1e6*sin(t)) == 0
        ytheta = y0*ones(size(x0:dxsig*cos(t):x1));
    else
        ytheta = y0:dysig*sin(t):y1;
    end
elseif x0>=x1 && y0>=y1
    if floor(1e6*cos(t+pi)) == 0   %Matlab rounding error
        xtheta = x0*ones(size(y0:dysig*sin(t+pi):y1));  
    else
        xtheta = x0:dxsig*cos(t+pi):x1;
    end
    if floor(1e6*sin(t+pi)) == 0
        ytheta = y0*ones(size(x0:dxsig*cos(t+pi):x1));
    else
        ytheta = y0:dysig*sin(t+pi):y1;
    end
end


interpcalc0 = zeros(size(xtheta));
interpcalc = zeros(size(xtheta));
interpobs = zeros(size(xtheta));
r = zeros(size(calc0));
% interpolate as weighted average for data
for interppoints = 1:length(xtheta)
    for points = 1:size(xyobs,1)
        r(points,1) = (sqrt((xtheta(interppoints)-xyobs(points,1))^2 + (ytheta(interppoints)-xyobs(points,2))^2))+0.00001;
    end
    interpcalc0(interppoints) = (sum(calc0./r.^3))/(sum(1./r.^3));
    interpcalc(interppoints) = (sum(calc./r.^3))/(sum(1./r.^3));
    interpobs(interppoints) = (sum(obs./r.^3))/(sum(1./r.^3));
end


subplot(3,1,1); 
plot3(xtheta,ytheta,interpcalc0,'- r'); 
hold on
plot3(xtheta,ytheta,interpobs,'- g'); 
hold off
if x1>=x0 && y1>=y0
    axis([x0-Lx,x1+Lx,y0-Ly,y1+Ly,plt1min, plt1max]); view(t*180/pi, 0)
elseif x0>=x1 && y1>=y0
    axis([x1-Lx,x0+Lx,y0-Ly,y1+Ly,plt1min, plt1max]); view(t*180/pi, 0)
elseif x1>=x0 && y0>=y1
    axis([x0-Lx,x1+Lx,y1-Ly,y0+Ly,plt1min, plt1max]); view(t*180/pi, 0)
elseif x0>=x1 && y0>=y1
    axis([x1-Lx,x0+Lx,y1-Ly,y0+Ly,plt1min, plt1max]); view(t*180/pi, 0)
end

% legends greened out due to bug that collapses the plot
title('INITIAL ESTIMATE (red) and OBSERVED (green) Gravity Profile')
xlabel('Distance from beginning of profile in km')
zlabel('gz in mGal')
% legend('location','Best')
% legend('Initial Estimate', 'Observed')




subplot(3,1,2); 
plot3(xtheta,ytheta,interpcalc,'- r'); 
hold on
plot3(xtheta,ytheta,interpobs,'- g'); 
hold off
if x1>=x0 && y1>=y0
    axis([x0-Lx,x1+Lx,y0-Ly,y1+Ly,plt2min, plt2max]); view(t*180/pi, 0)
elseif x0>=x1 && y1>=y0
    axis([x1-Lx,x0+Lx,y0-Ly,y1+Ly,plt2min, plt2max]); view(t*180/pi, 0)
elseif x1>=x0 && y0>=y1
    axis([x0-Lx,x1+Lx,y1-Ly,y0+Ly,plt2min, plt2max]); view(t*180/pi, 0)
elseif x0>=x1 && y0>=y1
    axis([x1-Lx,x0+Lx,y1-Ly,y0+Ly,plt2min, plt2max]); view(t*180/pi, 0)
end
title(['Gravity profile due to ITERATION (red) number ', num2str(iter)])
xlabel('Distance from beginning of profile in km')
zlabel('gz in mGal')
% legend('location','Best')
% legend('Iterated Plot', 'Observed')



subplot(3,1,3); grid on;
for xprismshow = 1:length(xpos)
    for yprismshow = 1:length(ypos)
        for obsdata = 1:length(xtheta)
            if xtheta(obsdata)<xpos(xprismshow)+Lx/2 && xtheta(obsdata)>xpos(xprismshow)-Lx/2 && ytheta(obsdata)<ypos(yprismshow)+Ly/2 && ytheta(obsdata)>ypos(yprismshow)-Ly/2

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
                
                patch(X1,Y1,Z1,[0.6,0.3,0.7]);
                patch(X2,Y2,Z2,[0.6,0.3,0.7]);
                patch(X3,Y3,Z3,[0.6,0.3,0.7]);
                patch(X4,Y4,Z4,[0.6,0.3,0.7]);
                patch(X5,Y5,Z5,[0.6,0.3,1]);
                patch(X6,Y6,Z6,[0.6,0.3,0.7]);
                
                hold on
            end
        end
    end
end

if x1>=x0 && y1>=y0
    axis([x0-Lx,x1+Lx,y0-Ly,y1+Ly,-1.1*max(max(depthmap)),0]); view(t*180/pi, 0)
elseif x0>=x1 && y1>=y0
    axis([x1-Lx,x0+Lx,y0-Ly,y1+Ly,-1.1*max(max(depthmap)),0]); view(t*180/pi, 0)
elseif x1>=x0 && y0>=y1
    axis([x0-Lx,x1+Lx,y1-Ly,y0+Ly,-1.1*max(max(depthmap)),0]); view(t*180/pi, 0)
elseif x0>=x1 && y0>=y1
    axis([x1-Lx,x0+Lx,y1-Ly,y0+Ly,-1.1*max(max(depthmap)),0]); view(t*180/pi, 0)
end
title(['Position plot of prisms with iteration number =  ', num2str(iter)])
xlabel('Distance from beginning of profile in km')
zlabel('Depth in km')
text(x0,y0,0.01,['(', num2str(x0), ', ', num2str(y0), ')']);
text(x1,y1,0.01,['(', num2str(x1), ', ', num2str(y1), ')']);
plot3([x0,x1],[y0,y1],[0.01,0.01],'k','LineWidth',2);
hold off



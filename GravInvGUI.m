function varargout = GravInvGUI(varargin)
% By Roman Hartmann
% Honours Project 2013
% % 
% GravInvGUI  generates the GUI and manages all subroutines involved in 
% the 3D inversion program
% % 
% % GRAVINVGUI M-file for GravInvGUI.fig
%      GRAVINVGUI, by itself, creates a new GRAVINVGUI or raises the existing
%      singleton*.
%
%      H = GRAVINVGUI returns the handle to a new GRAVINVGUI or the handle to
%      the existing singleton*.
%
%      GRAVINVGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRAVINVGUI.M with the given input arguments.
%
%      GRAVINVGUI('Property','Value',...) creates a new GRAVINVGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GravInvGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GravInvGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GravInvGUI

% Last Modified by GUIDE v2.5 19-Jan-2014 14:15:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GravInvGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GravInvGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before GravInvGUI is made visible.
function GravInvGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GravInvGUI (see VARARGIN)

% Choose default command line output for GravInvGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GravInvGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GravInvGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% 
% This function computes the inversion
% 
function pushbuttonInvert_Callback(hObject, eventdata, handles)

tic

% The inversion runs based in km.  All handles.km are there to change the
% interface to show in meters if checked, but the program still runs in km.
dxsig = handles.dxsig *handles.km;
dysig = handles.dysig *handles.km;
Lx = handles.Lx *handles.km;
Ly = handles.Ly *handles.km;
bhyes = handles.Boreholeyes;
k1 = handles.kay1;
k2 = handles.kay2;
n = handles.n;
x0 = handles.x0*handles.km;
x1 = handles.x1*handles.km;
y0 = handles.y0*handles.km;
y1 = handles.y1*handles.km;
azi = handles.azi;
elev = handles.elev;
MsftD = handles.msftdisp;
MapD = handles.dispmap;
ProfD = handles.dispprof;
rh = handles.rh; %reference height

try

    % create progress bars
    waitinv = waitbar(0,'0%', 'position', [400, 400, 280, 50]);
    waitit = waitbar(0, ['1/' num2str(n), ' iterations'], 'position', [400, 480, 280, 50]);

    % 
    % inversion code begin

    % Call the observed model.
    xobs = handles.xobs;
    yobs = handles.yobs;
    obsxls = handles.obsmap;

    obs = reshape(obsxls,(xobs*yobs),1);   % turn 2D data matrix into vector

    % Original guess of prism depth distribution 

    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;


    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

    % creates permutations of x and y coordinates as x-y pairs.
    count = 0;
    permxy = zeros([length(xpos)*length(ypos) 2]);
    for valx = xpos
        for valy = ypos
             count = count + 1;
             permxy(count,1) = valx;
             permxy(count,2) = valy;
        end
    end

    % Initial depth and density information added to permxy matrix
    depthdist = handles.depthdist*handles.km;

    depth = reshape(depthdist,(length(xpos)*length(ypos)),1);
    origdepth = depth;

    
    densitydist = handles.densitydist;
    density = reshape(densitydist,(length(xpos)*length(ypos)),1);


    [nrprisms,~] = size(permxy);
    dz = 0.01; % change in depth to create dP


    % Borehole info
    if bhyes == 1
            try
                [BHFileName,BHPathName] = uigetfile('*.xlsx','Please select Borehole Excel Data');
                BHData = strcat(BHPathName,BHFileName);
                bhxls = xlsread(BHData);
                [bhrow,~]=size(bhxls);
                bhx = bhxls(:,1); 
                bhy = bhxls(:,2);
                bhd = bhxls(:,3)+0.00005;    % add 5cm to depth to prevent NaN errors
            catch
                BHerr = msgbox('Import of Borehole data aborted.  Inversion will continue', 'Import Error' ,'warn');
                movegui(BHerr, 'center')
                bhyes = 0;
            end
    end


    % preallocation
    A = zeros(length(obs),nrprisms);
    msft(1:n) = 0;
    dP = zeros(nrprisms,1);



    for iter = 1:n
        calc0 = zeros(size(obs));
        calc = zeros(size(obs));

        for pn = 1:nrprisms
            xloc = permxy(pn,1);
            yloc = permxy(pn,2);
            D0 = origdepth(pn);
            D = depth(pn);
            p = density(pn);

            c0 = Fmodel(xloc,yloc,D0,p,xobs,yobs,Lx,Ly,dxsig,dysig,rh);  % original forward model
            c0rs = reshape(c0,(xobs*yobs),1);       %c0 reshaped
            newD = D;

            newD = newD + dP(pn);

            % Borehole data is fixed depth data
            if bhyes == 1
                for bhid = 1:bhrow
                    if bhx(bhid)>xloc-Lx/2 && bhx(bhid)<xloc+Lx/2 && bhy(bhid)>yloc-Ly/2 && bhy(bhid)<yloc+Ly/2
                        newD = bhd(bhid);
                    end
                end
            end


    %       iterated forward models
            calcone = Fmodel(xloc,yloc,newD,p,xobs,yobs,Lx,Ly,dxsig,dysig,rh);
            calctwo = Fmodel(xloc,yloc,(newD+dz),p,xobs,yobs,Lx,Ly,dxsig,dysig,rh);
            calconers = reshape(calcone,(xobs*yobs),1);
            calctwors = reshape(calctwo,(xobs*yobs),1);
            A(:,pn) = (calctwors - calconers)/dz;

            c1rs = reshape(calcone,(xobs*yobs),1);

            calc0 = calc0 + c0rs;
            calc = calc + c1rs;
            depth(pn) = newD;

            % creates a progress bar
            waitbar((pn/nrprisms),waitinv, [num2str(100*(pn/nrprisms)), '%']);
            waitbar((iter/n),waitit, [num2str(iter), '/', num2str(n), ' iterations']);

        end


        E = obs - calc;


        % k is the inversion dampening factor
        k = k1*10^k2;



        dP = ((A' * A)+k*eye(size(A' * A)))\(A' * E);


        for pn2 = 1:nrprisms
            if dP(pn2) + depth(pn2) < 0
                while dP(pn2) + depth(pn2) < 0;        %  preventing negative depths
                   dP(pn2) = dP(pn2)/2 ;
                end
            end
        end


    %   transform vector back into a 2d matrix for plotting and mapping
        obsmap = reshape(obs,xobs,yobs);
        calc0map = reshape(calc0,xobs,yobs);
        calcmap = reshape(calc,xobs,yobs);
        depthmap = (reshape(depth,length(ypos),length(xpos))'/handles.km);

        proferrid = 0;  %to display an error only once at the end
        if ProfD == 1
            if xobs ==1 || yobs ==1
                f1 = figure(1);
                set(f1,'Name','Position Plot')
                movegui(f1,'northwest')
                pause(0.1); cla; 
                profileplots(xobs,yobs,obs,calc0,calc,xpos/handles.km,ypos/handles.km,Lx/handles.km,Ly/handles.km,depthmap,dxsig/handles.km,dysig/handles.km,iter,x0/handles.km,y0/handles.km,x1/handles.km,y1/handles.km);
            else
            proferrid = 1;
            end
        end

        maperrid = 0;
        if MapD == 1
            if xobs > 1 && yobs > 1 
                f2 = figure (2); 
                set(f2,'Name','Map Plot')
                movegui(f2,'southwest')
                pause(0.1); cla;
                mapplot(obsmap,calc0map,calcmap,xpos/handles.km,ypos/handles.km,Lx/handles.km,Ly/handles.km,depthmap,iter,x0/handles.km,y0/handles.km,x1/handles.km,y1/handles.km,dxsig/handles.km,dysig/handles.km,xobs,yobs,azi,elev)
            else
            maperrid = 1;
            end
        end

        msft(iter) = misfit(obs,calc);

    %     changes k1 every iteration by a factor of 10 depending on misfit
        if handles.dynk == 1;
            if iter>1
                 if msft(iter) > msft(iter-1)
                     k1 = 10*k1;
                 else
                     k1 = 0.1*k1;
                 end
            end
        end

    end


    if MsftD ==1
        f3 = figure (3);
        set(f3,'Name','Misfit Plot')
        movegui(f3,'north')
        semilogy((1:n),msft);
        axis([0.9,n,0,max(msft)])
        title(['MISFIT due to ', num2str(n), ' iterations'])
        xlabel('Iteration number')
        ylabel('Misfit value')
        [~,mi]=min(msft);
        text(mi-0.5,1.1*min(msft),num2str(min(msft)));
    end
    % 
    % Inversion code end
    % Display code begin
    % 

    close(waitinv);
    close(waitit);

    if handles.dispobs == 1;
        f4 = figure (4);
        set(f4,'Name','Depth distributions of inverted prisms')
        movegui(f4,'southeast')
        uitable('Data', flipud(depthmap),'Position',[10 10 500 350]);
    end
    % Export and display Depths and maps
    % Exports the Depth values into an excel sheet
    if handles.expobs ==1;
        try
            [expobsFileName,expobsPathName] = uiputfile('*.xlsx','Save Depthmap Data');
            expdepthpath = strcat(expobsPathName,expobsFileName);
            xlswrite(expdepthpath,flipud(depthmap));
        catch
            expoerr = msgbox('Inverted depths data export aborted', 'Export aborted' ,'warn');
            movegui(expoerr, 'center') 
        end
    end

    % Exports the inverted signal into an excel sheet
    if handles.expinv ==1;
        try
            [expinvFileName,expinvPathName] = uiputfile('*.xlsx','Save Inverted Signal Data');
            expinvpath = strcat(expinvPathName,expinvFileName);
            xlswrite(expinvpath,calcmap);
        catch
            expierr = msgbox('Inverted signal data export aborted', 'Export aborted', 'warn');
            movegui(expierr, 'center')
        end
    end

    % Plots the Observed data in a seperate figure via checkbox
    % converts the map into a profile if the data are 1D
    if xobs > 1 && yobs > 1
        if handles.expobsmap ==1;
            f5 = figure (5);
            set(f5,'Name','OBSERVED Data Map')
            movegui(f5,'west')
            imagesc(0:dysig:dysig*yobs,0:dxsig:dxsig*xobs,obsmap); 
            cb1 = colorbar('peer',gca);
            set(get(cb1,'xlabel'),'String', 'mGal');
            title('OBSERVED Gravity map')
            xlabel('x-length (km)')
            ylabel('y-length (km)')
            zlabel('gz in mGal')
            axis xy; axis equal; axis tight;
        end
    else
        if xobs == 1 && handles.expobsmap ==1;
            f5 = figure (5);
            set(f5,'Name','OBSERVED Data Map')
            movegui(f5,'west')
            plot(0:dysig:dysig*(yobs-1), obsmap);
            title('OBSERVED Data Profile');
            xlabel('y-length (km)');
            ylabel('gz in mGal');
            axis xy; 
        elseif yobs == 1 && handles.expobsmap ==1;
            f5 = figure (5);
            set(f5,'Name','OBSERVED Data Map')
            movegui(f5,'west')
            plot(0:dxsig:dxsig*(xobs-1), obsmap);
            title('OBSERVED Data Profile');
            xlabel('y-length (km)');
            ylabel('gz in mGal');
            axis xy; 
        end
    end

    
    % Plots the Modelled data in a seperate figure via checkbox
    % converts the map into a profile if the data are 1D
    if xobs > 1 && yobs > 1
        if handles.expmodmap ==1;
            f6 = figure(6);
            set(f6,'Name','MODELLED Data Map')
            movegui(f6,'center')
            imagesc(0:dysig:dysig*yobs,0:dxsig:dxsig*xobs,calcmap); 
            cb2 = colorbar('peer',gca);
            set(get(cb2,'xlabel'),'String', 'mGal');
            title(['Gravity profile due to ITERATION number ', num2str(iter)])
            xlabel('x-length (km)')
            ylabel('y-length (km)')
            zlabel('gz in mGal')
            axis xy; axis equal; axis tight;
        end
    else
        if xobs == 1 && handles.expmodmap ==1;
            f6 = figure (6);
            set(f6,'Name','MODELLED Data Map')
            movegui(f6,'center')
            plot(0:dysig:dysig*(yobs-1), calcmap);
            title('MODELLED Data Profile');
            xlabel('y-length (km)');
            ylabel('gz in mGal');
            axis xy; 
        elseif yobs == 1 && handles.expmodmap ==1;
            f6 = figure (6);
            set(f6,'Name','MODELLED Data Map')
            movegui(f6,'center')
            plot(0:dxsig:dxsig*(xobs-1), calcmap);
            title('MODELLED Data Profile');
            xlabel('y-length (km)');
            ylabel('gz in mGal');
            axis xy; 
        end
    end

    
    %     plots an interpolated surface plot of the depth estimates
    if handles.dispsurfmap == 1;
      if xobs > 1 && yobs > 1
            f7 = figure (7);
            set(f7,'Name','Interpolated depth map')
            movegui(f7,'east')
            dm=-depthmap;
            interpdm=interp2(dm,handles.interpsurf);
            

            dxsig = handles.dxsig *handles.km;
            dysig = handles.dysig *handles.km;
            Lx = handles.Lx *handles.km;
            Ly = handles.Ly *handles.km;
            
            
            intf = handles.interpsurf;
            % calculates the interpolation factor
            count = 0;
            fn = 0;
            while count < intf
                fn = (2*fn) + 1;
                count = count + 1;
            end
            
            
            intx = -xpadding*Lx:Lx/(1+fn):(floor(((xobs-1)*dxsig)/Lx))*Lx + xpadding*Lx;
            inty = -ypadding*Ly:Ly/(1+fn):(floor(((yobs-1)*dysig)/Ly))*Ly + ypadding*Ly;
            
            
            surf(inty, intx, interpdm);
            view(azi+75,elev+20); %Defined by the azi and elev under the map plot options
            cb3 = colorbar('peer',gca);
            if handles.km == 1;
                set(get(cb3,'xlabel'),'String', 'km');
            else
                set(get(cb3,'xlabel'),'String', 'm');
            end
            title('Interpolated depth map')
            xlabel('x-length (km)')
            ylabel('y-length (km)')
            zlabel('Depth (km)')
            axis([-ypadding*Ly (dysig*(yobs-1)+Ly*ypadding) -xpadding*Lx (dxsig*(xobs-1)+Lx*xpadding) min(min(dm)) 0])
            axis xy; 

        else % display an error if the dataset is a profile
            surferr = msgbox('The interpolated surface plot is 1 dimensional', 'Map interpolation not successful', 'error');
            movegui(surferr, 'center')
      end
    end
    
%   In order to display errors only once 
    if proferrid == 1;
        proferr = msgbox('Profiles can only be created for 1D data', 'Profile not plot' ,'warn');
        movegui(proferr, 'center')
    end
    if maperrid == 1;
        maperr = msgbox('Maps cannot be created for 1D data', 'Map not plot', 'warn');
        movegui(maperr, 'center')
    end
    
    
    catch  % Display error message instead of crashing program
    try
    close(waitinv);
    close(waitit);
    catch
        toc;
        inverr = msgbox(['Inversion  was aborted after ', num2str(toc), ' seconds'], 'Operation ended' ,'warn');
        movegui(inverr, 'center')
    end
end
    
tocinv = (round(100*toc))/100;
set(handles.textinvMinutes, 'String', num2str(floor((tocinv)/60)))
set(handles.textinvSeconds, 'String', num2str((tocinv) - 60*floor((tocinv)/60)))



guidata(hObject,handles)
% 
% 
% 



% Creates the misfit plot
function [misfit] = misfit(obs,calc)
misfit = 0;
for ndata = 1:length(calc);
    dmisfit = (obs(ndata)-calc(ndata))^2;
    misfit = misfit + dmisfit;
end


% Toggles depth table as editable
function togglebutton2_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    set(handles.uitable2, 'ColumnEditable', true)
else
    set(handles.uitable2, 'ColumnEditable', false)
end


% Toggles density table as editable
function togglebutton1_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    set(handles.uitable1, 'ColumnEditable', true)
else
    set(handles.uitable1, 'ColumnEditable', false)
end


% x0 in 'Split Dataset'
function editsplitx0_Callback(hObject, eventdata, handles)

x0split = str2double(get(hObject,'String'));

if x0split < 0
    x0split = 0;
    set(handles.editsplitx0,'String',x0split)
end
if x0split > handles.x1split
    x0split = handles.x1split;
end

set(handles.editsplitx0,'String',x0split);
handles.x0split = x0split;

yobs = handles.yobs;
dysig = handles.dysig *handles.km;
axes(handles.axesObs); 


plot([0,dysig*(yobs-1)/handles.km],[x0split,x0split],'k','LineWidth',2);
hold on;

% set the maximum profile plot length
handles.x1 = handles.x1split - handles.x0split;
set(handles.edit13,'String',handles.x1)


guidata(hObject,handles)


function editsplitx0_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.x0split = 0;
set(hObject,'String',handles.x0split)

guidata(hObject, handles)

% x1 in 'Split Dataset'
function editsplitx1_Callback(hObject, eventdata, handles)

x1split = str2double(get(hObject,'String'));

xobs = handles.xobs;
yobs = handles.yobs;
dxsig = handles.dxsig *handles.km;
dysig = handles.dysig *handles.km;
axes(handles.axesObs); 

if x1split > dxsig*(xobs-1)/handles.km
    x1split = dxsig*(xobs-1)/handles.km;
end
if x1split < 0
    x1split = 0;
end
if x1split < handles.x0split
    x1split = handles.x0split;
end
set(handles.editsplitx1,'String',x1split)


handles.x1split = x1split;


plot([0,dysig*(yobs-1)/handles.km],[x1split,x1split],'k','LineWidth',2);
hold on;

% set the maximum profile plot length
handles.x1 = handles.x1split - handles.x0split;
set(handles.edit13,'String',handles.x1)

guidata(hObject,handles)


function editsplitx1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% y0 in 'Split Dataset'
% it is y and not x because of the matrix mismatches
function editsplity0_Callback(hObject, eventdata, handles)

y0split = str2double(get(hObject,'String'));

if y0split < 0
    y0split = 0;
    set(handles.editsplity0,'String',handles.y0split)
end
if y0split > handles.y1split
    y0split = handles.y1split;
end

set(handles.editsplity0,'String',y0split);
handles.y0split = y0split;

xobs = handles.xobs;
dxsig = handles.dxsig *handles.km;
axes(handles.axesObs); 


plot([y0split,y0split],[0,dxsig*(xobs-1)/handles.km],'k','LineWidth',2);
hold on;

% set the maximum profile plot length
handles.y1 = handles.y1split - handles.y0split;
set(handles.edit15,'String',handles.y1)

guidata(hObject,handles)


function editsplity0_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.y0split = 0;
set(hObject,'String',handles.y0split)

guidata(hObject, handles)

% y1 in 'Split Dataset'
function editsplity1_Callback(hObject, eventdata, handles)

y1split = str2double(get(hObject,'String'));


xobs = handles.xobs;
yobs = handles.yobs;
dxsig = handles.dxsig *handles.km;
dysig = handles.dysig *handles.km;
axes(handles.axesObs); 

if y1split > dysig*(yobs-1)/handles.km
    y1split = dysig*(yobs-1)/handles.km;
end
if y1split < 0
    y1split = 0;
end
if y1split < handles.y0split
    y1split = handles.y0split;
end
set(handles.editsplity1,'String',y1split)


handles.y1split = y1split;

plot([y1split,y1split],[0,dxsig*(xobs-1)/handles.km],'k','LineWidth',2);
hold on;

% set the maximum profile plot length
handles.y1 = handles.y1split - handles.y0split;
set(handles.edit15,'String',handles.y1)

guidata(hObject,handles)


function editsplity1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Split Dataset checkbox
function checkboxSplit_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    set(handles.editsplitx0, 'visible', 'on')
    set(handles.editsplitx1, 'visible', 'on')
    set(handles.editsplity0, 'visible', 'on')
    set(handles.editsplity1, 'visible', 'on')
    set(handles.textsplitx0, 'visible', 'on')
    set(handles.textsplitx1, 'visible', 'on')
    set(handles.textsplity0, 'visible', 'on')
    set(handles.textsplity1, 'visible', 'on')
    set(handles.popupmenuObs, 'visible', 'on')
    set(handles.axesObs, 'visible', 'on')
    set(handles.pushbuttoncalcobs, 'visible', 'on')
    set(handles.pushbutton9, 'visible', 'on')
    handles.split = 1;
else
    set(handles.editsplitx0, 'visible', 'off')
    set(handles.editsplitx1, 'visible', 'off')
    set(handles.editsplity0, 'visible', 'off')
    set(handles.editsplity1, 'visible', 'off')
    set(handles.textsplitx0, 'visible', 'off')
    set(handles.textsplitx1, 'visible', 'off')
    set(handles.textsplity0, 'visible', 'off')
    set(handles.textsplity1, 'visible', 'off')
    set(handles.pushbutton9, 'visible', 'off')
    handles.split = 0;
end


% aliasing multiple in x direction
function edit22_Callback(hObject, eventdata, handles)

aliasx = str2double(get(hObject,'String'));

if aliasx <1
    aliasx = 1; %negative number of iterations not allowed
end
aliasx = round(aliasx);
set(hObject,'String',aliasx)

handles.aliasx = aliasx;


guidata(hObject, handles)



function edit22_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.aliasx = 1;
set(hObject,'String',num2str(handles.aliasx))
guidata(hObject, handles)


% 3D depth inversion figure checkbox
function checkbox7_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    handles.dispobs = 1;
else
    handles.dispobs = 0;
end
guidata(hObject, handles)



% Export 3D depth values as xcel checkbox
function checkbox8_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    handles.expobs = 1;
else
    handles.expobs = 0;
end
guidata(hObject, handles)


% Dsiplay profile checkbox
function checkbox11_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    set(handles.edit12, 'visible', 'on')
    set(handles.edit13, 'visible', 'on')
    set(handles.edit14, 'visible', 'on')
    set(handles.edit15, 'visible', 'on')
    set(handles.text12, 'visible', 'on')
    set(handles.text13, 'visible', 'on')
    set(handles.text14, 'visible', 'on')
    set(handles.text15, 'visible', 'on')
    handles.dispprof = 1;
else
    set(handles.edit12, 'visible', 'off')
    set(handles.edit13, 'visible', 'off')
    set(handles.edit14, 'visible', 'off')
    set(handles.edit15, 'visible', 'off')
    set(handles.text12, 'visible', 'off')
    set(handles.text13, 'visible', 'off')
    set(handles.text14, 'visible', 'off')
    set(handles.text15, 'visible', 'off')
    handles.dispprof = 0;
end
guidata(hObject, handles)



% x0 for profile plotting
function edit12_Callback(hObject, eventdata, handles)

x0 = str2double(get(hObject,'String'));
handles.x0 = x0;
guidata(hObject, handles)



function edit12_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% x1 for profile plotting
function edit13_Callback(hObject, eventdata, handles)

x1 = str2double(get(hObject,'String'));
handles.x1 = x1;
guidata(hObject, handles)



function edit13_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% y0 for profile plotting
function edit14_Callback(hObject, eventdata, handles)

y0 = str2double(get(hObject,'String'));
handles.y0 = y0;
guidata(hObject, handles)



function edit14_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% y1 for profile plotting
function edit15_Callback(hObject, eventdata, handles)

y1 = str2double(get(hObject,'String'));
handles.y1 = y1;
guidata(hObject, handles)



function edit15_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Checkbox for Map display
function checkbox12_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    set(handles.edit16, 'visible', 'on')
    set(handles.edit17, 'visible', 'on')
    set(handles.text16, 'visible', 'on')
    set(handles.text17, 'visible', 'on')
    handles.dispmap = 1;
else
    set(handles.edit16, 'visible', 'off')
    set(handles.edit17, 'visible', 'off')
    set(handles.text16, 'visible', 'off')
    set(handles.text17, 'visible', 'off')
    handles.dispmap = 0;
end
handles.azi = 290;
guidata(hObject, handles)


% Azimuth for map plot
function edit16_Callback(hObject, eventdata, handles)

azi = str2double(get(hObject,'String'));
handles.azi = azi;
guidata(hObject, handles)



function edit16_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.azi = 290;
set(hObject,'String',handles.azi)
guidata(hObject, handles)


% Elevation for map plot
function edit17_Callback(hObject, eventdata, handles)

elev = str2double(get(hObject,'String'));
handles.elev = elev;
guidata(hObject, handles)



function edit17_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.elev = 35;
set(hObject,'String',num2str(handles.elev))
guidata(hObject, handles)


% Checkbox for Borehole information
function checkboxBH_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    handles.Boreholeyes = 1;
else
    handles.Boreholeyes = 0;
end
guidata(hObject,handles)


% Fetches the Observed data and displays path
function pushbuttonObs_Callback(hObject, eventdata, handles)


handles.dxsig = 1;  %resets the data values such that changes do not carry
handles.dysig = 1;
handles.Lx = 1;
handles.Ly = 1;
handles.aliasx = 1;
handles.aliasy = 1;

set(handles.edit8,'String',handles.dxsig)
set(handles.edit1,'String',handles.dysig)
set(handles.edit2,'String',handles.Lx)
set(handles.edit3,'String',handles.Ly)
set(handles.edit22,'String',handles.aliasx)
set(handles.edit30,'String',handles.aliasy)

% hide the invert button
set(handles.pushbuttonInvert, 'visible', 'off')
set(handles.text37, 'visible', 'on')
set(handles.pushbutton6, 'visible', 'off')
set(handles.text35, 'visible', 'on')
handles.dispbuttdens = 0;
handles.dispbuttdepth = 0;

try
    [ObsFileName,ObsPathName] = uigetfile('*.xlsx','Please select Observed Excel Data');
    ObsData = strcat(ObsPathName,ObsFileName);

    handles.observeddata = ObsData;
    set(handles.editObs,'String',ObsData)

    
    obsmap = xlsread(ObsData);


    handles.obsmap = obsmap;
    [xobs, yobs] = size(obsmap);
    handles.xobs = xobs;
    handles.yobs = yobs;

    dxsig = handles.dxsig*handles.km;
    dysig = handles.dysig*handles.km;
    Lx = handles.Lx*handles.km;
    Ly = handles.Ly*handles.km;
    
   
    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;
    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

    %  computation time estimate
    estFmTime = (10^((0.625*log10(xobs*yobs*length(xpos)*length(ypos)))-2.2))/2.75;
    estFmTime3sf = (round(100*estFmTime)/100);

    tmin = floor((handles.n*2.75*estFmTime3sf)/60);
    tsec = (handles.n*2.75*estFmTime3sf) - 60*tmin;
    handles.tmin = tmin;
    handles.tsec = tsec;
    
    set(handles.textMinutes, 'String', num2str(tmin))
    set(handles.textSeconds, 'String', num2str(tsec))
    set(handles.text32, 'String', estFmTime3sf)
    
    
    handles.x0split = 0;
    handles.y0split = 0;
    handles.x1split = dxsig*(xobs-1)/handles.km;
    handles.y1split = dysig*(yobs-1)/handles.km;
    % create duplicate to prevent value escalation
    handles.x1split0 = dxsig*(xobs-1)/handles.km;
    handles.y1split0 = dysig*(yobs-1)/handles.km;
    set(handles.editsplitx1,'String',handles.x1split)
    set(handles.editsplity1,'String',handles.y1split)
    set(handles.editsplitx0,'String',handles.x0split)
    set(handles.editsplity0,'String',handles.y0split)

    % Profile plot Standard values
    set(handles.edit12,'String',0)
    set(handles.edit14,'String',0)
    set(handles.edit13,'String',handles.x1split)
    set(handles.edit15,'String',handles.y1split)
    handles.x0 = 0;
    handles.x1 = handles.x1split;
    handles.y0 = 0;
    handles.y1 = handles.y1split;

    if xobs > 1 && yobs > 1
        obs_list{1} = 'Scaled Image';
        obs_list{2} = 'Filled Contour (10)';
        obs_list{3} = '3D Surface';
        set(handles.popupmenuObs,'String',obs_list)
    else
        obs_list{1} = 'Profile';
        set(handles.popupmenuObs,'String',obs_list)
    end

catch  % Display error message instead of crashing program
    patherr = msgbox('Please enter a valid path', 'Path Error' ,'error');
    movegui(patherr, 'center')
    
end

guidata(hObject,handles)


% Observed data path is displayed and text editable
function editObs_Callback(hObject, eventdata, handles)

handles.dxsig = 1;  %resets the data values such that changes do not carry
handles.dysig = 1;
handles.Lx = 1;
handles.Ly = 1;
handles.aliasx = 1;
handles.aliasy = 1;

set(handles.edit8,'String',handles.dxsig)
set(handles.edit1,'String',handles.dysig)
set(handles.edit2,'String',handles.Lx)
set(handles.edit3,'String',handles.Ly)
set(handles.edit22,'String',handles.aliasx)
set(handles.edit30,'String',handles.aliasy)

% hide the invert button
set(handles.pushbuttonInvert, 'visible', 'off')
set(handles.text37, 'visible', 'on')
set(handles.pushbutton6, 'visible', 'off')
set(handles.text35, 'visible', 'on')
handles.dispbuttdens = 0;
handles.dispbuttdepth = 0;

try

    ObsData = get(hObject,'String');
    handles.observeddata = ObsData;

    obsmap = xlsread(ObsData);
        
    [xobs, yobs] = size(obsmap);

    handles.obsmap = obsmap;
    handles.xobs = xobs;
    handles.yobs = yobs;

    % for the time estimate
    dxsig = handles.dxsig*handles.km;
    dysig = handles.dysig*handles.km;
    Lx = handles.Lx*handles.km;
    Ly = handles.Ly*handles.km;
    
    xobs = handles.xobs;
    yobs = handles.yobs;
    
    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;
    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

    estFmTime = (10^((0.625*log10(xobs*yobs*length(xpos)*length(ypos)))-2.2))/2.75;
    estFmTime3sf = (round(100*estFmTime)/100);

    tmin = floor((handles.n*2.75*estFmTime3sf)/60);
    tsec = (handles.n*2.75*estFmTime3sf) - 60*tmin;
    handles.tmin = tmin;
    handles.tsec = tsec;
    
    set(handles.textMinutes, 'String', num2str(tmin))
    set(handles.textSeconds, 'String', num2str(tsec))
    set(handles.text32, 'String', estFmTime3sf) 
    
    
    handles.x0split = 0;
    handles.y0split = 0;
    handles.x1split = dxsig*(xobs-1)/handles.km;
    handles.y1split = dysig*(yobs-1)/handles.km;
    % create duplicate to prevent value escalation
    handles.x1split0 = dxsig*(xobs-1)/handles.km;
    handles.y1split0 = dysig*(yobs-1)/handles.km;
    set(handles.editsplitx1,'String',handles.x1split)
    set(handles.editsplity1,'String',handles.y1split)
    set(handles.editsplitx0,'String',handles.x0split)
    set(handles.editsplity0,'String',handles.y0split)

    % Profile plot Standard values
    set(handles.edit12,'String',0)
    set(handles.edit14,'String',0)
    set(handles.edit13,'String',handles.x1split)
    set(handles.edit15,'String',handles.y1split)
    handles.x0 = 0;
    handles.x1 = handles.x1split;
    handles.y0 = 0;
    handles.y1 = handles.y1split;
    
    if xobs > 1 && yobs > 1
        obs_list{1} = 'Scaled Image';
        obs_list{2} = 'Filled Contour (10)';
        obs_list{3} = '3D Surface';
        set(handles.popupmenuObs,'String',obs_list)
    else
        obs_list{1} = 'Profile';
        set(handles.popupmenuObs,'String',obs_list)
    end 
    
    
catch  % Display error message instead of crashing program
    patherr = msgbox('Please enter a valid path', 'Path Error' ,'error');
    movegui(patherr, 'center')
end

guidata(hObject,handles)


function editBH_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% sets the y signal station spacing
function edit1_Callback(hObject, eventdata, handles)

dy_sig = str2double(get(hObject,'String'));

if dy_sig < 0
    dy_sig = -dy_sig; %negative distance not allowed
end
if dy_sig == 0
    dy_sig = 1;
end
if dy_sig > handles.Ly % prevents underdetermined problem
    handles.Ly = dy_sig;
    set(handles.edit3, 'String', dy_sig)
end

set(hObject,'String',dy_sig)

handles.dysig = dy_sig;

handles.y1split = handles.y1split*handles.dysig;
handles.y1 = handles.y1*dy_sig;
handles.y0split = handles.y0split*handles.dysig;
handles.y0 = handles.y0*dy_sig;

set(handles.edit15,'String',handles.y1)
set(handles.editsplity1,'String',handles.y1split)
set(handles.edit14,'String',handles.y0)
set(handles.editsplity0,'String',handles.y0split)

dysig = handles.dysig*handles.km;
Ly = handles.Ly*handles.km;
yobs = handles.yobs;
ypadding = handles.paddingprismsy;
ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

handles.totalprismsy = length(ypos);

% display the number of prisms in the textbox
set(handles.edit36, 'string', handles.totalprismsy);

% hide the invert button
set(handles.pushbuttonInvert, 'visible', 'off')
set(handles.text37, 'visible', 'on')
set(handles.pushbutton6, 'visible', 'off')
set(handles.text35, 'visible', 'on')
handles.dispbuttdens = 0;
handles.dispbuttdepth = 0;

guidata(hObject,handles)



function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.dysig= 1;
set(hObject,'String',handles.dysig)
guidata(hObject, handles)


% sets the x-length per prism
function edit2_Callback(hObject, eventdata, handles)

L_x = str2double(get(hObject,'String'));

if L_x < 0
    L_x = -L_x; %negative distance not allowed
end
if L_x == 0
    L_x = 1;
end

set(hObject,'String',L_x)


handles.Lx = L_x;

dxsig = handles.dxsig*handles.km;
Lx = handles.Lx*handles.km;
xobs = handles.xobs;
xpadding = handles.paddingprismsx;
xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;

handles.totalprismsx = length(xpos);

% display the number of prisms in the textbox
set(handles.edit35, 'string', handles.totalprismsx);

% hide the invert button
set(handles.pushbuttonInvert, 'visible', 'off')
set(handles.text37, 'visible', 'on')
set(handles.pushbutton6, 'visible', 'off')
set(handles.text35, 'visible', 'on')
handles.dispbuttdens = 0;
handles.dispbuttdepth = 0;

guidata(hObject,handles)



function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.Lx = 1;
set(hObject,'String',handles.Lx)

guidata(hObject, handles)


% seths the y-length per prism
function edit3_Callback(hObject, eventdata, handles)

L_y = str2double(get(hObject,'String'));

if L_y < 0
    L_y = -L_y; %negative distance not allowed
end
if L_y == 0
    L_y = 1;
end
set(hObject,'String',L_y)

handles.Ly = L_y;

dysig = handles.dysig*handles.km;
Ly = handles.Ly*handles.km;
yobs = handles.yobs;
ypadding = handles.paddingprismsy;
ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

handles.totalprismsy = length(ypos);

% display the number of prisms in the textbox
set(handles.edit36, 'string', handles.totalprismsy);

% hide the invert button
set(handles.pushbuttonInvert, 'visible', 'off')
set(handles.text37, 'visible', 'on')
set(handles.pushbutton6, 'visible', 'off')
set(handles.text35, 'visible', 'on')
handles.dispbuttdens = 0;
handles.dispbuttdepth = 0;

guidata(hObject,handles)


function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.Ly = 1;
set(hObject,'String',handles.Ly)

guidata(hObject, handles)



% checkbox for enabaling dynamic dampening factor k
function checkbox5_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    handles.dynk = 1;
else
    handles.dynk = 0;
end
guidata(hObject, handles)


% defines the number of iterations for the inversion process
function edititer_Callback(hObject, eventdata, handles)

iter = str2double(get(hObject,'String'));


if iter <1
    iter = 1; %negative number of iterations not allowed
end
iter = round(iter);
set(hObject,'String',iter)

handles.n = iter;


%     calculation time estimations
% The computation time is obtained from the previous calculation 
% and then just multiplied by the number of iterations.
tmin = handles.tmin;
tsec = handles.tsec;
   
t = tmin*60 + tsec;

nt = handles.n*t;
ntm = floor(nt/60);
nts = nt - 60*ntm;

set(handles.textMinutes, 'String', ntm)
set(handles.textSeconds, 'String', nts)


guidata(hObject,handles)


function edititer_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.n = 1;
set(hObject,'String',handles.n)

handles.tocFm = 0; %for the timer calculations
guidata(hObject,handles)


% base number for k
function editk1_Callback(hObject, eventdata, handles)
k1 = str2double(get(hObject,'String'));
if k1 < 0
    k1 = 0; %negative dampening not allowed
end

handles.kay1 = k1;
set(hObject,'String',handles.kay1)
guidata(hObject, handles)


function editk1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.kay1 = 1;
set(hObject,'String',handles.kay1)
guidata(hObject,handles)


%  sets the x signal station spacing
function edit8_Callback(hObject, eventdata, handles)

dx_sig = str2double(get(hObject,'String'));


if dx_sig < 0
    dx_sig = -dx_sig; %negative distance not allowed
end
if dx_sig == 0
    dx_sig = 1;
end
if dx_sig > handles.Lx  %Prevents an underdetermined problem
    handles.Lx = dx_sig;
    set(handles.edit2, 'String', dx_sig)
end

set(hObject,'String',dx_sig)


handles.dxsig = dx_sig;

handles.x1split = handles.x1split*handles.dxsig;
handles.x1 = handles.x1*handles.dxsig;
handles.x0split = handles.x0split*handles.dxsig;
handles.x0 = handles.x0*handles.dxsig;

set(handles.edit13,'String',handles.x1split)
set(handles.editsplitx1,'String',handles.x1split)
set(handles.edit12,'String',handles.x0split)
set(handles.editsplitx0,'String',handles.x0split)


dxsig = handles.dxsig*handles.km;
Lx = handles.Lx*handles.km;
xobs = handles.xobs;
xpadding = handles.paddingprismsx;
xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;

handles.totalprismsx = length(xpos);

% display the number of prisms in the textbox
set(handles.edit35, 'string', handles.totalprismsx);

% hide the invert button
set(handles.pushbuttonInvert, 'visible', 'off')
set(handles.text37, 'visible', 'on')
set(handles.pushbutton6, 'visible', 'off')
set(handles.text35, 'visible', 'on')
handles.dispbuttdens = 0;
handles.dispbuttdepth = 0;

guidata(hObject,handles)


function edit8_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.dxsig = 1*handles.km;
set(hObject,'String',handles.dxsig)

guidata(hObject, handles)


% exponent number for k
function editk2_Callback(hObject, eventdata, handles)

k2 = str2double(get(hObject,'String'));
handles.kay2 = k2;
guidata(hObject, handles)



function editk2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.kay2 = -9;
strk2 = num2str(handles.kay2);
set(hObject,'String',strk2)
guidata(hObject, handles)


% checkbox for displaying observed Data Map and options
function checkbox1_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    set(handles.popupmenuObs, 'visible', 'on')
    set(handles.axesObs, 'visible', 'on')
    set(handles.pushbuttoncalcobs, 'visible', 'on')
else
    set(handles.popupmenuObs, 'visible', 'off')
    set(handles.axesObs, 'visible', 'off')
    set(handles.pushbuttoncalcobs, 'visible', 'off')
end
guidata(hObject, handles)


% Popup menu for chosing type of Observed data map plot
% 1D data create profile plots, 2D data have map plot options
function popupmenuObs_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
str = get(hObject, 'String');
switch str{val};
case 'Scaled Image'
    handles.Obsplottype = 1;
case 'Filled Contour (10)'
    handles.Obsplottype = 2;
case '3D Surface' 
    handles.Obsplottype = 3;
end

guidata(hObject,handles)



function popupmenuObs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% checkbox for displaying Forward model and options
function checkbox2_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    set(handles.popupmenuFmodel, 'visible', 'on')
    set(handles.axesFmodel, 'visible', 'on')
    set(handles.pushbuttoncalcfmodel, 'visible', 'on')
else
    set(handles.popupmenuFmodel, 'visible', 'off')
    set(handles.axesFmodel, 'visible', 'off')
    set(handles.pushbuttoncalcfmodel, 'visible', 'off')
end


% Popup menu for chosing type of forward model map types
function popupmenuFmodel_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
str = get(hObject, 'String');
switch str{val};
case 'Scaled Image'
	handles.Fmplottype = 1;
case 'Filled Contour (10)'
	handles.Fmplottype = 2;
case '3D Surface' 
	handles.Fmplottype = 3;
end

guidata(hObject,handles)



function popupmenuFmodel_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkboxBH_CreateFcn(hObject, eventdata, handles)

handles.Boreholeyes = 0;
guidata(hObject, handles);


% Plots the Observed Data as map or profile if Data only 1D
function pushbuttoncalcobs_Callback(hObject, eventdata, handles)

try

    obsmap = handles.obsmap;
    dxsig = handles.dxsig;
    dysig = handles.dysig;

    axes(handles.axesObs);  %#ok<*MAXES>
    hold off;
    xsmin=handles.x0split;
    xsmax=handles.x1split;
    ysmin=handles.y0split;
    ysmax=handles.y1split;

    [xobs, yobs] = size(obsmap);

    % Plots one of several types of plots
    % Converts map data into profile data if data are 1D

    if xobs > 1 && yobs > 1
        if handles.Obsplottype == 1
            imagesc(ysmin:dysig:ysmax,xsmin:dxsig:xsmax,obsmap); hold on; axis xy; axis equal; axis tight;
            colorbar('location','south');
        elseif handles.Obsplottype == 2
            contourf(obsmap, 10); axis xy; colorbar('location','south'); axis equal; axis tight; hold on;
        elseif handles.Obsplottype == 3
            surf(obsmap); axis xy; shading flat; colorbar('location','south');axis equal; axis tight; hold on;
        end
    else
        if xobs == 1
            plot(ysmin:dysig:ysmax, obsmap);
            dt1 = title('Observed Data profile');
            set(dt1, 'Position', [ysmax-0.5*(ysmax-ysmin),max(obsmap)-0.05*(max(obsmap)-min(obsmap)),0]);
        elseif yobs == 1
            plot(xsmin:dxsig:xsmax, obsmap);
            dt2 = title('Observed Data profile');
            set(dt2, 'Position', [xsmax-0.5*(xsmax-xsmin),max(obsmap)-0.05*(max(obsmap)-min(obsmap)),0]);
        end
    end
    
% Display error message instead of crashing
catch
    imperr = msgbox('Please enter an observed data set first', 'Import Error', 'warn');
    movegui(imperr, 'center')
end

guidata(hObject,handles)


% Plots the Forward Modelled as map or profile if Data only 1D
function pushbuttoncalcfmodel_Callback(hObject, eventdata, handles)

try

    dxsig = handles.dxsig;
    dysig = handles.dysig;
    xsmin=handles.x0split;
    xsmax=handles.x1split;
    ysmin=handles.y0split;
    ysmax=handles.y1split;

    axes(handles.axesFmodel); 

    calc0map = handles.c0;

    xobs = handles.xobs;
    yobs = handles.yobs;

    if xobs > 1 && yobs > 1
        if handles.Fmplottype == 1
            imagesc(ysmin:dysig:ysmax,xsmin:dxsig:xsmax,calc0map); axis xy; axis equal; axis tight;
            colorbar('location','south');
        elseif handles.Fmplottype == 2
            contourf(calc0map', 10); axis xy; colorbar('location','south');
        elseif handles.Fmplottype == 3
            surf(calc0map'); axis xy; shading flat; colorbar('location','south');
        end
    else
        if xobs == 1
            plot(ysmin:dysig:ysmax, calc0map);
            dt1 = title('Calculated Data profile');
            set(dt1, 'Position', [ysmax-0.5*(ysmax-ysmin),max(calc0map)-0.05*(max(calc0map)-min(calc0map)),0]);
        elseif yobs == 1
            plot(xsmin:dxsig:xsmax, calc0map);
            dt2 = title('Calculated Data profile');
            set(dt2, 'Position', [xsmax-0.5*(xsmax-xsmin),max(calc0map)-0.05*(max(calc0map)-min(calc0map)),0]);
        end
    end

% Display error message instead of crashing
catch
    imperr = msgbox('Please enter an observed data set first', 'Import Error', 'warn');
    movegui(imperr, 'center')
end

guidata(hObject,handles)



% displays Profile options and populates defaults
function checkbox11_CreateFcn(hObject, eventdata, handles)

handles.dispprof = 0;
handles.x0 = 0;
handles.x1 = 0;
handles.y0 = 0;
handles.y1 = 0;
guidata(hObject, handles)


% checkbox to display inverted depths in a table
function checkbox7_CreateFcn(hObject, eventdata, handles)

handles.dispobs = 0;
guidata(hObject, handles)


% checkbox to export inverted depths to excel sheet 
function checkbox8_CreateFcn(hObject, eventdata, handles)

handles.expobs = 0;
guidata(hObject, handles)


% toggles between km and m for all calculations
function radiobutton18_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    handles.km = 1;
else
    handles.km = 1/1000;
end
guidata(hObject,handles)


% toggles between km and m for all calculations
function radiobutton17_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    handles.km = 1/1000;
else
    handles.km = 1;
end
guidata(hObject,handles)



function radiobutton18_CreateFcn(hObject, eventdata, handles)

handles.km = 1;
guidata(hObject,handles)


% checkbox to display the misfit curve
function checkbox14_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    handles.msftdisp = 1;
else
    handles.msftdisp = 0;
end
guidata(hObject,handles)


function checkbox14_CreateFcn(hObject, eventdata, handles)

handles.msftdisp = 0;
guidata(hObject,handles)


function checkbox12_CreateFcn(hObject, eventdata, handles)

handles.dispmap = 0;
handles.azi = 290;
handles.elev = 35;
guidata(hObject, handles)



function pushbuttoncalcobs_CreateFcn(hObject, eventdata, handles)

handles.Obsplottype = 1;
handles.dxsig = 1*handles.km;
handles.dysig= 1*handles.km;
guidata(hObject, handles)



function checkbox1_CreateFcn(hObject, eventdata, handles)



function axesObs_CreateFcn(hObject, eventdata, handles)



function axesFmodel_CreateFcn(hObject, eventdata, handles)



% Does an initial calculation to create the Forward model map in the GUI
function pushbutton6_Callback(hObject, eventdata, handles)

tic

try
    
    waitFm = waitbar(0,'0%');

    dxsig = handles.dxsig *handles.km;
    dysig = handles.dysig *handles.km;
    Lx = handles.Lx *handles.km;
    Ly = handles.Ly *handles.km;

    xobs = handles.xobs;
    yobs = handles.yobs;


    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;

    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

    % creates permutations of x and y coordinates as x-y pairs.
    count = 0;
    permxy = zeros([length(xpos)*length(ypos) 2]);
    for valx = xpos
        for valy = ypos
             count = count + 1;
             permxy(count,1) = valx;
             permxy(count,2) = valy;
        end
    end

    % Initial depth and density information added to permxy matrix
    depthdist = handles.depthdist;
    depth = reshape(depthdist,(length(xpos)*length(ypos)),1);
    origdepth = depth;

    densitydist = handles.densitydist;
    density = reshape(densitydist,(length(xpos)*length(ypos)),1);


    [nrprisms,~] = size(permxy);

    calc0 = zeros((xobs*yobs),1);


    for pn = 1:nrprisms
        xloc = permxy(pn,1);
        yloc = permxy(pn,2);
        D0 = origdepth(pn);
        p = density(pn);

        c0 = Fmodel(xloc,yloc,D0,p,xobs,yobs,Lx,Ly,dxsig,dysig,handles.rh);  % original forward model
        c0rs = reshape(c0,(xobs*yobs),1);       %c0 reshaped

        calc0 = calc0 + c0rs;

        waitbar(pn/nrprisms,waitFm, [num2str(100*pn/nrprisms), '%']);
    end


    calc0map = reshape(calc0,xobs,yobs);

    handles.c0 = calc0map;
    
    handles.tocFm = (round(100*toc))/100;
    
    set(handles.text28, 'visible', 'off')
    set(handles.text29, 'visible', 'on')
    set(handles.text30, 'visible', 'on')
    set(handles.text30, 'String', num2str(handles.tocFm))
    
	% set all relevant buttons to visible
    set(handles.popupmenuFmodel, 'visible', 'on')
    set(handles.axesFmodel, 'visible', 'on')
    set(handles.pushbuttoncalcfmodel, 'visible', 'on')
    
%   calculation time estimations
    mapdelay = 0;
    if handles.dispmap == 1
        mapdelay = 0.3;
    end
    pdelay = 0;
    if handles.dispprof == 1
        pdelay = 0.25;
    end
    
%   The computation during inversion takes about 2.75 times longer than a
%    simple forward modelling. Found empirically
    tmin = floor((handles.n*(2.75*handles.tocFm + mapdelay + pdelay))/60);
    tsec = (handles.n*(2.75*handles.tocFm + mapdelay + pdelay)) - 60*floor((handles.n*(2.75*handles.tocFm + mapdelay + pdelay))/60);
    
    set(handles.textMinutes, 'String', num2str(tmin))
    set(handles.textSeconds, 'String', num2str(tsec))

    close(waitFm);

    catch  %#ok<*CTCH> % Display error message instead of crashing program
    toc
    moderr = msgbox(['Forward modelling was aborted after ', num2str(toc), ' seconds'], 'Operation cancelled' ,'warn');
    movegui(moderr, 'center')
    
end
    
guidata(hObject, handles)



function pushbuttoncalcfmodel_CreateFcn(hObject, eventdata, handles)

handles.Fmplottype = 1;
handles.dxsig = 1;
handles.dysig= 1;
handles.Lx = 1;
handles.Ly = 1;

guidata(hObject, handles)



function text28_CreateFcn(hObject, eventdata, handles)



function text30_CreateFcn(hObject, eventdata, handles)


% Padding checkbox.  
function checkbox15_Callback(hObject, eventdata, handles)

try

    isDown = get(hObject,'Value');
    if isDown
        set(handles.edit26, 'visible', 'on')
        set(handles.text33, 'visible', 'on')
        set(handles.edit31, 'visible', 'on')
        set(handles.text45, 'visible', 'on')
        handles.padding = 1;
    else
        set(handles.edit26, 'visible', 'off')
        set(handles.text33, 'visible', 'off')
        set(handles.edit31, 'visible', 'off')
        set(handles.text45, 'visible', 'off')
        handles.padding = 0;
        handles.paddingprismsx = 0;
        handles.paddingprismsy = 0;
        set(handles.edit31,'String','0')
        set(handles.edit26,'String','0')
        % hide the invert button
        set(handles.pushbuttonInvert, 'visible', 'off')
        set(handles.text37, 'visible', 'on')
        set(handles.pushbutton6, 'visible', 'off')
        set(handles.text35, 'visible', 'on')
        handles.dispbuttdens = 0;
        handles.dispbuttdepth = 0;
    end

    dxsig = handles.dxsig*handles.km;
    dysig = handles.dysig*handles.km;
    Lx = handles.Lx*handles.km;
    Ly = handles.Ly*handles.km;
    xobs = handles.xobs;
    yobs = handles.yobs;
    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;
    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

    handles.totalprismsx = length(xpos);
    handles.totalprismsy = length(ypos);

    % display the number of prisms in the textbox
    set(handles.edit35, 'string', handles.totalprismsx);
    set(handles.edit36, 'string', handles.totalprismsy);
    
catch
    imperr = msgbox('Please enter an observed dataset first', 'Import Error', 'warn');
    movegui(imperr, 'center')
end

guidata(hObject, handles)


% Padding text box.  States number of padding prisms in x direction
function edit26_Callback(hObject, eventdata, handles)

if handles.padding == 1;
    paddingprismsx = str2double(get(hObject,'String'));
    
    if paddingprismsx < 0
    paddingprismsx = 0; %negative number of prisms not allowed
    end
    paddingprismsx = round(paddingprismsx);
    set(hObject,'String',paddingprismsx)
    
    handles.paddingprismsx = paddingprismsx;
else
    handles.paddingprismsx = 0;

end

dxsig = handles.dxsig*handles.km;
Lx = handles.Lx*handles.km;
xobs = handles.xobs;
xpadding = handles.paddingprismsx;
xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;

handles.totalprismsx = length(xpos);

% display the number of prisms in the textbox
set(handles.edit35, 'string', handles.totalprismsx);

% hide the invert button
set(handles.pushbuttonInvert, 'visible', 'off')
set(handles.text37, 'visible', 'on')
set(handles.pushbutton6, 'visible', 'off')
set(handles.text35, 'visible', 'on')
handles.dispbuttdens = 0;
handles.dispbuttdepth = 0;

guidata(hObject, handles)



function edit26_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.paddingprismsx = 0;
strpaddingprismsx = num2str(handles.paddingprismsx);
set(hObject,'String',strpaddingprismsx)

guidata(hObject,handles)


% Padding text box.  States number of padding prisms in y direction
function edit31_Callback(hObject, eventdata, handles)
if handles.padding == 1;
    paddingprismsy = str2double(get(hObject,'String'));
    
    if paddingprismsy < 0
    paddingprismsy = 0; %negative number of prisms not allowed
    end
    paddingprismsy = round(paddingprismsy);
    set(hObject,'String',paddingprismsy)
    
    handles.paddingprismsy = paddingprismsy;
else
    handles.paddingprismsy = 0;
end

dysig = handles.dysig*handles.km;
Ly = handles.Ly*handles.km;
yobs = handles.yobs;
ypadding = handles.paddingprismsy;
ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

handles.totalprismsy = length(ypos);

% display the number of prisms in the textbox
set(handles.edit36, 'string', handles.totalprismsy);

% hide the invert button
set(handles.pushbuttonInvert, 'visible', 'off')
set(handles.text37, 'visible', 'on')
set(handles.pushbutton6, 'visible', 'off')
set(handles.text35, 'visible', 'on')
handles.dispbuttdens = 0;
handles.dispbuttdepth = 0;

guidata(hObject, handles)


function edit31_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.paddingprismsy = 0;
strpaddingprismsy = num2str(handles.paddingprismsy);
set(hObject,'String',strpaddingprismsy)

guidata(hObject,handles)




function edit27_Callback(hObject, eventdata, handles)



function edit27_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Populates the Depth table and sends the information to the GUI
function pushbuttonDepth_Callback(hObject, eventdata, handles)

try

    dxsig = handles.dxsig*handles.km;
    dysig = handles.dysig*handles.km;
    Lx = handles.Lx*handles.km;
    Ly = handles.Ly*handles.km;

    xobs = handles.xobs;
    yobs = handles.yobs;

    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;
    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

    handles.totalprismsx = length(xpos);
    handles.totalprismsy = length(ypos);

    depthdist(1:handles.totalprismsx, 1:handles.totalprismsy) = 1/handles.km;
    handles.depthdist = depthdist;

    % display the number of prisms in the textbox
    set(handles.edit35, 'string', handles.totalprismsx);
    set(handles.edit36, 'string', handles.totalprismsy);

    set(handles.table2, 'Data', handles.depthdist);


    estFmTime = (10^((0.625*log10(xobs*yobs*length(xpos)*length(ypos)))-2.2))/2.75;
    estFmTime3sf = (round(100*estFmTime)/100);

    tmin = floor((handles.n*2.75*estFmTime3sf)/60);
    tsec = (handles.n*2.75*estFmTime3sf) - 60*tmin;
    handles.tmin = tmin;
    handles.tsec = tsec;
    
    set(handles.textMinutes, 'String', num2str(tmin))
    set(handles.textSeconds, 'String', num2str(tsec))
    set(handles.text32, 'String', estFmTime3sf) 
    
    
    % Show the invert buttons

    handles.dispbuttdepth = 1;

    if handles.dispbuttdens ==1 && handles.dispbuttdepth ==1
        set(handles.pushbuttonInvert, 'visible', 'on')
        set(handles.text37, 'visible', 'off')
        set(handles.pushbutton6, 'visible', 'on')
        set(handles.text35, 'visible', 'off')
    end

% Display error message instead of crashing
catch
    imperr = msgbox('Please enter an observed data set first', 'Import Error', 'warn');
    movegui(imperr, 'center')
end

guidata(hObject, handles);


% Populates the density table and sends the information to the GUI
function pushbuttonDensity_Callback(hObject, eventdata, handles)

try

    dxsig = handles.dxsig*handles.km;
    dysig = handles.dysig*handles.km;
    Lx = handles.Lx*handles.km;
    Ly = handles.Ly*handles.km;

    xobs = handles.xobs;
    yobs = handles.yobs;

    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;
    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

    handles.totalprismsx = length(xpos);
    handles.totalprismsy = length(ypos);

    % create and fill the table with the right number of prisms
    densitydist(1:handles.totalprismsx, 1:handles.totalprismsy) = 1000;
    handles.densitydist = densitydist;

    % display the number of prisms in the textbox
    set(handles.edit35, 'string', handles.totalprismsx);
    set(handles.edit36, 'string', handles.totalprismsy);

    set(handles.table1, 'Data', handles.densitydist);

    estFmTime = (10^((0.625*log10(xobs*yobs*length(xpos)*length(ypos)))-2.2))/2.75;
    estFmTime3sf = (round(100*estFmTime)/100);

    tmin = floor((handles.n*2.75*estFmTime3sf)/60);
    tsec = (handles.n*2.75*estFmTime3sf) - 60*tmin;
    handles.tmin = tmin;
    handles.tsec = tsec;
    
    set(handles.textMinutes, 'String', num2str(tmin))
    set(handles.textSeconds, 'String', num2str(tsec))
    set(handles.text32, 'String', estFmTime3sf) 

    % Show the invert buttons

    handles.dispbuttdens = 1;

    if handles.dispbuttdens ==1 && handles.dispbuttdepth ==1
        set(handles.pushbuttonInvert, 'visible', 'on')
        set(handles.text37, 'visible', 'off')
        set(handles.pushbutton6, 'visible', 'on')
        set(handles.text35, 'visible', 'off')
    end
catch
    imperr = msgbox('Please enter an observed data set first', 'Import Error', 'warn');
    movegui(imperr, 'center')
end

guidata(hObject, handles);



function pushbuttonDepth_CreateFcn(hObject, eventdata, handles)
handles.dispbuttdepth = 0;
guidata(hObject, handles)

function pushbuttonInvert_CreateFcn(hObject, eventdata, handles)

function uipanel18_CreateFcn(hObject, eventdata, handles)

function radiobutton17_CreateFcn(hObject, eventdata, handles)

function uipanel3_CreateFcn(hObject, eventdata, handles)

function pushbutton6_CreateFcn(hObject, eventdata, handles)


% Creates depth table
function uitable2_CreateFcn(hObject, eventdata, handles)

handles.table2 = gcbo;
guidata(hObject, handles);



function uitable2_CellEditCallback(hObject, eventdata, handles)

NewData = get(hObject,'Data');
handles.depthdist = flipud(NewData);
guidata(hObject, handles);



% enables visibility of the depth custom constant value
function checkbox17_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    set(handles.edit29, 'visible', 'on')
else
    set(handles.edit29, 'visible', 'off')
end
guidata(hObject, handles)


% populates all depth values with a constant number
function edit29_Callback(hObject, eventdata, handles)

try

    dxsig = handles.dxsig*handles.km;
    dysig = handles.dysig*handles.km;
    Lx = handles.Lx*handles.km;
    Ly = handles.Ly*handles.km;

    xobs = handles.xobs;
    yobs = handles.yobs;

    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;
    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

    depthconst = str2double(get(hObject,'String'));

    if depthconst == 0
        depthconst = 0.00005;
    end

    handles.totalprismsx = length(xpos);
    handles.totalprismsy = length(ypos);

    depthdist(1:handles.totalprismsx, 1:handles.totalprismsy) = depthconst;
    handles.depthdist =depthdist;

    % display the number of prisms in the textbox
    set(handles.edit35, 'string', handles.totalprismsx);
    set(handles.edit36, 'string', handles.totalprismsy);

    estFmTime = (10^((0.625*log10(xobs*yobs*length(xpos)*length(ypos)))-2.2))/2.75;
    estFmTime3sf = (round(100*estFmTime)/100);

    tmin = floor((handles.n*2.75*estFmTime3sf)/60);
    tsec = (handles.n*2.75*estFmTime3sf) - 60*tmin;
    handles.tmin = tmin;
    handles.tsec = tsec;
    
    set(handles.textMinutes, 'String', num2str(tmin))
    set(handles.textSeconds, 'String', num2str(tsec))
    set(handles.text32, 'String', estFmTime3sf) 

    set(handles.table2, 'Data', handles.depthdist);

    % Show the invert buttons

    handles.dispbuttdepth = 1;

    if handles.dispbuttdens ==1 && handles.dispbuttdepth ==1
        set(handles.pushbuttonInvert, 'visible', 'on')
        set(handles.text37, 'visible', 'off')
        set(handles.pushbutton6, 'visible', 'on')
        set(handles.text35, 'visible', 'off')
    end
catch
    imperr = msgbox('Please enter an observed data set first', 'Import Error', 'warn');
    movegui(imperr, 'center')
end

guidata(hObject, handles);


function edit29_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% enables visibility of the density custom constant values
function checkbox16_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    set(handles.edit28, 'visible', 'on')
else
    set(handles.edit28, 'visible', 'off')
end
guidata(hObject, handles)


% populates all density values with a constant number
function edit28_Callback(hObject, eventdata, handles)

try

    dxsig = handles.dxsig*handles.km;
    dysig = handles.dysig*handles.km;
    Lx = handles.Lx*handles.km;
    Ly = handles.Ly*handles.km;

    xobs = handles.xobs;
    yobs = handles.yobs;

    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;
    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;

    handles.totalprismsx = length(xpos);
    handles.totalprismsy = length(ypos);

    densityconst = str2double(get(hObject,'String'));
    densitydist(1:handles.totalprismsx, 1:handles.totalprismsy) = densityconst;
    handles.densitydist = densitydist;

    % display the number of prisms in the textbox
    set(handles.edit35, 'string', handles.totalprismsx);
    set(handles.edit36, 'string', handles.totalprismsy);

    set(handles.table1, 'Data', handles.densitydist);


    estFmTime = (10^((0.625*log10(xobs*yobs*length(xpos)*length(ypos)))-2.2))/2.75;
    estFmTime3sf = (round(100*estFmTime)/100);

    tmin = floor((handles.n*2.75*estFmTime3sf)/60);
    tsec = (handles.n*2.75*estFmTime3sf) - 60*tmin;
    handles.tmin = tmin;
    handles.tsec = tsec;
    
    set(handles.textMinutes, 'String', num2str(tmin))
    set(handles.textSeconds, 'String', num2str(tsec))
    set(handles.text32, 'String', estFmTime3sf) 

    % Show the invert buttons

    handles.dispbuttdens = 1;

    if handles.dispbuttdens ==1 && handles.dispbuttdepth ==1
        set(handles.pushbuttonInvert, 'visible', 'on')
        set(handles.text37, 'visible', 'off')
        set(handles.pushbutton6, 'visible', 'on')
        set(handles.text35, 'visible', 'off')
    end

catch
    imperr = msgbox('Please enter an observed data set first', 'Import Error', 'warn');
    movegui(imperr, 'center')
end

guidata(hObject, handles);



function edit28_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_DeleteFcn(hObject, eventdata, handles)


function uitable1_CreateFcn(hObject, eventdata, handles)

handles.table1 = gcbo;
guidata(hObject, handles);



function uitable1_CellEditCallback(hObject, eventdata, handles)

NewData = get(hObject,'Data');
handles.densitydist = flipud(NewData);
guidata(hObject, handles);



function pushbuttonObs_CreateFcn(hObject, eventdata, handles)



function checkbox5_CreateFcn(hObject, eventdata, handles)

handles.dynk = 0;
guidata(hObject, handles);



function checkboxSplit_CreateFcn(hObject, eventdata, handles)

handles.split = 0;
guidata(hObject, handles)


% calculates and replaces the Observed dataset according to split borders
function pushbutton9_Callback(hObject, eventdata, handles)

try

    obsmap = (handles.obsmap);
    dxsig = handles.dxsig;
    dysig = handles.dysig;
    
    [xobs, yobs] = size(obsmap);

    xsmin=handles.x0split;
    xsmax=handles.x1split;
    ysmin=handles.y0split;
    ysmax=handles.y1split;

    
    if xsmax > dxsig*(xobs-1)/handles.km
        xsmax = dxsig*(xobs-1)/handles.km;
    end
    if ysmax > dysig*(yobs-1)/handles.km
        ysmax = dysig*(yobs-1)/handles.km;
    end
    
    xsplitmin = floor(xsmin/dxsig)+1;
    xsplitmax = ceil(xsmax/dxsig)+1;
    ysplitmin = floor(ysmin/dysig)+1;
    ysplitmax = ceil(ysmax/dysig)+1;

    try  %x and y swapped because axis is wrong way round
        obsmap = obsmap(xsplitmin:xsplitmax, ysplitmin:ysplitmax);
    catch %in case there is a matrix mismatch
        ploterr = msgbox('Please enter a value lower than the current limits', 'Dimension error', 'warn');
        movegui(ploterr, 'center')
    end
    
    [xobs, yobs] = size(obsmap);
    
    handles.obsmap = obsmap;
    handles.xobs = xobs;
    handles.yobs = yobs;

    if xobs > 1 && yobs > 1
        obs_list{1} = 'Scaled Image';
        obs_list{2} = 'Filled Contour (10)';
        obs_list{3} = '3D Surface';
        set(handles.popupmenuObs,'String',obs_list)
    else
        obs_list{1} = 'Profile';
        set(handles.popupmenuObs,'String',obs_list)
    end
    
    handles.x0split = 0;
    handles.x1split = xsmax-xsmin;
    handles.y0split = 0;
    handles.y1split = ysmax-ysmin;

    set(handles.editsplitx0,'String',handles.x0split);
    set(handles.editsplitx1,'String',handles.x1split);
    set(handles.editsplity0,'String',handles.y0split);
    set(handles.editsplity1,'String',handles.y1split);

    Lx = handles.Lx;
    Ly = handles.Ly;

    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;
    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;


    estFmTime = (10^((0.625*log10(xobs*yobs*length(xpos)*length(ypos)))-2.2))/2.75;
    estFmTime3sf = (round(100*estFmTime)/100);

    tmin = floor((handles.n*2.75*estFmTime3sf)/60);
    tsec = (handles.n*2.75*estFmTime3sf) - 60*tmin;
    handles.tmin = tmin;
    handles.tsec = tsec;
    
    set(handles.textMinutes, 'String', num2str(tmin))
    set(handles.textSeconds, 'String', num2str(tsec))
    set(handles.text32, 'String', estFmTime3sf) 

    length(xpos);
    length(ypos);

    handles.totalprismsx = length(xpos);
    handles.totalprismsy = length(ypos);

    % display the number of prisms in the textbox
    set(handles.edit35, 'string', handles.totalprismsx);
    set(handles.edit36, 'string', handles.totalprismsy);

    % hide the invert button
    set(handles.pushbuttonInvert, 'visible', 'off')
    set(handles.text37, 'visible', 'on')
    set(handles.pushbutton6, 'visible', 'off')
    set(handles.text35, 'visible', 'on')
    handles.dispbuttdens = 0;
    handles.dispbuttdepth = 0;
    
%   stops holding for plotted lines for cropping
    hold off;
    
catch
    dataerr = msgbox('Please enter an observed data set first', 'Import error', 'warn');
    movegui(dataerr, 'center')
end
    

guidata(hObject, handles)


% Checkbox to display Observed data in a seperate figure
function checkbox18_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    handles.expobsmap = 1;
else
    handles.expobsmap = 0;
end
guidata(hObject, handles)


% Checkbox to display Modelled data in a seperate figure
function checkbox19_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    handles.expmodmap = 1;
else
    handles.expmodmap = 0;
end
guidata(hObject, handles)



function checkbox18_CreateFcn(hObject, eventdata, handles)

handles.expobsmap = 0;
guidata(hObject, handles)



function checkbox19_CreateFcn(hObject, eventdata, handles)

handles.expmodmap = 0;
guidata(hObject, handles)


% Diplays options for Aliasing and editing observed data
function checkbox20_Callback(hObject, eventdata, handles)

isDown = get(hObject,'Value');
if isDown
    set(handles.pushbutton10, 'visible', 'on')
    set(handles.edit22, 'visible', 'on')
    set(handles.edit30, 'visible', 'on')
    set(handles.text44, 'visible', 'on')
    set(handles.text43, 'visible', 'on')
    set(handles.popupmenuObs, 'visible', 'on')
    set(handles.axesObs, 'visible', 'on')
    set(handles.pushbuttoncalcobs, 'visible', 'on')
else
    set(handles.pushbutton10, 'visible', 'off')
    set(handles.edit22, 'visible', 'off')
    set(handles.edit30, 'visible', 'off')
    set(handles.text44, 'visible', 'off')
    set(handles.text43, 'visible', 'off')
end
guidata(hObject, handles)


% Calculates and replaces the observed signal by a downsampled signal
function pushbutton10_Callback(hObject, eventdata, handles) %#ok<*DEFNU,*INUSL>

try

    % swopped because [y,x]=size(obsmap)
    aliasy = handles.aliasx;
    aliasx = handles.aliasy;

    obsmap = handles.obsmap;
    dxsig = handles.dxsig;
    dysig = handles.dysig;
    [xobs, yobs] = size(obsmap);

    dxsig = dxsig*aliasx;
    dysig = dysig*aliasy;
    handles.Lx = handles.Lx*aliasx;
    handles.Ly = handles.Ly*aliasy;
    set(handles.edit8,'String', dxsig)
    set(handles.edit1,'String', dysig)
    set(handles.edit2,'String', handles.Lx)
    set(handles.edit3,'String', handles.Ly)

    obsmap = obsmap(1:aliasx:xobs, 1:aliasy:yobs);

    [xobs, yobs] = size(obsmap);

    handles.obsmap = obsmap;
    handles.xobs = xobs;
    handles.yobs = yobs;
    handles.dxsig = dxsig;
    handles.dysig = dysig;

    if xobs > 1 && yobs > 1
        obs_list{1} = 'Scaled Image';
        obs_list{2} = 'Filled Contour (10)';
        obs_list{3} = '3D Surface';
        set(handles.popupmenuObs,'String',obs_list)
    else
        obs_list{1} = 'Profile';
        set(handles.popupmenuObs,'String',obs_list)
    end

    Lx = handles.Lx;
    Ly = handles.Ly;

    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;
    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;


    handles.totalprismsx = length(xpos);
    handles.totalprismsy = length(ypos);

    % redefine the downsized grid size
    x1split = handles.x1split;
    y1split = handles.y1split;
    
    if x1split > dxsig*(xobs-1)/handles.km
        x1split = dxsig*(xobs-1)/handles.km;
    end
    if y1split > dysig*(yobs-1)/handles.km
    
    y1split = dysig*(yobs-1)/handles.km;
    end
    
    set(handles.editsplity1,'String',y1split)
    set(handles.editsplitx1,'String',x1split)
    

    % display the number of prisms in the textbox
    set(handles.edit35, 'string', handles.totalprismsx);
    set(handles.edit36, 'string', handles.totalprismsy);

    estFmTime = (10^((0.625*log10(xobs*yobs*length(xpos)*length(ypos)))-2.2))/2.75;
    estFmTime3sf = (round(100*estFmTime)/100);

    tmin = floor((handles.n*2.75*estFmTime3sf)/60);
    tsec = (handles.n*2.75*estFmTime3sf) - 60*tmin;
    handles.tmin = tmin;
    handles.tsec = tsec;
    
    set(handles.textMinutes, 'String', num2str(tmin))
    set(handles.textSeconds, 'String', num2str(tsec))
    set(handles.text32, 'String', estFmTime3sf) 


    % hide the invert button
    set(handles.pushbuttonInvert, 'visible', 'off')
    set(handles.text37, 'visible', 'on')
    set(handles.pushbutton6, 'visible', 'off')
    set(handles.text35, 'visible', 'on')
    handles.dispbuttdens = 0;
    handles.dispbuttdepth = 0;
    
    handles.aliasx = 1;
    handles.aliasy = 1;
    
    set(handles.edit22, 'string', handles.aliasx);
    set(handles.edit30, 'string', handles.aliasy);
    
% Display error message instead of crashing
catch
    imperr = msgbox('Please enter an observed data set first', 'Import Error', 'warn');
    movegui(imperr, 'center')
end

guidata(hObject, handles)



function pushbutton11_Callback(hObject, eventdata, handles) %#ok<*INUSD>


% Aliasing factor in the y direction
function edit30_Callback(hObject, eventdata, handles)

aliasy = str2double(get(hObject,'String'));

if aliasy <1
    aliasy = 1; %negative number of iterations not allowed
end
aliasy = round(aliasy);
set(hObject,'String',aliasy)

handles.aliasy = aliasy;


guidata(hObject, handles)



function edit30_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.aliasy = 1;
set(hObject,'String',num2str(handles.aliasy))
guidata(hObject, handles)



function checkbox20_CreateFcn(hObject, eventdata, handles)


% imports excel sheet for initial depth values and updates station spacing
function pushbutton13_Callback(hObject, eventdata, handles)

try

    [ObsFileName,ObsPathName] = uigetfile('*.xlsx','Please select Initial Depth Distribution Data');
    IniDepthData = strcat(ObsPathName,ObsFileName);

    IniDepth = xlsread(IniDepthData);

    [xldepth,yldepth] = size(IniDepth);


    dxsig = handles.dxsig*handles.km;
    dysig = handles.dysig*handles.km;

    xobs = handles.xobs;
    yobs = handles.yobs;

    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;

    % create prism length from imported depth values
    Lx = ((xobs - 1)*(dxsig))/(xldepth - 2*xpadding - 1)/handles.km;
    Ly = ((yobs - 1)*(dysig))/(yldepth - 2*ypadding - 1)/handles.km;

    handles.Lx = Lx;
    set(handles.edit2,'String',handles.Lx)
    handles.Ly = Ly;
    set(handles.edit3,'String',handles.Ly)

    handles.totalprismsx = xldepth;
    handles.totalprismsy = yldepth;

    % display the number of prisms in the textbox
    set(handles.edit35, 'string', handles.totalprismsx);
    set(handles.edit36, 'string', handles.totalprismsy);

    handles.depthdist = IniDepth;

    set(handles.table2, 'Data', handles.depthdist);

    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;
    
    estFmTime = (10^((0.625*log10(xobs*yobs*length(xpos)*length(ypos)))-2.2))/2.75;
    estFmTime3sf = (round(100*estFmTime)/100);

    tmin = floor((handles.n*2.75*estFmTime3sf)/60);
    tsec = (handles.n*2.75*estFmTime3sf) - 60*tmin;
    handles.tmin = tmin;
    handles.tsec = tsec;
    
    set(handles.textMinutes, 'String', num2str(tmin))
    set(handles.textSeconds, 'String', num2str(tsec))
    set(handles.text32, 'String', estFmTime3sf) 

    % Show the invert buttons

    handles.dispbuttdepth = 1;

    if handles.dispbuttdens ==1 && handles.dispbuttdepth ==1
        set(handles.pushbuttonInvert, 'visible', 'on')
        set(handles.text37, 'visible', 'off')
        set(handles.pushbutton6, 'visible', 'on')
        set(handles.text35, 'visible', 'off')
    end
% Display error message instead of crashing
catch
    imperr = msgbox('Please enter an observed data set first or or enter a valid path for Depth values', 'Import Error', 'warn');
    movegui(imperr, 'center')
end

guidata(hObject, handles);


% imports excel sheet for initial density values and updates station spacin
function pushbutton12_Callback(hObject, eventdata, handles)

try

    [ObsFileName,ObsPathName] = uigetfile('*.xlsx','Please select Initial Density Distribution Data');
    IniDensData = strcat(ObsPathName,ObsFileName);

    IniDens = xlsread(IniDensData);

    [xldens,yldens] = size(IniDens);


    dxsig = handles.dxsig*handles.km;
    dysig = handles.dysig*handles.km;

    xobs = handles.xobs;
    yobs = handles.yobs;

    xpadding = handles.paddingprismsx;
    ypadding = handles.paddingprismsy;

    % create prism length from imported density values
    Lx = ((xobs - 1)*(dxsig))/(xldens - 2*xpadding - 1)/handles.km;
    Ly = ((yobs - 1)*(dysig))/(yldens - 2*ypadding - 1)/handles.km;

    handles.Lx = Lx;
    set(handles.edit2,'String',handles.Lx)
    handles.Ly = Ly;
    set(handles.edit3,'String',handles.Ly)

    handles.totalprismsx = xldens;
    handles.totalprismsy = yldens;

    % display the number of prisms in the textbox
    set(handles.edit35, 'string', handles.totalprismsx);
    set(handles.edit36, 'string', handles.totalprismsy);

    xpos = -xpadding*Lx:Lx:(xobs-1)*dxsig+xpadding*Lx;
    ypos = -ypadding*Ly:Ly:(yobs-1)*dysig+ypadding*Ly;
    
    handles.densitydist = IniDens;

    set(handles.table1, 'Data', handles.densitydist);

    estFmTime = (10^((0.625*log10(xobs*yobs*length(xpos)*length(ypos)))-2.2))/2.75;
    estFmTime3sf = (round(100*estFmTime)/100);

    tmin = floor((handles.n*2.75*estFmTime3sf)/60);
    tsec = (handles.n*2.75*estFmTime3sf) - 60*tmin;
    handles.tmin = tmin;
    handles.tsec = tsec;
    
    set(handles.textMinutes, 'String', num2str(tmin))
    set(handles.textSeconds, 'String', num2str(tsec))
    set(handles.text32, 'String', estFmTime3sf) 

    % Show the invert buttons

    handles.dispbuttdens = 1;

    if handles.dispbuttdens ==1 && handles.dispbuttdepth ==1
        set(handles.pushbuttonInvert, 'visible', 'on')
        set(handles.text37, 'visible', 'off')
        set(handles.pushbutton6, 'visible', 'on')
        set(handles.text35, 'visible', 'off')
    end
catch
    imperr = msgbox('Please enter an observed data set first or enter a valid path for Density values', 'Import Error', 'warn');
    movegui(imperr, 'center')
end

guidata(hObject, handles);


% Reference height
function edit32_Callback(hObject, eventdata, handles)
rh = str2double(get(hObject,'String'));
if rh <1
    rh = 1; %reference height of less than one meter not allowed
end
handles.rh = rh;
set(hObject,'String',handles.rh)
guidata(hObject,handles)


function edit32_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.rh = 1;
set(hObject,'String',handles.rh)

guidata(hObject, handles)

function checkbox15_CreateFcn(hObject, eventdata, handles)


% Exports the inverted data into an excel sheet
function checkbox22_Callback(hObject, eventdata, handles)
isDown = get(hObject,'Value');
if isDown
    handles.expinv = 1;
else
    handles.expinv = 0;
end
guidata(hObject, handles)

% checkbox to export inverted signal to excel sheet 
function checkbox22_CreateFcn(hObject, eventdata, handles)

handles.expinv = 0;
guidata(hObject, handles)



function checkbox23_Callback(hObject, eventdata, handles)
isDown = get(hObject,'Value');
if isDown
    set(handles.text47, 'visible', 'on')
    set(handles.edit34, 'visible', 'on')
    handles.dispsurfmap = 1;
else
    set(handles.text47, 'visible', 'off')
    set(handles.edit34, 'visible', 'off')
    handles.dispsurfmap = 0;
end

guidata(hObject, handles)



function checkbox23_CreateFcn(hObject, eventdata, handles)

handles.dispsurfmap = 0;
guidata(hObject, handles)



function edit34_Callback(hObject, eventdata, handles)
interpsurf = str2double(get(hObject,'String'));
handles.interpsurf = interpsurf;
guidata(hObject, handles)


function edit34_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.interpsurf = 0;
set(hObject,'String',handles.interpsurf)
guidata(hObject, handles)



function edit35_Callback(hObject, eventdata, handles)

totalprismsx = str2double(get(hObject,'String'));

if totalprismsx < 0
totalprismsx = 1; %negative number of prisms not allowed
end
totalprismsx = round(totalprismsx);
set(hObject,'String',totalprismsx)

handles.totalprismsx = totalprismsx;

% calculate the length of each prism
dxsig = handles.dxsig*handles.km;
xobs = handles.xobs;
xpadding = handles.paddingprismsx;

Lx = ((xobs - 1)*(dxsig))/(totalprismsx - 2*xpadding - 1)/handles.km;

handles.Lx = Lx;
set(handles.edit2,'String',handles.Lx)

% hide the invert button
set(handles.pushbuttonInvert, 'visible', 'off')
set(handles.text37, 'visible', 'on')
set(handles.pushbutton6, 'visible', 'off')
set(handles.text35, 'visible', 'on')
handles.dispbuttdens = 0;
handles.dispbuttdepth = 0;

guidata(hObject, handles)



function edit35_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.totalprismsx = 0;
strtotalprismsx = num2str(handles.totalprismsx);
set(hObject,'String',strtotalprismsx)

guidata(hObject,handles)


function edit36_Callback(hObject, eventdata, handles)
totalprismsy = str2double(get(hObject,'String'));

if totalprismsy < 0
totalprismsy = 1; %negative number of prisms not allowed
end
totalprismsy = round(totalprismsy);
set(hObject,'String',totalprismsy)

handles.totalprismsy = totalprismsy;

% calculate the length of each prism
dysig = handles.dysig*handles.km;
yobs = handles.yobs;
ypadding = handles.paddingprismsy;

Ly = ((yobs - 1)*(dysig))/(totalprismsy - 2*ypadding - 1)/handles.km;

handles.Ly = Ly;
set(handles.edit3,'String',handles.Ly)

% hide the invert button
set(handles.pushbuttonInvert, 'visible', 'off')
set(handles.text37, 'visible', 'on')
set(handles.pushbutton6, 'visible', 'off')
set(handles.text35, 'visible', 'on')
handles.dispbuttdens = 0;
handles.dispbuttdepth = 0;

guidata(hObject, handles)


function edit36_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.totalprismsy = 0;
strtotalprismsy = num2str(handles.totalprismsy);
set(hObject,'String',strtotalprismsy)

guidata(hObject,handles)



function pushbuttonDensity_CreateFcn(hObject, eventdata, handles)

handles.dispbuttdens = 0;
guidata(hObject,handles)


function pushbuttonInvert_ButtonDownFcn(hObject, eventdata, handles)

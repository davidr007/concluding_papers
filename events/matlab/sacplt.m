function varargout = sacplt(varargin)
%---------------------------------------------
%Coded by Erdinc Saygin-July 2004, RSES-ANU
%email:erdinc@rses.anu.edu.au
%---------------------------------------------


% SACPLT M-file for sacplt.fig
%      SACPLT, by itself, creates a new SACPLT or raises the existing
%      singleton*.
%
%      H = SACPLT returns the handle to a new SACPLT or the handle to
%      the existing singleton*.
%
%      SACPLT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SACPLT.M with the given input arguments.
%
%      SACPLT('Property','Value',...) creates a new SACPLT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sacplt_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sacplt_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sacplt

% Last Modified by GUIDE v2.5 04-Jul-2005 09:58:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sacplt_OpeningFcn, ...
                   'gui_OutputFcn',  @sacplt_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before sacplt is made visible.
function sacplt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sacplt (see VARARGIN)


% Choose default command line output for sacplt
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sacplt wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%reseting status to 0 at the beginning

handles.statetog1=0;
handles.pdfout=0;
handles.push1=0;
handles.push2=0;
handles.push3=0;
handles.sacwrite_st=0;
handles.sacbin1=0;
handles.sacbin2=0;
handles.sacbin3=0;
handles.mainpath=pwd;
handles.axislim=0;
handles.varpl=0;
handles.detrnd=0;
handles.markbegin=0;

guidata(hObject,handles);

% 
% f = uimenu('Label','File');
% uimenu(f,'Label','Load','Callback',uigetfile({'*.*', 'All Files (*.*)';
% '*.SAC;*.sac;*.mdl','SAC files (*.SAC,*.sac,*.*)'}, 'Choose a File'));
% uimenu(f,'Label','Write SAC file','Callback','save');



% --- Outputs from this function are returned to the command line.
function varargout = sacplt_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1

statetog=get(hObject,'Value');
handles.statetog1=statetog;
guidata(hObject,handles);

%--------------------------------------------------------------------------



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
push1=handles.push1;
push2=handles.push2;
push3=handles.push3;
sacwrite_st=handles.sacwrite_st;
filt1=handles.filt1;

if(filt1==1)
    w1=handles.w1;
    w2=handles.w2;
end


if (push1==1)
filename = handles.fil_1;
pathname = handles.path_1;
end

if (push2==1)
filename2 = handles.fil_2;
pathname2 = handles.path_2;
end

if (push3==1)
filename3 = handles.fil_3;
pathname3= handles.path_3;
end

pathnamecurr=pwd;





if (push1==1)
    sacbin1=handles.sacbin1;
    cd(pathname);
    
    if(sacbin1==1)
    typfil='b';
    else
    typfil='l';
    end
    
    [data1,kstnm1,delta,b,e,kcmnp,nzyear,nzjday,nzhour,nzmin,nzsec,stla,stlo,evla,evlo,cmpaz,cmpinc] = sacread(filename,typfil);
   
    
    handles.b=b;
    handles.delta=delta;
    handles.nzyear=nzyear;
    handles.nzjday=nzjday;
    handles.nzhour=nzhour;
    handles.nzmin=nzmin;
    handles.nzsec=nzsec;
    handles.kstnm1=kstnm1;
    handles.kcmnp=kcmnp;
    handles.stla=stla;
    handles.stlo=stlo;
    handles.evla=evla;
    handles.evlo=evlo;
    handles.cmpaz=cmpaz;
    handles.cmpinc=cmpinc;
    
    guidata(hObject,handles);
    
    cd(pathnamecurr);
end

if (push2==1)
    sacbin2=handles.sacbin2;
    cd(pathname2);
    
    if(sacbin2==1)
    typfil2='b';
    else
    typfil2='l';
    end
    
    [data2,kstnm2,delta,b2,e2,kcmnp2,nzyear2,nzjday2,nzhour2,nzmin2,nzsec2,stla2,stlo2,evla2,evlo2,cmpaz2,cmpinc2] = sacread(filename2,typfil2);
    
    handles.b2=b2;
    %handles.delta=delta;
    handles.nzyear2=nzyear2;
    handles.nzjday2=nzjday2;
    handles.nzhour2=nzhour2;
    handles.nzmin2=nzmin2;
    handles.nzsec2=nzsec2;
    handles.kstnm2=kstnm2;
    handles.kcmnp2=kcmnp2;
    handles.stla2=stla2;
    handles.stlo2=stlo2;
    handles.evla2=evla2;
    handles.evlo2=evlo2;
    handles.cmpaz2=cmpaz2;
    handles.cmpinc2=cmpinc2;
    guidata(hObject,handles);
    
    
    cd(pathnamecurr);
end

if (push3==1)
    sacbin3=handles.sacbin3;
    cd(pathname3);
    
    if(sacbin3==1)
    typfil3='b';
    else
    typfil3='l';
    end
    
    [data3,kstnm3,delta,b3,e3,kcmnp3,nzyear3,nzjday3,nzhour3,nzmin3,nzsec3,stla3,stlo3,evla3,evlo3,cmpaz3,cmpinc3] = sacread(filename3,typfil3);
    
    handles.b3=b3;
    %handles.delta=delta;
    handles.nzyear3=nzyear3;
    handles.nzjday3=nzjday3;
    handles.nzhour3=nzhour3;
    handles.nzmin3=nzmin3;
    handles.nzsec3=nzsec3;
    handles.kstnm3=kstnm3;
    handles.kcmnp3=kcmnp3;
    handles.stla3=stla3;
    handles.stlo3=stlo3;
    handles.evla3=evla3;
    handles.evlo3=evlo3;
    handles.cmpaz3=cmpaz3;
    handles.cmpinc3=cmpinc3;
    guidata(hObject,handles);
    
    cd(pathnamecurr);
end


statetog=handles.statetog1;

if (push1==1) %checking for the first graph if it is present
%checking if zoom button pressed and setting the values for ginput
if(statetog==1)
    [x,y]=ginput(2);
else
    x(1)=b;
    x(2)=e;
end


%saving x(1) and x(2)
handles.x1=x(1);
handles.x2=x(2);
guidata(hObject,handles);



%defining the indices (in points form, rounding the choosen points)
%if(statetog==1)
x(1);
x(2);
t0=fix((x(1)-b)./delta)+1; %make sure that indice starts from 1
t1=fix((x(2)-b)./delta)+1;
% else
% t0=((x(1)-b)./delta)+1    
% t1=((x(2)-b)./delta)+1
% end




%regenerating time axis

%t=b+(t0*delta):delta:b+(t1*delta); v0.1

t=b+((t0-1)*delta):delta:b+((t1-1)*delta); %v0.2
den1=b+((t0-1)*delta);
den2=b+((t1-1)*delta);




%This section writes the beginning (points) of the file and filename to
%separate file




end






%checking the axis lim
axslim=handles.axislim;
if(axslim==1)
%finding the maximum of y values for rearranging axis
if(push1==1)
    ymaxval(1)=max(data1);
    yminval(1)=min(data1);
else
   ymaxval(1)=0; 
   yminval(1)=0;
end

if(push2==1)
    ymaxval(2)=max(data2);
    yminval(2)=min(data2);
else
   ymaxval(2)=0;
   yminval(2)=0;
end

if(push3==1)
    ymaxval(3)=max(data3);
    yminval(3)=min(data3);
else
    ymaxval(3)=0;
    yminval(3)=0;
end
ymaxsing=max(ymaxval);
yminsing=min(yminval);

end
%--------------------------
varpl=handles.varpl; %getting the state of variable plot option
detrnd=handles.detrnd; %detrend option, obtraining value

if (push1==1)
    
%plotting with choosen time window
axes(handles.pl1)
    %filtering option
    if (filt1==1)
     w1n=2*delta*w1; %normalizing the frequency
     w2n=2*delta*w2;
     Wn=[w1n w2n];   
    [B,A] = BUTTER(2,Wn);
    %N1=fix((t1-t0+1)*0.1);
    
    %B=fir1(N1,Wn,hamming(N1+1));
    %A=1;
    data1_f=filtfilt(B,A,data1); 
    h=plot(t,data1_f(t0:t1));
    
    if(varpl==1)
    varplot_plain(data1_f(t0:t1),t(1),t(end),delta,0.03)
    end
    
%     if(axslim==1)
%     axis([t(1) t(end) (yminsing+0.1*yminsing) (ymaxsing+0.1*ymaxsing)])
%     end
    handles.data1f=data1_f;
    guidata(hObject,handles);
    %meandata=mean(data1_f(t0:t1));
    else
        if(detrnd==1)
            data1=detrend(data1,'constant');
        end
     

    h=plot(t,data1(t0:t1));

    
    if(varpl==1)
    varplot_plain(data1(t0:t1),t(1),t(end),delta,0.03)
    end

    if(axslim==1)
    axis([t(1) t(end) (yminsing+0.1*yminsing) (ymaxsing+0.1*ymaxsing)])
    end
    %h=stem(t,data1(t0:t1))
    %meandata=mean(data1(t0:t1));
    
    
    
    end


%labels and texts for plot
%if(t(end)>1000)
   %xlim([t(1)./3600 t(end)./3600])  
   
   xlabel('Time [hr]');
   %else
    
xlabel('Time [sec]');
ylabel('Amplitude');
%end


date1=strcat(num2str(nzjday),'-',num2str(nzyear),'/',num2str(nzhour),':',num2str(nzmin),':',num2str(nzsec));

t11=text(0,0,kstnm1,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom');
t12=text(0,1,kcmnp,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top');
t13=text(1,1,date1,'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
%t14=text(1,0,num2str(meandata),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom');
%t12=text(v(1),v(4),kcmnp,'HorizontalAlignment','left','VerticalAlignment','top');

set(t11,'Fontsize',10);
set(t12,'Fontsize',10,'FontWeight','bold');
set(t13,'Fontsize',10,'FontWeight','bold');
%set(t14,'Fontsize',10,'FontWeight','bold','Color','r');

%set(gcf,'PaperPositionMode','auto')
%set(gca,'PlotBoxAspectRatio',[2 0.5 1])


%get(t12)
grid on;
%     if(sacwrite_st==1)
%         filenamemod='modified.SAC';
%         sacmod(filename,filenamemod,'l',fix(x(1)),fix(x(2)));
%     end

%handles.data1=data1(t0:t1);
handles.t0=t0;
handles.t1=t1;
guidata(hObject,handles);



tete=(90-evla)*pi/180;
tets=(90-stla)*pi/180;

%---added 5 october 2004
if (evlo<0.0)
    evlo=360.0+evlo;
end

if (stlo<0.0)
    stlo=360.0+stlo;
end
%-----------------------------


fie=(evlo)*pi/180;
fis=(stlo)*pi/180;

delt=acos((cos(tete)*cos(tets))+(sin(tete)*sin(tets)*cos(fis-fie)));
delt=180*delt/pi;


set(handles.epcnt, 'String',num2str(delt));  


pdfout=handles.pdfout;
if(pdfout==1)
    %print -dpdf out.pdf
    print('-f1', '-ddf', 'out.pdf');
end





end

if (push2==1)
%PLOTTING THE SECOND GRAPH
axes(handles.pl2)
 %filtering option
    if (filt1==1)
     w1n=2*delta*w1; %normalizing the frequency
     w2n=2*delta*w2;
     Wn=[w1n w2n];   
    [B,A] = BUTTER(2,Wn);
    %N1=fix((t1-t0+1)*0.1);
    %[B,A]=fir1(N1,Wn,hamming(N1+1));
    data2_f=filtfilt(B,A,data2);
    handles.data2f=data2_f;
    guidata(hObject,handles);
    %meandata2=mean(data2_f(t0:t1));
    h2=plot(t,data2_f(t0:t1));
    
    if(varpl==1)
    varplot_plain(data2_f(t0:t1),t(1),t(end),delta,0.03)
    end
    
    else
    h2=plot(t,data2(t0:t1));
    if(axslim==1)
    axis([t(1) t(end) (yminsing+0.1*yminsing) (ymaxsing+0.1*ymaxsing)])
    end
    %meandata2=mean(data2(t0:t1));
    
    if(varpl==1)
    varplot_plain(data2(t0:t1),t(1),t(end),delta,0.03)
    end
    
    end



%h2=plot(t,data2(t0:t1));

xlabel('Time [sec]');
ylabel('Amplitude');
date2=strcat(num2str(nzjday2),'-',num2str(nzyear2),'/',num2str(nzhour2),':',num2str(nzmin2),':',num2str(nzsec2));
%v=axis; %values of axis

t21=text(0,0,kstnm2,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom');
t22=text(0,1,kcmnp2,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top');
t23=text(1,1,date2,'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
%t24=text(1,0,num2str(meandata2),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom');


set(t21,'Fontsize',10);
set(t22,'Fontsize',10,'FontWeight','bold');
set(t23,'Fontsize',10,'FontWeight','bold');
%set(t24,'Fontsize',10,'FontWeight','bold','Color','r');

grid on;

%handles.data2=data2(t0:t1);
handles.t0=t0;
handles.t1=t1;
guidata(hObject,handles);


end

if (push3==1)
%PLOTTING THE THIRD GRAPH
axes(handles.pl3);
    if (filt1==1)
     w1n=2*delta*w1; %normalizing the frequency
     w2n=2*delta*w2;
     Wn=[w1n w2n];   
    [B,A] = BUTTER(2,Wn);
    %N1=fix((t1-t0+1)*0.1);
    %[B,A]=fir1(N1,Wn,hamming(N1+1));
    data3_f=filtfilt(B,A,data3); 
    handles.data3f=data3_f;
    guidata(hObject,handles);
    %meandata3=mean(data3_f(t0:t1));
    h3=plot(t,data3_f(t0:t1));
    
    if(varpl==1)
    varplot_plain(data3_f(t0:t1),t(1),t(end),delta,0.03)
    end
    
    else
    h3=plot(t,data3(t0:t1));
    
    if(varpl==1)
    varplot_plain(data3(t0:t1),t(1),t(end),delta,0.03)
    end
    
    if(axslim==1)
    axis([t(1) t(end) (yminsing+0.1*yminsing) (ymaxsing+0.1*ymaxsing)])
    end
    %meandata3=mean(data3(t0:t1));
    end
    
%h3=plot(t,data3(t0:t1));
xlabel('Time [sec]');
ylabel('Amplitude');
%v=axis; %values of axis
date3=strcat(num2str(nzjday3),'-',num2str(nzyear3),'/',num2str(nzhour3),':',num2str(nzmin3),':',num2str(nzsec3));

t31=text(0,0,kstnm3,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom');
t32=text(0,1,kcmnp3,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top');
t33=text(1,1,date3,'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
%t34=text(1,0,num2str(meandata3),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom');

set(t31,'Fontsize',10);
set(t32,'Fontsize',10,'FontWeight','bold');
set(t33,'Fontsize',10,'FontWeight','bold');
%set(t34,'Fontsize',10,'FontWeight','bold','Color','r');
grid on;



%handles.data3=data3(t0:t1);
handles.t0=t0;
handles.t1=t1;
guidata(hObject,handles);

end



%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.*', 'All Files (*.*)';
'*.SAC;*.sac;*.mdl','SAC files (*.SAC,*.sac,*.*)'}, 'Choose a File');
if (filename~=0)
handles.fil_1 = filename;
handles.path_1= pathname;

push1=get(hObject,'Value');
handles.push1=push1;
guidata(hObject,handles);
else
push1=0
handles.push1=push1;
guidata(hObject,handles);  
end
    
% --- Executes on button press in pushl2.
function pushl2_Callback(hObject, eventdata, handles)
% hObject    handle to pushl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pathname = handles.path_1;

cd(pathname);

[filename2, pathname2] = uigetfile({'*.*', 'Alx Files (*.*)';
'*.SAC;*.sac;*.mdl','SAC files (*.SAC,*.sac,*.*)'}, 'Choose a File');

if (filename2~=0)
handles.fil_2 = filename2;
handles.path_2= pathname2;

push2=get(hObject,'Value');
handles.push2=push2;
guidata(hObject,handles);
else
push2=0;
handles.push2=push2;
guidata(hObject,handles);    
end

% --- Executes on button press in psuh_l3.
function psuh_l3_Callback(hObject, eventdata, handles)
% hObject    handle to psuh_l3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in checkbox1.
pathname = handles.path_1;

cd(pathname);

[filename3, pathname3] = uigetfile({'*.*', 'All Files (*.*)';
'*.SAC;*.sac;*.mdl','SAC files (*.SAC,*.sac,*.*)'}, 'Choose a File');

if (filename3~=0)
handles.fil_3 = filename3;
handles.path_3= pathname3;

push3=get(hObject,'Value');
handles.push3=push3;
guidata(hObject,handles);
else 
push3=0;
handles.push3=push3;
guidata(hObject,handles);
end


% --- Executes on button press in sacwrite.
function sacwrite_Callback(hObject, eventdata, handles)
% hObject    handle to sacwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sacwrite_st=get(hObject,'Value');
handles.sacwrite_st=sacwrite_st;
guidata(hObject,handles);

push1=handles.push1;
push2=handles.push2;
push3=handles.push3;

filt1=handles.filt1;



if(push1==1)
filename = handles.fil_1;
x(1)=handles.x1;
x(2)=handles.x2;

%data1=handles.data1;
pathname = handles.path_1;
t0=handles.t0
t1=handles.t1

%checking the file type according to check box
    sacbin1=handles.sacbin1;
    
    if(sacbin1==1)
    typfil='b';
    else
    typfil='l';
    end

cd(pathname);
%filename1_2=strcat(filename,'.',num2str(t0*0.04),'_',num2str(t1*0.04));%v0.1
    if(filt1==1) %checking for filter option
    data1f=handles.data1f;    
    delta=handles.delta;
    nzyear=handles.nzyear;
    nzjday=handles.nzjday;
    nzhour=handles.nzhour;
    nzmin=handles.nzmin;
    nzsec=handles.nzsec;
    kstnm1=handles.kstnm1;
    kcmnp=handles.kcmnp;
    b=handles.b;
    stla=handles.stla;
    stlo=handles.stlo;
    evla=handles.evla;
    evlo=handles.evlo;
    cmpaz=handles.cmpaz;
    cmpinc=handles.cmpinc;
    
    bf=b+((t0-1).*delta);
    ef=b+((t1-1).*delta);
    
    filename1_2=strcat(filename,'.',num2str((t0-1)*delta),'_',num2str((t1-1)*delta),'F');
    sacwrite(filename1_2,typfil,data1f(t0:t1),delta,(t1-t0)+1,bf,ef,nzyear,nzjday,nzhour,nzmin,nzsec,kstnm1,kcmnp,stla,stlo,evla,evlo,cmpaz,cmpinc);
    
    else
    delta=handles.delta;    
    filename1_2=strcat(filename,'.',num2str((t0-1)*delta),'_',num2str((t1-1)*delta));
    sacmodpo(filename,filename1_2,typfil,t0,t1);
    end
end

if(push2==1)
filename2 = handles.fil_2;
x(1)=handles.x1;
x(2)=handles.x2;

%data2=handles.data2;
pathname2 = handles.path_2;
t0=handles.t0;
t1=handles.t1;


%checking the file type according to check box
    sacbin2=handles.sacbin2;
    
    if(sacbin2==1)
    typfil2='b';
    else
    typfil2='l';
    end
cd(pathname2);
%filename2_2=strcat(filename2,'.',num2str(t0*0.04),'_',num2str(t1*0.04)); %v0.1
    if(filt1==1) %checking for filter option
    data2f=handles.data2f;    
    delta=handles.delta;
    nzyear2=handles.nzyear2;
    nzjday2=handles.nzjday2;
    nzhour2=handles.nzhour2;
    nzmin2=handles.nzmin2;
    nzsec2=handles.nzsec2;
    kstnm2=handles.kstnm2;
    kcmnp2=handles.kcmnp2;
    b2=handles.b2;
    stla2=handles.stla2
    stlo2=handles.stlo2;
    evla2=handles.evla2;
    evlo2=handles.evlo2;
    cmpaz2=handles.cmpaz2;
    cmpinc2=handles.cmpinc2;
    
    bf2=b2+((t0-1).*delta);
    ef2=b2+((t1-1).*delta);
    
    filename2_2=strcat(filename2,'.',num2str((t0-1)*delta),'_',num2str((t1-1)*delta),'F');
    sacwrite(filename2_2,typfil,data2f(t0:t1),delta,(t1-t0)+1,bf2,ef2,nzyear2,nzjday2,nzhour2,nzmin2,nzsec2,kstnm2,kcmnp2,stla2,stlo2,evla2,evlo2,cmpaz2,cmpinc2);
    
    else
    delta=handles.delta;    
    filename2_2=strcat(filename2,'.',num2str((t0-1)*delta),'_',num2str((t1-1)*delta));
    sacmodpo(filename2,filename2_2,typfil,t0,t1);
    end
% filename2_2=strcat(filename2,'.',num2str((t0-1)*0.04),'_',num2str((t1-1)*0.04));  %v0.2
% 
% sacmodpo(filename2,filename2_2,typfil2,t0,t1);
end

if(push3==1)
filename3 = handles.fil_3;
x(1)=handles.x1;
x(2)=handles.x2;

%data3=handles.data3;
pathname3= handles.path_3;
t0=handles.t0;
t1=handles.t1;

%checking the file type according to check box
    sacbin3=handles.sacbin3;
    
    if(sacbin3==1)
    typfil3='b';
    else
    typfil3='l';
    end
cd(pathname3);
% %filename3_2=strcat(filename3,'.',num2str(t0*0.04),'_',num2str(t1*0.04)); %v0.1
% filename3_2=strcat(filename3,'.',num2str((t0-1)*0.04),'_',num2str((t1-1)*0.04)); %v0.2
% sacmodpo(filename3,filename3_2,typfil3,t0,t1);

if(filt1==1) %checking for filter option
    data3f=handles.data3f;    
    delta=handles.delta;
    nzyear3=handles.nzyear3;
    nzjday3=handles.nzjday3;
    nzhour3=handles.nzhour3;
    nzmin3=handles.nzmin3;
    nzsec3=handles.nzsec3;
    kstnm3=handles.kstnm3;
    kcmnp3=handles.kcmnp3;
    b3=handles.b3;
    
    stla3=handles.stla3;
    stlo3=handles.stlo3;
    evla3=handles.evla3;
    evlo3=handles.evlo3;
    cmpaz3=handles.cmpaz3;
    cmpinc3=handles.cmpinc3;
    
    bf3=b3+((t0-1).*delta);
    ef3=b3+((t1-1).*delta);
    
    filename3_2=strcat(filename3,'.',num2str((t0-1)*delta),'_',num2str((t1-1)*delta),'F');
    sacwrite(filename3_2,typfil,data3f(t0:t1),delta,(t1-t0)+1,bf3,ef3,nzyear3,nzjday3,nzhour3,nzmin3,nzsec3,kstnm3,kcmnp3,stla3,stlo3,evla3,evlo3,cmpaz3,cmpinc3);
    
    else
    delta=handles.delta;    
    filename3_2=strcat(filename3,'.',num2str((t0-1)*delta),'_',num2str((t1-1)*delta));
    sacmodpo(filename3,filename3_2,typfil,t0,t1);
    end

    
end
mainpath=handles.mainpath;
cd(mainpath);

% --- Executes on button press in sacbin1.
function sacbin1_Callback(hObject, eventdata, handles)
% hObject    handle to sacbin1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sacbin1
sacbin1=get(hObject,'Value');
handles.sacbin1=sacbin1;
guidata(hObject,handles);

% --- Executes on button press in sacbin2.
function sacbin2_Callback(hObject, eventdata, handles)
% hObject    handle to sacbin2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sacbin2
sacbin2=get(hObject,'Value');
handles.sacbin2=sacbin2;
guidata(hObject,handles);

% --- Executes on button press in sacbin3.
function sacbin3_Callback(hObject, eventdata, handles)
% hObject    handle to sacbin3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sacbin3

sacbin3=get(hObject,'Value');
handles.sacbin3=sacbin3;
guidata(hObject,handles);







% --- Executes on button press in pushexit.
function pushexit_Callback(hObject, eventdata, handles)
% hObject    handle to pushexit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;






% --- Executes during object creation, after setting all properties.
function w1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
%handles.w1=w1
%guidata(hObject,handles);



function w1_Callback(hObject, eventdata, handles)
% hObject    handle to w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w1 as text
%        str2double(get(hObject,'String')) returns contents of w1 as a double
filt1=get(hObject,'Value')

w1=str2double(get(hObject,'String'))
if isnan(w1)
errordlg('You must enter a numeric value','Bad Input','modal')
%w1=1;
end

handles.w1=w1;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function w2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function w2_Callback(hObject, eventdata, handles)
% hObject    handle to w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w2 as text
%        str2double(get(hObject,'String')) returns contents of w2 as a double
w2=str2double(get(hObject,'String'))
if isnan(w2)
errordlg('You must enter a numeric value','Bad Input','modal')
%w2=2;
end

handles.w2=w2;
guidata(hObject,handles);



% --- Executes on button press in filt1.
function filt1_Callback(hObject, eventdata, handles)
% hObject    handle to filt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filt1

filt1=get(hObject,'Value');
handles.filt1=filt1;
guidata(hObject,handles);


%----Notes
%For every modified file, it may show one point less than the original. This is due to
%"b" and "e" values and fixing the division. This small obstacle won't hurt anything. 



% --- Executes on button press in pdf.
function pdf_Callback(hObject, eventdata, handles)
% hObject    handle to pdf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pdfout=get(hObject,'Value');
handles.pdfout=pdfout;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function epcnt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epcnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function epcnt_Callback(hObject, eventdata, handles)
% hObject    handle to epcnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epcnt as text
%        str2double(get(hObject,'String')) returns contents of epcnt as a double

set(handles.epcnt,'Enable','off');



% --- Executes on button press in axislim.
function axislim_Callback(hObject, eventdata, handles)
% hObject    handle to axislim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of axislim


axislim=get(hObject,'Value');
handles.axislim=axislim;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in varplot.
function varplot_Callback(hObject, eventdata, handles)
% hObject    handle to varplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of varplot

varpl=get(hObject,'Value');
handles.varpl=varpl;
guidata(hObject,handles);


% --- Executes on button press in detrnd.
function detrnd_Callback(hObject, eventdata, handles)
% hObject    handle to detrnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detrnd

detrnd=get(hObject,'Value');
handles.detrnd=detrnd;
guidata(hObject,handles);


% --- Executes on button press in markbegin.
function markbegin_Callback(hObject, eventdata, handles)
% hObject    handle to markbegin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of markbegin
%this function marks the beginning portion of the viewed file in points
%(not seconds)

markbegin=get(hObject,'Value');
handles.markbegin=markbegin;
guidata(hObject,handles);

fid2 = fopen('sactimemarked.dat','a'); %creation of file the append state
filenam=handles.fil_1;
t0=handles.t0;
fprintf(fid2,'%s %10d\n',filenam,t0);
fclose(fid2);


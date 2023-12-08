function varargout = my_mice_tracking_56_generalSize(varargin)
% MY_MICE_TRACKING_56_GENERALSIZE MATLAB code for my_mice_tracking_56_generalSize.fig
%      MY_MICE_TRACKING_56_GENERALSIZE, by itself, creates a new MY_MICE_TRACKING_56_GENERALSIZE or raises the existing
%      singleton*.
%
%      H = MY_MICE_TR   ACKING_56_GENERALSIZE returns the handle to a new MY_MICE_TRACKING_56_GENERALSIZE or the handle to
%      the existing singleton*.
%
%      MY_MICE_TRACKING_56_GENERALSIZE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MY_MICE_TRACKING_56_GENERALSIZE.M with the given input arguments.
%C:\Users\zrpgs\Documents\Adobe\Premiere Pro\14.0
%      MY_MICE_TRACKING_56_GENERALSIZE('Property','Value',...) creates a new MY_MICE_TRACKING_56_GENERALSIZE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before my_mice_tracking_56_generalSize_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to my_mice_tracking_56_generalSize_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% %  based on my_mice_tracking_53_TimeMap_NewFilter, mask part from
% my_mice_tracking_4a_mask_Xor_v2 
%
% Edit the above text to modify the response to help my_mice_tracking_56_generalSize

% Last Modified by GUIDE v2.5 13-Jan-2017 17:16:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @my_mice_tracking_56_generalSize_OpeningFcn, ...
    'gui_OutputFcn',  @my_mice_tracking_56_generalSize_OutputFcn, ...
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


% --- Executes just before my_mice_tracking_56_generalSize is made visible.
function my_mice_tracking_56_generalSize_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to my_mice_tracking_56_generalSize (see VARARGIN)

% Choose default command line output for my_mice_tracking_56_generalSize
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
axes(handles.axes1);

text(0.15,0.55,'Mice Tracking','FontSize',30);
set(gca, 'XTick', []);
set(gca, 'YTick', []);

% axes(handles.axes2);
% text(0.15,0.55,'Raw Signal','FontSize',10);
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);

axes(handles.axes3);
text(0.15,0.55,'Baseline Correction Signal','FontSize',10);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
% UIWAIT makes my_mice_tracking_56_generalSize wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = my_mice_tracking_56_generalSize_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonLoadAVI.
function pushbuttonLoadAVI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadAVI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);
axes(handles.axes1);

%load an AVI


set(handles.text13,'String','loading ...');
[filename, pathname, filterindex] = uigetfile('*.avi', 'Select an AVI file.');
if filterindex > 0
    
    VideoObj = VideoReader([pathname,filename]);
    set(handles.editVideoPath,'String', [pathname filename]);
    set(handles.text13,'String','file loaded');
    % Extract video information
    handles.VideoObj=VideoObj;
    % Initialize the displaying size of video
    handles.Control.Xrange=[1, VideoObj.Width];
    handles.Control.Yrange=[1, VideoObj.Height];
    % Display the first frame
    ImageFrame = read(VideoObj, 1);
    imagesc(ImageFrame);
    axis image;
    axis off;
    
    axes(handles.axes3); cla;
    
    NumberofFrames = handles.VideoObj.NumberofFrames;
    % Initialize the sliding control of the video
    SliderStep(1) = 1/(NumberofFrames-1);
    SliderStep(2) = 2/(NumberofFrames-1);
    set(handles.ImageIndexSlider,'sliderstep',SliderStep,'max',NumberofFrames,'min',1,'Value',1);
    set(handles.FrameIndex,'String', '1');
    
    handles.pathname=pathname;
    handles.filename=filename;
    
    tempoffset = 0;% default is 5
    handles.Xrange = [tempoffset+1, VideoObj.Width-tempoffset];
    handles.Yrange = [tempoffset+1, VideoObj.Height-tempoffset];
    
    set(handles.StartTimeFrame,'String',1);
    set(handles.EndTimeFrame,'String',NumberofFrames);
    
    handles.tempmask =0;
    handles.trackOk = false;
    set(handles.edit_lengthofpath,'String','');
end
set(handles.text13,'String',' ');
guidata(hObject,handles);
axes(handles.axes1);



% --- Executes on button press in pushbuttonPlayAVI.
function pushbuttonPlayAVI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlayAVI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1);
cla;

% Check whether the Slider button is pushed
SliderButtonValue=get(handles.ImageIndexSlider,'Value');
FrameRate = handles.VideoObj.FrameRate;
NumberOfFrames=handles.VideoObj.NumberOfFrames;
% Initialize the sliding control of the video
SliderStep(1) = 1/(NumberOfFrames-1);
SliderStep(2) = 2/(NumberOfFrames-1);
set(handles.ImageIndexSlider,'sliderstep',SliderStep,'max',NumberOfFrames,'min',1,'Value',SliderButtonValue);

%check the tracking file on the current fold
if ~handles.trackOk
    trackfilename = [handles.pathname, handles.filename(1:end-4) '_track.mat'];
    Xoffset = 1; Yoffset = 1;
    if exist(trackfilename,'file')==2
        load(trackfilename);
        handles.xypos = xypos ;
        handles.trackOk = true;
        
        TempMaskfile = [handles.pathname, handles.filename(1:end-4) '_TempMask.mat'];
        if exist(TempMaskfile,'file')==2
            load(TempMaskfile);
            handles.tempmask = tempmask ;
        end
        
        if  exist('threshold','var')
            handles.threshold  = threshold ;
        end
        if  exist('TestRadium','var')
            handles.TestRadium = TestRadium;
        end
        
        
        %     else
        %          Xoffset = 0;
        %          Yoffset = 0;
    end
else
    Xoffset = handles.Xrange(1);
    Yoffset = handles.Yrange(1);
    xypos = handles.xypos;
end
Xoffset = Xoffset- 1;
Yoffset = Yoffset -1;
%
% prepare  mark circle
ang=0:0.01:2*pi;
r=10;
xp=r*cos(ang)+Xoffset; yp=r*sin(ang)+Yoffset;
%     pb=patch(xydata(1)+xp,xydata(2)+yp,'r','edgecolor','r');
%     alpha(pb,0.5);


%Displaying each frame of the video
for ImageIndex=SliderButtonValue:2:NumberOfFrames
    %   SliderButtonValue=handles.SliderValue;
    ImageDataFrame = read(handles.VideoObj, ImageIndex);
    imagesc(ImageDataFrame);
    
    axis image;
    axis off;
    
    if handles.trackOk
        x0=xypos(1,ImageIndex);
        y0=xypos(2,ImageIndex);
        pb=patch(x0+xp,y0+yp,'g','edgecolor','g');
        alpha(pb,0.4);
    end
    %caculate the play time
    Time=roundn(ImageIndex/FrameRate,-3);
    StringofTime=num2str(Time);
    set(handles.editTime,'String',StringofTime);
    % Set the frame index and Image Index Slider
    set(handles.ImageIndexSlider,'Value',ImageIndex);
    set(handles.FrameIndex,'String',num2str(ImageIndex));
    % Check whether the pause button is pushed
    PauseButtonValue=get(handles.pushbuttonPause ,'Value');
    if PauseButtonValue
        text = get(handles.pushbuttonPause, 'String');
        if strcmp(text, 'Pause') == 1
            set(handles.pushbuttonPause,'String','Continue');
            uiwait;
            %waitforbuttonpress;
            %set(handles.pushbuttonPause,'Value',0);
            %set(handles.pushbuttonPause,'String','Pause');
        end
    end
    
    % Check whether the stop button is pushed
    StopButtonValue=get(handles.pushbuttonStop,'Value');
    if StopButtonValue
        ImageDataFrame = read(handles.VideoObj, ImageIndex);
        imagesc(ImageDataFrame);
        axis  image off;
        set(handles.ImageIndexSlider,'Value',ImageIndex);
        set(handles.FrameIndex,'String',num2str(ImageIndex));
        set(handles.pushbuttonStop,'Value',0);
        set(handles.editTime,'String',StringofTime);
        % Quit the displaying
        break;
    end
    pause(0.01);
end

guidata(hObject,handles);
axes(handles.axes1);

% function pushbuttonCreatObjects_Callback
% guidata(hObject,handles);
% axes(handles.axes1);
%
% %creat the PrtScn range
% [Xpositions,Ypositions] = ginput(4);
% XpositionsSorted=sort(Xpositions);
% YpositionsSorted=sort(Ypositions);
% Xrange=ceil([(XpositionsSorted(2)), XpositionsSorted(3)]);
% Yrange=ceil([YpositionsSorted(2), YpositionsSorted(3)]);
% handles.Xrange=Xrange;
% handles.Yrange=Yrange;

% --- Executes on button press in pushbuttonCreatObjects.
function pushbuttonCreatObjects_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCreatObjects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);
axes(handles.axes1);

%creat the PrtScn range
[Xpositions,Ypositions] = ginput(4);
XpositionsSorted=sort(Xpositions);
YpositionsSorted=sort(Ypositions);
Xrange=ceil([(XpositionsSorted(2)), XpositionsSorted(3)]);
Yrange=ceil([YpositionsSorted(2), YpositionsSorted(3)]);
handles.Xrange=Xrange;
handles.Yrange=Yrange;

guidata(hObject,handles);

% --- Executes on button press in pushbuttonMiceTracking.
function pushbuttonMiceTracking_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMiceTracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1); colormap(gray);
% cla;
nFrames = handles.VideoObj.NumberOfFrames;
FrameRate = handles.VideoObj.FrameRate;
Yrange=handles.Yrange;
Xrange=handles.Xrange;
yst=Yrange(1):Yrange(2);
xst=Xrange(1):Xrange(2);
% Extract the information for the screen
StratFrame =  str2double(get(handles.StartTimeFrame,'String'));
EndFrame= str2double(get(handles.EndTimeFrame,'String'));

if EndFrame >= nFrames
    EndFrame = nFrames;
end


% prepare the working enviroment
clear ContourStoreSt;
ContourStoreSt{1}=[];
xypos = zeros(2,nFrames);    %xypos(:,k) = xydata;
se1=strel('disk',6);
ang=0:0.01:2*pi;    r=10;  xp=r*cos(ang);yp=r*sin(ang) ;
threshold = 0.90;
Error_r = str2double(get(handles.edit_r,'String')) ;

%   setup initial mice location manually
InitalFrame = round( get(handles.ImageIndexSlider,'Value'));
if InitalFrame >= EndFrame
    set(handles.text13,'String','Initial frame out of range');
    return;
end
tempmask = handles.tempmask;
if size(tempmask) < 100
    ht = msgbox('No Mask data available, stopped here','CreateMode','modal');
    return;
end
% % %     %     go to prepare the mask frame
% % %     set(handles.text13,'String','setup tempMask');
% % %     NumberofFrames = EndFrame - InitalFrame ;
% % %     orderRandom = randperm(NumberofFrames)+InitalFrame;
% % %     colormap(gray);
% % %     for k = 1:NumberofFrames
% % %         i = orderRandom(k);
% % %         ImageFrameIndex = read(handles.VideoObj, i);
% % %         I1=rgb2gray(ImageFrameIndex);
% % %         axes(handles.axes5); imshow(I1);
% % %         I2 = uint8( tempmask > I1) ;
% % %         tempmask = tempmask.*I2 + I1.*(1-I2);
% % %         axes(handles.axes1);imagesc(tempmask);
% % %         axis image off;
% % %         set(handles.ImageIndexSlider,'Value',i);
% % %         set(handles.FrameIndex,'String', i);
% % %         drawnow;
% % %         pause(3);
% % %         
% % %         StopButtonValue=get(handles.pushbuttonStop,'Value');
% % %         if StopButtonValue
% % %             set(handles.pushbuttonStop,'Value',0);
% % %             break; % Quit the displaying
% % %         end
% % %     end
% % %     handles.tempmask = tempmask;
% % %     save ([handles.pathname, handles.filename(1:end-4) '_TempMask'], 'tempmask');
    
    % % restore the initial mask;
%     set(handles.ImageIndexSlider,'Value',InitalFrame);
% end
I1=rgb2gray(read(handles.VideoObj, InitalFrame));
I2 = imsubtract(tempmask,I1);
imagesc(I2);         axis image off;

% % % % %
if handles.trackOk
    if isfield(handles,'xypos')
        if length(handles.xypos) == EndFrame
            xypos = handles.xypos ;
        end
    end
    if isfield(handles,'TestRadium')
        TestRadium = handles.TestRadium ;
        threshold = handles.threshold;
    end
    %      startframe = InitalFrame + 1;
    xydata = xypos(:,InitalFrame);
else
    %     I1 = rgb2gray( read(handles.VideoObj, 1) ); %intial the mask with 1st frame
    %     imagesc(I1); colormap(gray);     axis off;
    ht = msgbox('Select a rectangle area to hold the mouse','CreateMode','modal');
    axes(handles.axes1);
    ht = imrect;
    setColor(ht,'red');
    pos = getPosition(ht);
    TestRadium = max(pos(3),pos(4));
    
    
    %     set(handles.text13,'String','select mice location');
    %     [x0,y0] = ginput(1);
    x0 = pos(1) + pos(3)/2;
    y0 = pos(2) + pos(4)/2;
    
    xydata = round([x0;y0]);
    xypos(:,InitalFrame)= xydata;
    
    I2 = imcrop(I2, pos);
    set(handles.text13,'String','setup polygon');
    ht = msgbox('setup polygon to hold the mouse','CreateMode','modal');
    axes(handles.axes5); colormap(gray);
    imagesc(I2); axis image off;
    
    A1 = uint8( roipoly() );
    threshold = 0.9;
    TestRadium = 4* round(sqrt(sum(sum(A1)/pi*2)));
    
end
startframe = InitalFrame + 1;
TrackCompleted = true;
statusString = 'finished';
% % % % %
set(handles.text11,'String', num2str(InitalFrame));
set(handles.text13,'String','running ...');  pause(0.2);
axes(handles.axes1); colormap(gray);
tic;
try
    for k=startframe : nFrames
        I1 = imsubtract(tempmask,rgb2gray( read(handles.VideoObj, k) ));
        SearchRect = checkSearchRect(Xrange,Yrange,xydata,TestRadium);
        %I1(480:590,1270:1400)=0;
        %I1(590:610,1300:1360)=0;
        %I1(460:480,1300:1360)=0;
        %I1(500:600,520:670)=0;
        %I1(480:500,570:630)=0;
        %I1(600:620,560:640)=0;
        I2 = mat2gray( imcrop(I1, SearchRect) );
% %        
%         SE=strel('square',20); %pengzhang
%         I2=imdilate(I2,SE);
%         I2 = imerode (I2,SE);
%         I1=imdilate(I1,SE);
%         I1=imerode(I1,SE);
%%

         threshold = graythresh(I2)*1.5;
         A1 =im2bw(I2,threshold);
        %         A1= 1- imerode(A1,se1);
        
        t = bwboundaries(A1, 'noholes');
        boundaryArray =  cellfun('length',t) ;
        
        r = length(boundaryArray) ;
        if  r < 1
            imagesc(I1); axis image off;
            set(handles.text13,'String','select mice location');
            
            [x0,y0] = ginput(1);
            xydata = round([x0;y0]);
            
            axes(handles.axes5);
            imagesc(A1); axis image off;  title('A1');
            pause(3);
            imagesc(I2); axis image off; title('I2');
            pause(2);
            
            axes(handles.axes1);
            
            set(handles.text13,'String','running...');
            pause(0.1);
        else
            
            %         if r >1
            boundaryarea=zeros(length(t),1);
            for nn=1:length(t)
                xxx=t{nn};
                x1=[xxx(:,1);xxx(1,1)];
                y1=[xxx(:,2);xxx(1,2)];
                boundaryarea(nn)=polyarea(x1,y1);                
            end
            [val, i] =max(boundaryarea);
            t = t(i);
            %         end
            contourss =  t{:};
            [y0,x0] = centroid(contourss);
            x0 = round(x0); y0 = round(y0);
            
            xydata =  SearchRect(1:2) + [x0;y0];
        end
        Lengthofmove = 3*norm( xypos(:,k)- xydata);
        %         if Lengthofmove < Error_r
        %             statusString = 'Repeated track, stop by program';
        %             break;
        %         end
        xypos(:,k)=  xydata ;
        
        
        if mod(k,500)== 0
            if get(handles.radiobutton3,'Value')
                %             set(handles.text13,'String','running...');
                
                imagesc(I1); axis image off;
                if r >= 1
                    hold on; plot(contourss(:,2)+SearchRect(1),contourss(:,1)+SearchRect(2),'y'); hold off;
                end
                pb=patch(xydata(1)+xp,xydata(2)+yp,'r','edgecolor','r');  alpha(pb,0.5);
                % Set the frame index
                set(handles.ImageIndexSlider,'Value',k);
                set(handles.FrameIndex,'String', num2str(k));
                %Set the Video Time
                Time=roundn(k/FrameRate,-3);
                StringofTime=num2str(Time);
                set(handles.editTime,'String',StringofTime);
                %                 set(handles.text11,'String', num2str(threshold));
                pause(0.0001);
            end % of if get(handles.radiobutton3,'Value')
            
            % Check whether the pause button is pushed
            PauseButtonValue=get(handles.pushbuttonPause ,'Value');
            if PauseButtonValue
                text = get(handles.pushbuttonPause, 'String');
                if strcmp(text, 'Pause') == 1
                    set(handles.pushbuttonPause,'String','Continue');
                    uiwait;
                end
            end
            
            StopButtonValue=get(handles.pushbuttonStop,'Value');
            if StopButtonValue
                TrackCompleted = false;
                statusString = 'Stop by user';
                set(handles.pushbuttonStop,'Value',0);
                break; % Quit the displaying
            end
        end %of mod(k,500)== 0
        
    end
catch
    %do nothinbg here and continue to save the result.
    display(['tracking error occurs @ ' num2str(k)]);
    errmsg = lasterr;
    TrackCompleted = false;
    statusString = 'Computing error';
    if(strfind(errmsg, 'Bogus marker length'))
        disp('** Bogus marker length');
    end
end
toc;
%   delete(h);
handles.data_frame=k;


handles.xst=xst;
handles.yst=yst;
pathname=handles.pathname;
filename=handles.filename;
% if TrackCompleted
if TrackCompleted
    handles.xypos = xypos ;
    save([pathname,filename(1:end-4) '_track.mat'],'xypos','threshold','TestRadium');
%     handles.tempmask = tempmask;
    save ([handles.pathname, handles.filename(1:end-4) '_TempMask'], 'tempmask');        
end
handles.trackOk = TrackCompleted;

set(handles.text13,'String',statusString);
handles.threshold  = threshold ;
handles.TestRadium = TestRadium;

guidata(hObject,handles);

% % % % % %

% % % % % %
function SearchRect = checkSearchRect(Xrange,Yrange,xydata,TestRadium)
SearchRect = [xydata(1)-TestRadium; xydata(2)-TestRadium; TestRadium*2; TestRadium*2] ;
if SearchRect(1) <1  SearchRect(1) = 1;end;
if SearchRect(2) <1  SearchRect(2) = 1;end;
if (SearchRect(1) + SearchRect(3)) >= Xrange(2)   SearchRect(3) = Xrange(2) - SearchRect(1) -1; end;
if (SearchRect(2) + SearchRect(4)) >= Yrange(2)   SearchRect(4) = Yrange(2) - SearchRect(2) -1; end;

% --- Executes on button press in pushbuttonEllipsePlot.
function pushbuttonEllipsePlot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonEllipsePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1);

%import the mouse position data xypos
xypostion=handles.xypos;
% xypostion=xypostion(:,1:end-60);
data_frame=handles.data_frame;
%
xst=size(data_frame,1);
yst=size(data_frame,2);
imagesc(data_frame);
axis image off;
% xst=handles.xst;
% yst=handles.yst;

h = imellipse(gca, [20 20 200 200]);
% addNewPositionCallback(h,@(p) title(mat2str(p,3)));
fcn = makeConstrainToRectFcn('imellipse',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn);
block= wait(h);

pos = getPosition(h);
% hold on;

center_x=pos(1)+pos(3)/2;
center_y=pos(2)+pos(4)/2;
radius_xy(1)=pos(3)/2;
radius_xy(2)=pos(4)/2;
sita=0:pi/20:2*pi;fi=pi/4;
xdata=center_x+radius_xy(1)*cos(sita+fi);
ydata=center_y+radius_xy(2)*sin(sita+fi);
% plot(xdata,ydata,'g'); %椭圆，中心点在（x0,y0），长轴为a，短轴为b，方向角为fi。
delete(h);
% xdata=xdata-xst(1);
% ydata=ydata-yst(1);
xc=xdata';
yc=ydata';


% %in the center
% radius=300;
% center_x=(xst(end)-xst(1))/2;%%%%%% the center of the frame :x
% center_y=(yst(end)-yst(1))/2;%%%%%% the center of the frame :y
% radius_xy(1)=radius/2;
% radius_xy(2)=radius/2;
% coord = plot_rectangle(center_x,center_y,radius_xy,'hollow')';

%caculate the center

center_xst=xst/2;
center_yst=yst/2;
radius_xyst(1)=xst/2;
radius_xyst(2)=yst/2;
coord_xyst = plot_rectangle(center_xst,center_yst,radius_xyst,'k')';


% in_image = inpolygon(xypostion(1,:),xypostion(2,:),coord(1,:),coord(2,:));
in_image = inpolygon(xypostion(1,:),xypostion(2,:),xc,yc);
image_pixels=xypostion(:,in_image==1);
time_in_zone=size(image_pixels,2)*(1/30);
total_time=size(xypostion,2)*(1/30);
percentage=(time_in_zone/total_time)*100;
%
axes(handles.axes1);
p=figure;
set(gca,'FontSize',25);
% imagesc(data_frame);
axis image;
% axis square;
axis ij;
axis off;
hold on;
plot(coord_xyst(1,:),coord_xyst(2,:),'k','linewidth',2);
plot(xypostion(1,:),xypostion(2,:),'r');
plot(xc,yc,'b','linewidth',2);
title(['Time in the zone: ' num2str(time_in_zone) ' s. Total: '  num2str(total_time) ' s. Percentage: ' num2str(percentage) ' (%).']);
pathname=handles.pathname;
filename=handles.filename;
save([pathname,filename(1:end-4) '_centertime.mat'],'time_in_zone','total_time','percentage');
% close(p);
axis image;
saveas(p,[pathname,filename(1:end-4) '_centertime.fig']);
hold off;
guidata(hObject,handles);
axes(handles.axes1);

function editVideoPath_Callback(hObject, eventdata, handles)
% hObject    handle to editVideoPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editVideoPath as text
%        str2double(get(hObject,'String')) returns contents of editVideoPath as a double

% --- Executes during object creation, after setting all properties.
function editVideoPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVideoPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonPause.
function pushbuttonPause_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ButtonValue=get(handles.pushbuttonPause,'Value');
% switch ButtonValue
%     case 1
%          set(handles.pushbuttonPause,'String','Continue');
%     case 0
%          set(handles.pushbuttonPause,'String','Pause');
% end
text = get(handles.pushbuttonPause, 'String');
if strcmp(text, 'Continue') == 1
    uiresume(gcf);
    set(handles.pushbuttonPause,'String','Pause');
end

% --- Executes on button press in pushbuttonStop.
function pushbuttonStop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function StartTimeFrame_Callback(hObject, eventdata, handles)
% hObject    handle to StartTimeFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StartTimeFrame as text
%        str2double(get(hObject,'String')) returns contents of StartTimeFrame as a double

% handles.Control.InputMouseTime=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function StartTimeFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StartTimeFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonExit.
function pushbuttonExit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ss=questdlg('Are you sure you want to exit?','Exit Message Window?','No, I still want to use it!','Yes, I want to quit!','Yes, I want to quit!');
switch ss
    case 'Yes, I want to quit!'
        delete(handles.figure1);
end

% --- Executes on slider movement.
function ImageIndexSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ImageIndexSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
guidata(hObject,handles);
SliderValue=round(get(handles.ImageIndexSlider,'Value'));


axes(handles.axes1);

NumberofFrames = handles.VideoObj.NumberofFrames;
FrameRate=handles.VideoObj.FrameRate;
% Set the sliding control of the video
SliderStep(1) = 1/(NumberofFrames-1);
SliderStep(2) = 10/(NumberofFrames-1);
set(handles.ImageIndexSlider,'sliderstep',SliderStep,'max',NumberofFrames,'min',1,'Value',1);

% Display the image
ImageDataFrame = read(handles.VideoObj, SliderValue);
imagesc(ImageDataFrame);
axis image;
axis off;
if handles.trackOk
    Xoffset = handles.Xrange(1);
    Yoffset = handles.Yrange(1);
    ang=0:0.01:2*pi;
    r=5;
    xp=r*cos(ang)+Xoffset; yp=r*sin(ang)+Yoffset;
    x0=handles.xypos(1,SliderValue);
    y0=handles.xypos(2,SliderValue);
    pb=patch(x0+xp,y0+yp,'g','edgecolor','g');
    alpha(pb,0.5);
end


% Set the frame index
set(handles.ImageIndexSlider,'Value',SliderValue);
set(handles.FrameIndex,'String', num2str(SliderValue));
%Set the Video Time
Time=roundn(SliderValue/FrameRate,-3);
StringofTime=num2str(Time);
set(handles.editTime,'String',StringofTime);

guidata(hObject,handles);
axes(handles.axes1);
%
% --- Executes during object creation, after setting all properties.
function ImageIndexSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageIndexSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function FrameIndex_Callback(hObject, eventdata, handles)
% hObject    handle to FrameIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameIndex as text
%        str2double(get(hObject,'String')) returns contents of FrameIndex as a double
guidata(hObject,handles);
axes(handles.axes1);

% Set the frame according to the input data
FrameofIndex=round( str2double(get(hObject,'String')));
if FrameofIndex <   handles.VideoObj.NumberofFrames
    
    % handles.Control.FrameIndex=FrameofIndex;
    %
    % NumberofFrames = handles.VideoObj.NumberofFrames;
    % FrameRate=handles.VideoObj.FrameRate;
    % % Set the sliding control of the video
    % SliderStep(1) = 1/(NumberofFrames-1);
    % SliderStep(2) = 2/(NumberofFrames-1);
    % set(handles.ImageIndexSlider,'sliderstep',SliderStep,'max',NumberofFrames,'min',1,'Value',1);
    % ImageDataFrame = read(handles.VideoObj, FrameofIndex);
    % imagesc(ImageDataFrame);
    % axis image;
    % axis off;
    % Set the slider
    set(handles.ImageIndexSlider,'Value',FrameofIndex);
    ImageIndexSlider_Callback(hObject, eventdata, handles);
    % %Set the Video Time
    % Time=roundn(FrameofIndex/FrameRate,-3);
    % StringofTime=num2str(Time);
    % set(handles.editTime,'String',StringofTime);
end
guidata(hObject,handles);
% axes(handles.axes1);

% --- Executes during object creation, after setting all properties.
function FrameIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonAnyRegionTime.
function pushbuttonAnyRegionTime_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAnyRegionTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
axes(handles.axes1);

data_frame=handles.data_frame;
xst=size(data_frame,1);
yst=size(data_frame,2);
imagesc(data_frame);
axis image off;
% xst=handles.xst;
% yst=handles.yst;

%plot an any region by hand
f = imfreehand;
hc = get(f,'Children');
XData = []; YData = [];
for ii=1:length(hc)
    x = get(hc(ii),'XData');
    y = get(hc(ii),'YData');
    XData = [XData; x(:)];
    YData = [YData; y(:)];
end
delete(f);
% XData = XData%-xst(1);
% YData = YData-yst(1);
xc=XData';yc=YData';

%caculate the center
center_xst=xst/2;
center_yst=yst/2;
radius_xyst(1)=xst/2;
radius_xyst(2)=yst/2;
coord_xyst = plot_rectangle(center_xst,center_yst,radius_xyst,'k')';

%import the mouse position data xypos
xypostion=handles.xypos;
FrameRate=handles.VideoObj.FrameRate;
% xypostion=xypostion(:,1:end-60);
%caculate the percent
in_image = inpolygon(xypostion(1,:),xypostion(2,:),xc,yc);
image_pixels=xypostion(:,in_image==1);
time_in_zone=size(image_pixels,2)*(1/FrameRate);
total_time=size(xypostion,2)*(1/FrameRate);
percentage=(time_in_zone/total_time)*100;
%
p=figure;
set(gca,'FontSize',25);
axis image;
% axis square;
axis ij;
axis off;
hold on;
plot(xc,yc,'b','LineWidth',2);
plot(coord_xyst(1,:),coord_xyst(2,:),'k','linewidth',2);
plot(xypostion(1,:),xypostion(2,:),'r');
title(['Time in the zone: ' num2str(time_in_zone) ' s. Total: '  num2str(total_time) ' s. Percentage: ' num2str(percentage) ' (%).']);
pathname=handles.pathname;
filename=handles.filename;
save([pathname,filename(1:end-4) '_AnyRegionTime.mat'],'time_in_zone','total_time','percentage');
% close(p);
axis image;
saveas(p,[pathname,filename(1:end-4) '_AnyRegionTime.fig']);
hold off;
guidata(hObject,handles);
axes(handles.axes1);

% --- Executes on button press in pushbuttonRectPlot.
function pushbuttonRectPlot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRectPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

pathname=handles.pathname;
filename=handles.filename;

%caculate the center
center_xst=floor(handles.Xrange(2)/2);
center_yst=floor(handles.Yrange(2)/2);
radius_xyst(1)=floor(handles.Xrange(2)/2);
radius_xyst(2)=floor(handles.Yrange(2)/2);

coord_xyst = plot_rectangle(center_xst,center_yst,radius_xyst,'k')';

%import the mouse position data xypos
startFrame = round(str2double(get(handles.StartTimeFrame,'String')));
EndTimeFrame = round(str2double(get(handles.EndTimeFrame,'String')));
xypostion=handles.xypos(:,startFrame:EndTimeFrame);

% to plot the figure
p=figure;
set(gca,'FontSize',25);
axis image;
% axis square;
axis ij;
axis off;
hold on;
plot(coord_xyst(1,:),coord_xyst(2,:),'k','linewidth',2);
% % plot(coord(1,:),coord(2,:),'b','LineWidth',2);
plot(xypostion(1,:),xypostion(2,:),'b');%plot corner rectangles on the button

% save([pathname,filename(1:end-4) '_Rect.mat'],'time_in_zone','total_time','time_in_corner','percentage');
% close(p);
axis image;
set(gcf,'Position',[85 500 810 590]);
hold off;
saveas(p,[pathname,filename(1:end-4) '_Rect.fig']);
saveas(p,[pathname,filename(1:end-4) '_Rect.eps']);

guidata(hObject,handles);
axes(handles.axes1);


% --- Executes during object creation, after setting all properties.
function editTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EndTimeFrame_Callback(hObject, eventdata, handles)
% hObject    handle to EndTimeFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EndTimeFrame as text
%        str2double(get(hObject,'String')) returns contents of EndTimeFrame as a double
% handles.Control.EndTime=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function EndTimeFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EndTimeFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MoveBack.
function MoveBack_Callback(hObject, eventdata, handles)
% hObject    handle to MoveBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
axes(handles.axes1);

xst=handles.xst;
yst=handles.yst;
% data_frame=handles.data_frame;
xypos=handles.xypos;
FrameRate=handles.VideoObj.FrameRate;
% Extract the information for the screen
% handles.Control.InputMouseTime=str2double(get(handles.StartTimeFrame,'String'));
% handles.Control.EndTimeFrame=str2double(get(handles.EndTimeFrame,'String'));
% StratFrame=1;
% EndFrame=handles.Control.EndTime*FrameRate;
% move back
nFrames=handles.VideoObj.NumberOfFrames;
for k = 1 : nFrames
    % Display the image
    movdata = read(handles.VideoObj, k);
    data_frame=movdata(yst,xst,:);
    imagesc(data_frame);
    axis image off;
    hold on;
    %plot a circle to show how it works
    x1=xypos(1,k);
    y1=xypos(2,k);
    %       viscircles([x1 y1],r);
    sita=0:0.01:2*pi;         %角度[0,2*pi]
    R=10;                                 %半径
    x=R*cos(sita)+x1;
    y=R*sin(sita)+y1;
    %      plot(x,y,'r-');
    pb=patch(x,y,'r','edgecolor','r');
    alpha(pb,0.5);
    axis image off;
    hold off;
    
    IndexofFrame=k;
    % Set the frame index
    set(handles.ImageIndexSlider,'Value',IndexofFrame);
    set(handles.FrameIndex,'String', num2str(IndexofFrame));
    %Set the Video Time
    Time=roundn(IndexofFrame/FrameRate,-3);
    StringofTime=num2str(Time);
    set(handles.editTime,'String',StringofTime);
    
    % Check whether the pause button is pushed
    PauseButtonValue=get(handles.pushbuttonPause ,'Value');
    if PauseButtonValue
        text = get(handles.pushbuttonPause, 'String');
        if strcmp(text, 'Pause') == 1
            set(handles.pushbuttonPause,'String','Continue');
            uiwait;
        end
    end
    
    % Check whether the stop button is pushed
    StopButtonValue=get(handles.pushbuttonStop,'Value');
    if StopButtonValue
        movdata = read(handles.VideoObj, IndexofFrame);
        data_frame=movdata(yst,xst,:);
        imagesc(data_frame);
        axis image;
        axis off;
        hold on;
        pb=patch(x,y,'g','edgecolor','g');
        alpha(pb,0.5);
        hold off;
        set(handles.ImageIndexSlider,'Value',IndexofFrame);
        set(handles.FrameIndex,'String',num2str(IndexofFrame));
        set(handles.editTime,'String',StringofTime);
        set(handles.pushbuttonStop,'Value',0);
        % Quit the displaying
        break;
    end
    %
    %      % Check whether the slider button is pushed
    %     SliderValue=get(handles.ImageIndexSlider ,'Value');
    %     % Display the image
    %     ImageDataFrame = read(handles.VideoObj, SliderValue);
    %     imagesc(ImageDataFrame);
    %     axis off;
    %     % Set the frame index
    % %     set(handles.ImageIndexSlider,'Value',SliderValue);
    %     set(handles.FrameIndex,'String', num2str(SliderValue));
    %     %Set the Video Time
    %     Time=floor(SliderValue/30);
    %     StringofTime=num2str(Time);
    %     set(handles.editTime,'String',StringofTime);
end


function Window_Callback(hObject, eventdata, handles)
% hObject    handle to Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Window as text
%        str2double(get(hObject,'String')) returns contents of Window as a double


% --- Executes during object creation, after setting all properties.
function Window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadSignal.
function LoadSignal_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
FrameRate = handles.VideoObj.FrameRate;
nFrames = handles.VideoObj.NumberOfFrames;
AVIFileName = handles.filename; %downsample

% to select the signal file and do downsample
[filename, pathname] = uigetfile([AVIFileName(1:(end-4)) '*Sectiondata*.mat'], 'Select Signal Txt file');
if isequal(filename,0)
    disp('User selected Cancel')
else
    axes(handles.axes3); cla;
    if filename(1:5) == AVIFileName(1:5)
        filematchControl = true;
    end
    
    startFrameString = get(handles.StartTimeFrame,'String');
    if filematchControl
        if isempty( strfind(AVIFileName,'-')) %F for complete avi file, T for sectioned avi file
            
            filematchControl = ~isempty(strfind(filename,['-' startFrameString '-'])) ;
            if ~filematchControl
                set(handles.text13,'String','Wrong StartFrame !' );
            else
                set(handles.text13,'String','' );
            end
        else
            %           pass the startFrame test for sectioned avi
        end
    end
    
    
    if filematchControl
        startFrame = round(str2double(startFrameString));
        
        %    disp(['User selected ', fullfile(pathname, filename)])
        FolderNameFull = fullfile(pathname, filename) ; %[pathname '\' filename];
        set(handles.text13,'String',filename );
        % Signal = importdata(FolderNameFull);
        load(FolderNameFull);
        %     BaseLineCorrected= DataImport.BaseLineCorrected;
        if  exist('fs','var')
            handles.fs = fs ;
        else
            handles.fs = length(SelectionTimeData)/SelectionTimeData(end);
        end
        Number=fix(handles.fs/FrameRate);
        
        
        DataLength = length(SelectionTimeData);
        indexRow = 1:Number:DataLength;%     start to ReSample the Data from 2K into video frequency, 50/30/25 Hz
        
        
        if (length(indexRow) ~= (nFrames -startFrame + 1  ) )
            %         msgbox('Array size not match after ReSample','ERROR');
            if length(indexRow) < (nFrames - startFrame)
                %             data is shorter than the video length, choose part of the video
                endframe  =  startFrame + length(indexRow)-1;
                set(handles.EndTimeFrame,'String',num2str(endframe));
            else
                %             data is longer than the video length, choose part of the data
                indexRow  = indexRow(1,1:(nFrames -startFrame + 1 ));
            end
        else
            endframe = round(str2double(get(handles.EndTimeFrame,'String')));
        end
        handles.NewSignal = SelectedDataadjusted(indexRow);
        
        
        set(handles.Threshold,'String',num2str(BaseLineCorrected));
        %     SelectionTimeData = DataImport.SelectionTimeData;
        plot(SelectionTimeData,SelectedDataadjusted,'r' );
        hold on;
        
        YLine = BaseLineCorrected*ones(1,DataLength);
        handles.baselineHandle = plot(SelectionTimeData,YLine,'--b','LineWidth',2 );
        hold off;
        axis tight
        
        handles.SelectionTimeData = SelectionTimeData ;
        handles.FolderNameFull=FolderNameFull;
        handles.BaseLineCorrected = BaseLineCorrected ;
    else
        set(handles.text13,'String','Unmatched avi & data file names' );
    end
end

guidata(hObject,handles);




% --- Executes on button press in MergePlot.
function MergePlot_Callback(hObject, eventdata, handles)
% hObject    handle to MergePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);

z=handles.NewSignal;
set(handles.edit11,'String',num2str(max(z)));
xypos = handles.xypos;
FolderNameFull = handles.FolderNameFull;
startFrame = round(str2double(get(handles.StartTimeFrame,'String')));
endFrame = round(str2double(get(handles.EndTimeFrame,'String')));
x= xypos(1,startFrame:endFrame);
y=xypos(2,startFrame:endFrame);
% to get the information of points above the threshold
threshold=str2double(get(handles.Threshold,'String'));%0.064; %the threshold for signal
SelectedIndex = find(z >= threshold) ;
zAdjust = z; zAdjust(:)= NaN;
zAdjust(SelectedIndex) = z(SelectedIndex);

%
% xypostion=handles.xypos;
xypostionCount = zeros(handles.Xrange(2),handles.Yrange(2));
xypostionFire = xypostionCount ;
for i = 1: size(xypostionCount,1)
    k = find(x(:)==i);
    for j = 1: size(xypostionCount,2)
        ValueAtcurrentPosition = zAdjust(k( find(y(k)==j)));
        xypostionFire(i,j) = sum( ValueAtcurrentPosition(find(~isnan(ValueAtcurrentPosition))));
        xypostionCount(i,j) = length(ValueAtcurrentPosition);
    end
    xypostionCount(i,find(xypostionCount(i,:)==0)) = NaN;
end
xypostionFireAverage = xypostionFire./xypostionCount ;

% hf=fspecial('disk',handles.TestRadium/2);
% xypostionMap = filter2(hf,xypostionFireAverage,'same');
r = handles.TestRadium/2 ; UsedR = 2*r;
xlimit = size(xypostionCount,1);  ylimit  = size(xypostionCount,2);
% increase the size of result array to fit the boundary requirement
xypostionFireAverageVirtual = NaN(xlimit+2*r,ylimit+2*r);
xypostionFireAverageVirtual ((r+1):(xlimit+r),(r+1):(ylimit+r)) = xypostionFireAverage;
% prepare the average Algorithm
%     hf=fspecial('gaussian',2*r+1, 4);
hf=fspecial('disk',r);

% prepare  mark circle
xypostionMap = zeros(xlimit, ylimit);
for i=1: xlimit
    for j= 1: ylimit
        tempArray = xypostionFireAverageVirtual(i:(i+UsedR),j:(j+UsedR));
        %         tempArray = tempArray.*hf;
        xypostionMap(i,j) = mean(tempArray(find(~isnan(tempArray))));
    end
end
% % % % % % % % % % % % % % % % % % % % % % % %
xypostionMap = xypostionMap';
p = figure;
imagesc(xypostionMap);axis ij;%plot corner rectangles on the button
title(get(handles.text13,'String'));
axis image;
axis off;

saveas(p,[FolderNameFull(1:end-4) '_averageFireMapV53n.fig']);
saveas(p,[FolderNameFull(1:end-4) '_averageFireMapV53n.png']);
%
guidata(hObject,handles);
% axes(handles.axes3);



function Threshold_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Threshold as text
%        str2double(get(hObject,'String')) returns contents of Threshold as a double
guidata(hObject,handles);

% axes(handles.baselineHandle);
YLine = ones(1,length(handles.SelectionTimeData))*str2double(get(hObject,'String'));
set(handles.baselineHandle,'YData',YLine);
% updateFig2Threshold(handles);
% MergePlot_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function Threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in CalculatePath.
function CalculatePath_Callback(hObject, eventdata, handles)
% hObject    handle to CalculatePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);

xypos = handles.xypos;
FrameRate=handles.VideoObj.FrameRate;
% Extract the information for the screen
StartFrame = round(str2double(get(handles.StartTimeFrame,'String')));
EndFrame = round( str2double(get(handles.EndTimeFrame,'String')));
LengthofPath=0;
for j=StartFrame+1:EndFrame
    lengthpath=sqrt((xypos(1,j)-xypos(1,j-1))^2+(xypos(2,j)-xypos(2,j-1))^2);
    LengthofPath=LengthofPath+lengthpath;
end
handles.LengthofPath = LengthofPath;
set(handles.edit_lengthofpath,'String',num2str(LengthofPath));
guidata(hObject,handles);



function edit_lengthofpath_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lengthofpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lengthofpath as text
%        str2double(get(hObject,'String')) returns contents of edit_lengthofpath as a double


% --- Executes during object creation, after setting all properties.
function edit_lengthofpath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lengthofpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in aigo.
function aigo_Callback(hObject, eventdata, handles)
% hObject    handle to aigo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aigo


% --- Executes on button press in canon.
function canon_Callback(hObject, eventdata, handles)
% hObject    handle to canon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of canon


% --- Executes on button press in pushbutton_checkError.
function pushbutton_checkError_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_checkError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);

set(handles.text10,'String','total errors');
if isfield(handles,'xypos')
    xypos = handles.xypos;
    FrameRate=handles.VideoObj.FrameRate;
    nFrames = handles.VideoObj.NumberOfFrames;
    % Extract the information for the screen
    
    StartFrame = round( str2double(get(handles.StartTimeFrame,'String')) );
    EndFrame = round( str2double(get(handles.EndTimeFrame,'String')) );
    if StartFrame < 1
        StartFrame = 1;
    end
    if EndFrame >= nFrames
        EndFrame = nFrames;
    end
    
    Lengthofmove=zeros(1, EndFrame-StartFrame);
    for j=StartFrame+1:EndFrame
        Lengthofmove(j) = norm( xypos(:,j)- xypos(:,j-1));
        %     lengthpath=sqrt((xypos(1,j)-xypos(1,j-1))^2+(xypos(2,j)-xypos(2,j-1))^2);
        %     LengthofPath=LengthofPath+lengthpath;
    end
    Error_r = str2double(get(handles.edit_r,'String')) ;
    ErrorList = find(Lengthofmove > Error_r);
    set(handles.text10,'String',num2str(length(ErrorList)));
    set(handles.text11,'String','0');
    
    handles.ErrorList = ErrorList;
    
    guidata(hObject,handles);
end



function edit_r_Callback(hObject, eventdata, handles)
% hObject    handle to edit_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_r as text
%        str2double(get(hObject,'String')) returns contents of edit_r as a double


% --- Executes during object creation, after setting all properties.
function edit_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2NextError.
function pushbutton2NextError_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2NextError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
currentError = round(str2double(get(handles.text11,'String')));
totalErr = str2double(get(handles.text10,'String'));
currentError = mod(currentError, totalErr) + 1;

set(handles.text11,'String',num2str(currentError));

errorFrame = handles.ErrorList(currentError);
if errorFrame < handles.VideoObj.NumberOfFrames
    errorFrame = errorFrame;
end

set(handles.ImageIndexSlider,'Value',errorFrame);
ImageIndexSlider_Callback(hObject, eventdata, handles);

guidata(hObject,handles);


% --- Executes on button press in pushbutton_manualTrack.
function pushbutton_manualTrack_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_manualTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
Xoffset = handles.Xrange(1)-1;
Yoffset = handles.Yrange(1)-1;

Cretiriaofmove = str2double(get(handles.edit_r,'String')) ;
errorFrame = round(get(handles.ImageIndexSlider,'Value'));
tempString = get(handles.text11,'String');

CheckContinuous = true;
while CheckContinuous
    displayPrevLocation(errorFrame-1,handles);
    [Lx,Ly] = ginput(1);
    %   xpos =  handles.xypos(1,:) ;
    %   ypos = handles.xypos(2,:) ;
    %   r = round (get(handles.slider1,'Value'));
    
    
    
    handles.xypos(1,errorFrame) = round(Lx) - Xoffset; %adjust back with offset
    handles.xypos(2,errorFrame) = round(Ly) - Yoffset;
    if get(handles.radiobutton8,'Value')
        Distance2judge = norm( handles.xypos(:,errorFrame)- handles.xypos(:,errorFrame+1));
        CheckContinuous = (errorFrame <  handles.VideoObj.NumberOfFrames) & (Distance2judge > Cretiriaofmove) ;
        
        if CheckContinuous
            errorFrame = errorFrame +1;
            set(handles.ImageIndexSlider,'Value',errorFrame);
            ImageIndexSlider_Callback(hObject, eventdata, handles);
        else
            set(handles.text11,'String',num2str(Distance2judge));
            pause(2);
            set(handles.text11,'String', tempString);
            break;
        end
        
    else
        CheckContinuous = false;
    end
end
ImageIndexSlider_Callback(hObject, eventdata, handles);
guidata(hObject,handles);

function displayPrevLocation(prevFrame,handles)
% guidata(hObject,handles);
if prevFrame > 0
    Xoffset = handles.Xrange(1);
    Yoffset = handles.Yrange(1);
    ang=0:0.01:2*pi;
    r=10;
    xp=r*cos(ang)+Xoffset; yp=r*sin(ang)+Yoffset;
    x0=handles.xypos(1,prevFrame);
    y0=handles.xypos(2,prevFrame);
    pb=patch(x0+xp,y0+yp,'b','edgecolor','b');
    alpha(pb,0.5);
end

% --- Executes on button press in pushbutton2RewriteTrackFile.
function pushbutton2RewriteTrackFile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2RewriteTrackFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);

pathname=handles.pathname;
filename=handles.filename;
xypos = handles.xypos;
% Xrange = handles.Xrange; Yrange = handles.Yrange;
threshold = handles.threshold ;
TestRadium = handles.TestRadium ;
save([pathname,filename(1:end-4) '_track.mat'],'xypos','threshold','TestRadium');
% xypos = transpose(xypos);
% save([pathname,filename(1:end-4) '_trajectory.txt'],'xypos','-ascii' );

fid = fopen([pathname,filename(1:end-4) '_trajectory.txt'],'wt');
for i=1:size(xypos,2)
    fprintf(fid,'%d',xypos(1,i));
    fprintf(fid,'% d \n',xypos(2,i));
end
fclose(fid);
set(handles.text13,'String','track file saved');

% --- Executes on button press in pushbuttonCheckCurrentF.
function pushbuttonCheckCurrentF_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCheckCurrentF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);
if isfield(handles,'tempmask')
    tempmask = handles.tempmask ;
else
    TempMaskfile = [handles.pathname, handles.filename(1:end-4) '_TempMask.mat'];
    if exist(TempMaskfile,'file')==2
        load(TempMaskfile);
        handles.tempmask = tempmask ;
    else
        set(handles.text13,'String','No TempMask loaded, returned');
        return;
    end
end
% cla;
Yrange=handles.Yrange;
Xrange=handles.Xrange;
CurrentFrame = round(get(handles.ImageIndexSlider,'Value'));
% movdata = read(handles.VideoObj, CurrentFrame);
I1 = imsubtract(tempmask,rgb2gray( read(handles.VideoObj, CurrentFrame) ));
% I1=rgb2gray(movdata);

if isfield(handles,'xypos')
    
    xydata = handles.xypos(:,CurrentFrame) ;
    if isfield(handles,'TestRadium')
        TestRadium = handles.TestRadium ;
        SearchRect = checkSearchRect(Xrange,Yrange,xydata,TestRadium);
        I2 = mat2gray( imcrop(I1, SearchRect) );
    end
else
    return;
end

threshold = graythresh(I2)*1.3;
A1 =im2bw(I2,threshold);
axes(handles.axes5); imshow(A1);
%         threshold = graythresh(I2)*1.5;
%         A1 =im2bw(I2,threshold);
            se1=strel('disk',6);
%         A2= imerode(A1,se1);
% A2 = imclose(A1,se1);
t = bwboundaries(A1, 'noholes');
boundaryArray =  cellfun('length',t) ;

r = length(boundaryArray) ;

if r >=1
    [val, i] =max(boundaryArray);
    t = t(i);
    contourss =  t{:};
    [y0,x0] = centroid(contourss);
    hold on; plot(contourss(:,2),contourss(:,1),'g'); hold off;
    ang=0:0.01:2*pi;    r=10;  xp=r*cos(ang);yp=r*sin(ang) ;
    pb=patch(x0+xp,y0+yp,'r','edgecolor','r');  alpha(pb,0.5);
end

%         x0 = round(x0); y0 = round(y0);
%          A1(y0,x0) = 0.5;
%         imagesc(A1); axis image off;
title(num2str([x0 y0]));
pause(2);

% guidata(hObject,handles);
%     dbstop;



% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in pushbutton2CheckMotion.
function pushbutton2CheckMotion_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2CheckMotion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
xypos = handles.xypos;

nFrames = handles.VideoObj.NumberOfFrames;
FrameRate=handles.VideoObj.FrameRate;
% Extract the information for the screen
EndFrame = length(xypos);
if get(handles.radiobutton_MotionCenter,'Value') %calcute the distance relative to the center of the image
    temparray = 1: EndFrame;
    Lengthofmove=zeros(1, length(temparray));
    %     Tcenter = [round(mean(handles.Xrange));round(mean(handles.Yrange))];
    %     for j=1: length(Lengthofmove)
    %         Lengthofmove(j) = norm( xypos(:,temparray(j))- Tcenter);
    %     end
    %code above is for the motion distance relative to the center of the image;
    %code below is for the motion distance relative to the 1st point of mouse;
    for j=1: length(Lengthofmove)
        Lengthofmove(j) = norm( xypos(:,temparray(j))-  xypos(:,temparray(1)));
    end
    
    Titlewords = 'Moves by frames';
else
    MoveSteps = 0.2*FrameRate ; % every second
    MoveSteps = 1;
    temparray = 1: MoveSteps: EndFrame;
    Lengthofmove=zeros(1, length(temparray)-1);
    for j=1: length(Lengthofmove)
        Lengthofmove(j) = norm( xypos(:,temparray(j+1))- xypos(:,temparray(j)));
    end
    Titlewords = 'Moves every seconds';
end
timemotion=[1:1:length(Lengthofmove)]/25/60;
figure; plot(timemotion,Lengthofmove);
title(Titlewords);
xlabel('Time (min)')

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function updateFig2Threshold(handles)

PeakCValue = str2double(get(handles.edit11,'String'));
threshold = str2double(get(handles.Threshold,'String'));

th0 = get(handles.pSignalThresholdScatter,'Children');
axes((th0(2)));
% axes(handles.pSignalThresholdScatter);
% ThresholdRange = get(gca,'Clim');
set (gca,'Clim',[threshold PeakCValue]);

% tempTitle = ['Calcium Fluorescence Intensity Plot ( threshold= ' num2str(threshold) ' )'] ;
% set(get(th0(1),'Title'),'String',tempTitle);
% here it's not necessary to update the threshold




function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
guidata(hObject,handles);

updateFig2Threshold(handles);
% PeakCValue = str2double(get(hObject,'String'));
% AllRange = get(handles.gSignalAllScatter,'Clim');
% set (handles.gSignalAllScatter,'Clim',[AllRange(1) PeakCValue]);
% th0 = get(handles.pSignalThresholdScatter,'Children');
%  axes((th0(2)));
% % axes(handles.pSignalThresholdScatter);
% ThresholdRange = get(gca,'Clim');
% set (gca,'Clim',[ThresholdRange(1) PeakCValue]);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2save.
function pushbutton2save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);

% th0 = get(handles.pSignalThresholdScatter,'Children');
%  axes((th0(2)));
FolderNameFull = handles.FolderNameFull ;
tic;
saveas(handles.pSignalThresholdScatter,[FolderNameFull(1:end-4) '_Threshold_Scatter.png']);
saveas(handles.gSignalAllScatter,[FolderNameFull(1:end-4) '_Scatter.png']);
toc;

% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in pushbutton2Batch.
function pushbutton2Batch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2Batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton_MotionCenter.
function radiobutton_MotionCenter_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_MotionCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_MotionCenter


% --- Executes on button press in pushbutton_partCorrect.
function pushbutton_partCorrect_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_partCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
xypos = handles.xypos;
Yrange=handles.Yrange;
Xrange=handles.Xrange;
FrameRate=handles.VideoObj.FrameRate;
nFrames = handles.VideoObj.NumberOfFrames;
EndFrame = nFrames;
if EndFrame >= nFrames
    EndFrame = nFrames;
end
threshold = handles.threshold  ;
TestRadium = handles.TestRadium ;

% prepare the working enviroment
clear ContourStoreSt;
ContourStoreSt{1}=[];

se1=strel('disk',6);
k = get(handles.ImageIndexSlider,'Value') ;
xydata = xypos(:,k);
Lengthofmove  = 10000;
Cretiriaofmove = str2double(get(handles.edit_r,'String')) ;
set(handles.text13,'String','part running ...');
while Lengthofmove >= Cretiriaofmove
    
    k = k +1;
    I1 = rgb2gray( read(handles.VideoObj, k) );
    SearchRect = [xydata(1)-TestRadium; xydata(2)-TestRadium; TestRadium*2; TestRadium*2] ;
    if SearchRect(1) <1  SearchRect(1) = 1;end;
    if SearchRect(2) <1  SearchRect(2) = 1;end;
    
    if (SearchRect(1) + SearchRect(3)) >= Xrange(2)   SearchRect(3) = Xrange(2) - SearchRect(1) -1; end;
    if (SearchRect(2) + SearchRect(4)) >= Yrange(2)   SearchRect(4) = Yrange(2) - SearchRect(2) -1; end;
    
    I2 = mat2gray( imcrop(I1, SearchRect) );
    A1 =im2bw(I2,threshold);
    A1= 1- imerode(A1,se1);
    t = bwboundaries(A1, 'noholes');
    boundaryArray =  cellfun('length',t) ;
    
    r = length(boundaryArray) ;
    if  r < 1
        imagesc(I1); axis image off;
        set(handles.text13,'String','select mice location');
        [x0,y0] = ginput(1);
        xydata = round([x0;y0]);
        set(handles.text13,'String','part running...');
        pause(0.1);
    else
        if r >1
            [val, i] =max(boundaryArray);
            t = t(i);
        end
        contourss =  t{:};
        [y0,x0] = centroid(contourss);
        x0 = round(x0); y0 = round(y0);
        
        xydata =  SearchRect(1:2) + [x0;y0];
    end
    xypos(:,k) = xydata;
    
    % Set the frame index
    set(handles.ImageIndexSlider,'Value',k);
    set(handles.FrameIndex,'String', num2str(k));
    %Set the Video Time
    Time=roundn(k/FrameRate,-3);
    StringofTime=num2str(Time);
    set(handles.editTime,'String',StringofTime);
    pause(0.0001);
    
    Lengthofmove = norm( xypos(:,k)- xypos(:,k-1));
    %     if Lengthofmove > Cretiriaofmove
    %         set(handles.ImageIndexSlider,'Value',k);
    %         ImageIndexSlider_Callback(hObject, eventdata, handles);
    %         break; %current tracking couldn't solve the problem, break
    %     end
    if (k+1) >= EndFrame
        break;
    else
        Lengthofmove = norm( xypos(:,k+1)- xypos(:,k));
    end
end
I1 = rgb2gray( read(handles.VideoObj, k) );
imagesc(I1); axis image off;

handles.xypos = xypos ;
set(handles.text13,'String','part finished');
guidata(hObject,handles);


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton8


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

pathname=handles.pathname;
filename=handles.filename;
nFrames = handles.VideoObj.NumberOfFrames;

%caculate the center
center_xst=floor(handles.Xrange(2)/2);
center_yst=floor(handles.Yrange(2)/2);
radius_xyst(1)=floor(handles.Xrange(2)/2);
radius_xyst(2)=floor(handles.Yrange(2)/2);

coord_xyst = plot_rectangle(center_xst,center_yst,radius_xyst,'k')';

%import the mouse position data xypos
startFrame = round(str2double(get(handles.StartTimeFrame,'String')));
EndTimeFrame = round(str2double(get(handles.EndTimeFrame,'String')));
xypostion=handles.xypos(:,startFrame:EndTimeFrame);
% xypostion=handles.xypos;
xypostionCount = zeros(handles.Xrange(2),handles.Yrange(2));
for i = 1: size(xypostionCount,1)
    k = find(xypostion(1,:)==i);
    for j = 1: size(xypostionCount,2)
        xypostionCount(i,j) = length( find(xypostion(2,k)==j));
    end
end

hf=fspecial('gaussian',100,10);%social 200,20, NOPR 150,15
xypostionMap = filter2(hf,xypostionCount,'same');
xypostionMap = xypostionMap';

% figure; imagesc(xypostionMap);axis ij;

% figure; mesh(xypostionMap);axis ij;
% to plot the figure
p=figure;
set(gca,'FontSize',25);
colormap jet
% axis image;
% axis square;

% hold on;
% plot(coord_xyst(1,:),coord_xyst(2,:),'k','linewidth',2);
% % plot(coord(1,:),coord(2,:),'b','LineWidth',2);
imagesc(xypostionMap);axis ij;;%plot corner rectangles on the button
axis image;
axis off;
% save([pathname,filename(1:end-4) '_Rect.mat'],'time_in_zone','total_time','time_in_corner','percentage');
% close(p);
% axis image;
% set(gcf,'Position',[85 500 810 590]);

saveas(p,[pathname,filename(1:end-4) '_TimeDistribution.fig']);
saveas(p,[pathname,filename(1:end-4) '_TimeDistribution.png']);
% hold off;
guidata(hObject,handles);
axes(handles.axes1);


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
% axes(handles.axes1);
pathname=handles.pathname;
filename=handles.filename;
pix=getframe(handles.axes1);
imwrite(pix.cdata,[pathname,filename(1:end-4)  '_IMAGE.tif']);


% --- Executes on button press in pushbutton36.
function pushbutton36_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
handles.tempmask =0;
TempMaskfile = [handles.pathname, handles.filename(1:end-4) '_TempMask.mat'];
if exist(TempMaskfile,'file')==2
    delete(TempMaskfile)
end
set(handles.text13,'String','Mask data deleted');
guidata(hObject,handles);


% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
handles.trackOk = false;
trackfilename = [handles.pathname, handles.filename(1:end-4) '_track.mat'];
if exist(trackfilename,'file')==2
    delete(trackfilename)
end
set(handles.text13,'String','Track data deleted');
guidata(hObject,handles);


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);

tempmask = handles.tempmask;
CurrentFrame = round(get(handles.ImageIndexSlider,'Value'));
I1 = rgb2gray(read(handles.VideoObj, CurrentFrame));    

if size(tempmask) < 100
    handles.tempmask = I1 ;
else
    I2 = uint8( tempmask > I1) ;
    handles.tempmask = tempmask.*I2 + I1.*(1-I2);
end
axes(handles.axes5); imshow(handles.tempmask);

%     handles.tempmask = tempmask;
%     save ([handles.pathname, handles.filename(1:end-4) '_TempMask'], 'tempmask');
guidata(hObject,handles);


% --- Executes on button press in pushbutton39.
function pushbutton39_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles);
Xoffset = handles.Xrange(1)-1;
Yoffset = handles.Yrange(1)-1;
errorFrame = round(get(handles.ImageIndexSlider,'Value'));
displayPrevLocation(errorFrame-1,handles);

    [Lx,Ly] = ginput(1);    
    
    handles.xypos(1,errorFrame) = round(Lx) - Xoffset; %adjust back with offset
    handles.xypos(2,errorFrame) = round(Ly) - Yoffset;
    
ImageIndexSlider_Callback(hObject, eventdata, handles);

guidata(hObject,handles);

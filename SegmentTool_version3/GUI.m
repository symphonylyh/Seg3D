function varargout = GUI(varargin)
%GUI M-file for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('Property','Value',...) creates a new GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to GUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUI('CALLBACK') and GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 25-May-2016 15:29:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_inputimage.
function pushbutton_inputimage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_inputimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Input_top;
global Input_middle;
global Input_bottom;
global adj_top;
global adj_middle;
global adj_bottom;
global FileName;
global r1t;
global r2t;
global r3t;
global r1m;
global r2m;
global r3m;
global r1b;
global r2b;
global r3b;

%initialize ratio #1,#2,#3 that'll be used in calculating IBFI
r1t=0;
r2t=0;
r3t=0;
r1m=0;
r2m=0;
r3m=0;
r1b=0;
r2b=0;
r3b=0;

[FileName,PathName] = uigetfile({'*.CR2;*.png;*.jpg;*.jpeg;*.tiff;*.tif'},...
                                 'Select an Image');
if FileName == 0
    h = msgbox({'No Image Selected' 'Please Select an Image to Preceed'},'Error');
else
    I = imread(strcat(PathName, FileName));
    
    if size(I,3) == 3;
        I=rgb2gray(I);
    end
    
    numrows=size(I,1);%crop the image into three layers
    numcolumns=size(I,2);
    crop_rows=floor(numrows/3);
    discard_rows=numrows-3*crop_rows;
    im_split = mat2cell(I, [crop_rows crop_rows crop_rows discard_rows], numcolumns);
    Input_top=im_split{1};
    Input_middle=im_split{2};
    Input_bottom=im_split{3};
    [~, FileName, ~] = fileparts(FileName); %
    mkdir(['', FileName]);
    cd(FileName);
    imwrite(Input_top,'top.jpg');
    imwrite(Input_middle,'middle.jpg');
    imwrite(Input_bottom,'bottom.jpg');
    cd ..;
    axes(handles.axestop);
    imshow(Input_top);
    axes(handles.axesmiddle);
    imshow(Input_middle);
    axes(handles.axesbottom);
    imshow(Input_bottom);
    
    adj_top = Input_top;
    adj_middle = Input_middle;
    adj_bottom = Input_bottom;
    
    figure, imshow(I);
    h = imellipse;
    %using imellipse to determine ballsize
    position = wait(h);
    txt=clipboard('paste');
    h = msgbox({'Calibration Ball Position:' 'The 3rd and the 4th elements of this array are suggested ballsize values' txt},'Calibration Ball Position');
       
end


% --- Executes on button press in pushbuttonibfit.
function pushbuttonibfit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonibfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ballSize;
global trockSize;
global m;
global n;
global r1t;
global r2t;
global r3t;
global r1m;
global r2m;
global r3m;
global r1b;
global r2b;
global r3b;

a_lower=0.8;
a_upper=6.25;
b_lower=0.6;
b_upper=9;

ms = trockSize(trockSize(:,4)== 1,:);
r = ms(:,1)/ballSize/ballSize; % !!! make sure ballSize is updated for each input image

[x,y] = find(r >= a_lower & r <= a_upper);

p1 = sum(ms(x,1))/m/n;
[x1,y] = find(r < a_lower & r >= b_lower);
[x2,y] = find(r < b_upper & r >= a_upper);
p2 = (sum(ms(x1,1))+sum(ms(x2,1)))/m/n;

r1t=p1*100;
r2t=p2*100;
r3t=p1*100+p2*100;

if r1m==0 && r1b==0
    x=1;
end
if (r1m==0 && r1b~=0)||(r1m~=0 && r1b==0)
    x=2;
end
if r1m~=0 && r1b~=0
    x=3;
end
average1=(r1t+r1m+r1b)/x;
average2=(r2t+r2m+r2b)/x;
average3=(r3t+r3m+r3b)/x;
set(handles.editr1top,'String',num2str(r1t));
set(handles.editr2top,'String',num2str(r2t));
set(handles.editr3top,'String',num2str(r3t));
set(handles.editaver1,'String',num2str(average1));
set(handles.editaver2,'String',num2str(average2));
set(handles.editaver3,'String',num2str(average3));
set(handles.editIBFI,'String',num2str(100-average3));



function editr1top_Callback(hObject, eventdata, handles)
% hObject    handle to editr1top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editr1top as text
%        str2double(get(hObject,'String')) returns contents of editr1top as a double


% --- Executes during object creation, after setting all properties.
function editr1top_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editr1top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editr2top_Callback(hObject, eventdata, handles)
% hObject    handle to editr2top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editr2top as text
%        str2double(get(hObject,'String')) returns contents of editr2top as a double


% --- Executes during object creation, after setting all properties.
function editr2top_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editr2top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editr1middle_Callback(hObject, eventdata, handles)
% hObject    handle to editr1middle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editr1middle as text
%        str2double(get(hObject,'String')) returns contents of editr1middle as a double


% --- Executes during object creation, after setting all properties.
function editr1middle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editr1middle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editr2middle_Callback(hObject, eventdata, handles)
% hObject    handle to editr2middle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editr2middle as text
%        str2double(get(hObject,'String')) returns contents of editr2middle as a double


% --- Executes during object creation, after setting all properties.
function editr2middle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editr2middle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editr1bottom_Callback(hObject, eventdata, handles)
% hObject    handle to editr1bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editr1bottom as text
%        str2double(get(hObject,'String')) returns contents of editr1bottom as a double


% --- Executes during object creation, after setting all properties.
function editr1bottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editr1bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editr2bottom_Callback(hObject, eventdata, handles)
% hObject    handle to editr2bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editr2bottom as text
%        str2double(get(hObject,'String')) returns contents of editr2bottom as a double


% --- Executes during object creation, after setting all properties.
function editr2bottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editr2bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_convexthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to slider_convexthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global Convex_threshold;
slider_val=get(hObject,'Value')-get(hObject,'Min');
Convex_threshold = slider_val;
set(handles.editconvexthreshold,'String',num2str(Convex_threshold));


% --- Executes during object creation, after setting all properties.
function slider_convexthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_convexthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editconvexthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to editconvexthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editconvexthreshold as text
%        str2double(get(hObject,'String')) returns contents of editconvexthreshold as a double


% --- Executes during object creation, after setting all properties.
function editconvexthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editconvexthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_ballsize_Callback(hObject, eventdata, handles)
% hObject    handle to slider_ballsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global ballSize;
slider_val=get(hObject,'Value')-get(hObject,'Min');
ballSize = floor(345*slider_val+5);
set(handles.editballsize,'String',num2str(ballSize));
ballSize=ballSize/1.687;


% --- Executes during object creation, after setting all properties.
function slider_ballsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_ballsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editballsize_Callback(hObject, eventdata, handles)
% hObject    handle to editballsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editballsize as text
%        str2double(get(hObject,'String')) returns contents of editballsize as a double


% --- Executes during object creation, after setting all properties.
function editballsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editballsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonconvex.
function pushbuttonconvex_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonconvex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global trockSize;
    global mrockSize;
    global brockSize;
    global negCount;
    global posCount;
    global m;
    global n;
    global Convex_Layer;
    global topL;
    global middleL;
    global bottomL;
    global BF_top;
    global BF_middle;
    global BF_bottom;
    global ballSize;
    
    Convex_threshold = 0.73;
    
    %lower=0.11;
    %upper=7.069;
    
if Convex_Layer==2
    L=topL;
    outImg=BF_top;
    L_color_label=L;
end
if Convex_Layer==3
    L=middleL;
    outImg=BF_middle;
    L_color_label=L;
end
if Convex_Layer==4
    L=bottomL;
    outImg=BF_bottom;
    L_color_label=L;
end
    %folder_name = uigetdir;
    numRock = max(max(L));
    rockSize = [];
    rockID = 1;
    %threshold = 0.91;
    negCount = 0;
    posCount = 0;
    [m,n] = size(outImg);
    
    maxW = n/5; % upper bound
    maxH = maxW;
    %numRock=30;
    for i = 1 : numRock
        wholeMask = L;
        [rows,cols] = find(L == i);
        %if(i==109)
            %assignin('base','rows',rows)
            %assignin('base','cols',cols)
        %end
        
        wholeMask(L~=i) = 0;
        A = max(min(rows)- 1, 1);
        B = min(max(rows)+ 1, m);
        C = max(min(cols)- 1, 1);
        D = min(max(cols)+ 1, n);

        height = max(rows) - min(rows)+ 1;
        width = max(cols)- min(cols) + 1;
        %centerX = sum(rows)/size(rows,1);
        %centerY = sum(cols)/size(cols,1);

        [K,area] = convhull(rows,cols);
        numPix = size(rows,1); %number of pixels;
        %if(i==109)
            %assignin('base','numPix',numPix)
            %assignin('base','area',area)
        %end
        ratio = numPix/area;
        
        rockSize(rockID,1) = numPix;
        rockSize(rockID,2) = area; %the area of the convex hull
        rockSize(rockID,3) = numPix/area; %ratio
        %singleMask = double(wholeMask(A:B, C:D));
        %singleRock = double(outImg(A:B,C:D)).*double(singleMask);
        if ratio < Convex_threshold %&& height < maxH && width < maxW
            negCount = negCount + 1;
            rockSize(rockID,4) = 0; 
            L_color_label(L==i) = 0;
            %if(i>98) 
             %   L_color_label(L==i) = 54;
            %end
            %if(i==109) 
                %fprintf('109 is false');
            %end
            %imwrite(singleMask, [folder_name,'/FalseSeg',int2str(negCount),'_org_',int2str(rockID),'.jpg'],'jpg');
        else
            posCount = posCount + 1;
            
            rockSize(rockID,4) = 1;
            %assignin('base','rockSize',rockSize)
            %assignin('base','rockSize',rockSize);
            r=rockSize(rockID,1)/((ballSize/2)^2*pi);
            
            %if(i>12) 
                %L_color_label(L==i) = 2;
            %else
                %L_color_label(L==i) = 49;
            %end
            
            if r<0.6 
                %r<0.6 
                %r<0.36
                L_color_label(L==i) = ceil(r/0.2);
                %L_color_label(L==i) = ceil(r/0.2);
                %L_color_label(L==i) = ceil(r/0.12);
            end
            
            if r>12.49999 %r>7.49999
                L_color_label(L==i) = 64;%bigger than 64 is 64;
                %L_color_label(L==i) = 64;
            end
            
            if r>=9.0 && r<=12.49999 %r>=3 && r<=7.49999
                L_color_label(L==i) = ceil((r-9.0)/4.5*9)+55;
                %LL_color_label(L==i) = ceil((r-3.0)/4.5*9)+55;
            end
            
            if r >= 0.8 && r <= 6.25 
                %r >= 0.8 && r <= 2.5
                %r >= 0.64 && r <= 6.25 
                L_color_label(L==i) = ceil((r-0.8)/5.45*36)+9;
                %L_color_label(L==i) = ceil((r-0.8)/1.7*36)+9;
                %L_color_label(L==i) = ceil((r-0.8)/5.45*36)+9;
            end
            if r < 0.8 && r >= 0.6 
                %r < 0.8 && r >= 0.6
                %r < 0.64 && r >= 0.36
                L_color_label(L==i) = ceil((r-0.6)/0.2*6)+3;
                %L_color_label(L==i) = ceil((r-0.6)/0.2*6)+3;
                %L_color_label(L==i) = ceil((r-0.36)/0.28*6)+3;
            end
            if r < 9.0 && r >= 6.25 %r < 3.0 && r >= 2.5
                L_color_label(L==i) = ceil((r-6.25)/2.75*10)+45;%2.5-3.0
                %L_color_label(L==i) = ceil((r-2.5)/0.5*10)+45;
            end

            
        end
        %imwrite(singleMask, [folder_name,'/TrueSeg',int2str(posCount),'_org_',int2str(rockID),'.jpg'],'jpg');
        rockID = rockID + 1;
    end
    L_color_label(1,1)=56;
    for i = 1 : numRock
        L_color_label(L==0) = 0;
    end
    load('BallastCMap','mycmap');
    Lrgb = label2rgb(L_color_label, mycmap, [0.7 0.0 0.7], 'noshuffle');
    figure, imshow(Lrgb), title('color labeled true/false segments')
    imshow(outImg,[]);
    hold on;
    himage = imshow(Lrgb);
    set(himage, 'AlphaData', 0.3);
    hold off;
    
    if Convex_Layer==2
        trockSize=rockSize;
        %param_list = [posCount negCount strel_size1 Convex_threshold];
        %dlmwrite([folder_name,'/parambottom.txt'], param_list, '-append');
        axes(handles.axessegtop);
        imshow(Lrgb);
    end
    if Convex_Layer==3
        mrockSize=rockSize;
        %param_list = [posCount negCount strel_size1 Convex_threshold];
        %dlmwrite([folder_name,'/parambottom.txt'], param_list, '-append');
        axes(handles.axessegmiddle);
        imshow(Lrgb);
    end
    if Convex_Layer==4
        brockSize=rockSize;
        %param_list = [posCount negCount strel_size1 Convex_threshold];
        %dlmwrite([folder_name,'/parambottom.txt'], param_list, '-append');
        axes(handles.axessegbottom);
        imshow(Lrgb);
    end
    
    %cd(FileName);
    %imwrite(Lrgb, 'Water_shed.jpg');
    %cd ..;
    %figure; plot(sort(rockSize(:,3)));
        


% --- Executes on button press in pushbuttonapplystreltop.
function pushbuttonapplystreltop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonapplystreltop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global BF_top;
    global strel_size1;
    global topL;
    global FileName;
    
    L = marker_watershed(BF_top, strel_size1);
    topL=L;
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    imshow(BF_top,[]);
    hold on;
    himage = imshow(Lrgb);
    set(himage, 'AlphaData', 0.3);
    hold off;
    axes(handles.axessegtop);
    imshow(Lrgb);
    cd(FileName);
    imwrite(Lrgb, 'Water_shed_top.jpg');
    cd ..;




% --- Executes on slider movement.
function sliderstreltop_Callback(hObject, eventdata, handles)
% hObject    handle to sliderstreltop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global strel_size1;
slider_val=get(hObject,'Value')-get(hObject,'Min');
strel_size1 = floor(55*slider_val+5);
set(handles.editstreltop,'String',num2str(strel_size1));


% --- Executes during object creation, after setting all properties.
function sliderstreltop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderstreltop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editstreltop_Callback(hObject, eventdata, handles)
% hObject    handle to editstreltop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editstreltop as text
%        str2double(get(hObject,'String')) returns contents of editstreltop as a double


% --- Executes during object creation, after setting all properties.
function editstreltop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editstreltop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderstrelmiddle_Callback(hObject, eventdata, handles)
% hObject    handle to sliderstrelmiddle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global strel_size1;
slider_val=get(hObject,'Value')-get(hObject,'Min');
strel_size1 = floor(55*slider_val+5);
set(handles.editstrelmiddle,'String',num2str(strel_size1));


% --- Executes during object creation, after setting all properties.
function sliderstrelmiddle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderstrelmiddle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderstrelbottom_Callback(hObject, eventdata, handles)
% hObject    handle to sliderstrelbottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global strel_size1;
slider_val=get(hObject,'Value')-get(hObject,'Min');
strel_size1 = floor(55*slider_val+5);
set(handles.editstrelbottom,'String',num2str(strel_size1));


% --- Executes during object creation, after setting all properties.
function sliderstrelbottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderstrelbottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editstrelmiddle_Callback(hObject, eventdata, handles)
% hObject    handle to editstrelmiddle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editstrelmiddle as text
%        str2double(get(hObject,'String')) returns contents of editstrelmiddle as a double


% --- Executes during object creation, after setting all properties.
function editstrelmiddle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editstrelmiddle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editstrelbottom_Callback(hObject, eventdata, handles)
% hObject    handle to editstrelbottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editstrelbottom as text
%        str2double(get(hObject,'String')) returns contents of editstrelbottom as a double


% --- Executes during object creation, after setting all properties.
function editstrelbottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editstrelbottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonapplystrelmiddle.
function pushbuttonapplystrelmiddle_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonapplystrelmiddle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global BF_middle;
    global strel_size1;
    global middleL;
    global FileName;
    
    L = marker_watershed(BF_middle, strel_size1);
    middleL=L;
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    imshow(BF_middle,[]);
    hold on;
    himage = imshow(Lrgb);
    set(himage, 'AlphaData', 0.3);
    hold off;
    axes(handles.axessegmiddle);
    imshow(Lrgb);
    cd(FileName);
    imwrite(Lrgb, 'Water_shed_middle.jpg');
    cd ..;



% --- Executes on button press in pushbuttonapplystrelbottom.
function pushbuttonapplystrelbottom_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonapplystrelbottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global BF_bottom;
    global strel_size1;
    global bottomL;
    global FileName;
    
    L = marker_watershed(BF_bottom, strel_size1);
    bottomL=L;
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    imshow(BF_bottom,[]);
    hold on;
    himage = imshow(Lrgb);
    set(himage, 'AlphaData', 0.3);
    hold off;
    axes(handles.axessegbottom);
    imshow(Lrgb);
    cd(FileName);
    imwrite(Lrgb, 'Water_shed_bottom.jpg');
    cd ..;


% --- Executes on slider movement.
function sliderspatialgaussian_Callback(hObject, eventdata, handles)
% hObject    handle to sliderspatialgaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global SigmaS;
slider_val=get(hObject,'Value')-get(hObject,'Min');
SigmaS = floor(9*slider_val+1);
set(handles.editspatialgaussian,'String',num2str(SigmaS));


% --- Executes during object creation, after setting all properties.
function sliderspatialgaussian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderspatialgaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderrangegaussian_Callback(hObject, eventdata, handles)
% hObject    handle to sliderrangegaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global SigmaR;
slider_val=get(hObject,'Value')-get(hObject,'Min');
SigmaR = floor(19*slider_val+1);
set(handles.editrangegaussian,'String',num2str(SigmaR));

% --- Executes during object creation, after setting all properties.
function sliderrangegaussian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderrangegaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editspatialgaussian_Callback(hObject, eventdata, handles)
% hObject    handle to editspatialgaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editspatialgaussian as text
%        str2double(get(hObject,'String')) returns contents of editspatialgaussian as a double


% --- Executes during object creation, after setting all properties.
function editspatialgaussian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editspatialgaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editrangegaussian_Callback(hObject, eventdata, handles)
% hObject    handle to editrangegaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editrangegaussian as text
%        str2double(get(hObject,'String')) returns contents of editrangegaussian as a double


% --- Executes during object creation, after setting all properties.
function editrangegaussian_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editrangegaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonbft.
function pushbuttonbft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonbft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global SigmaS;
    global SigmaR;
    global BF_top;
    global adj_top;
    
    tol = 0.01;
    Img=adj_top;
    
    % make odd
    if (mod(SigmaS,2) == 0)
      w  = SigmaS + 1;
    else
      w  = SigmaS;
    end
    [outImg, param] =  shiftableBF(double(Img), SigmaS, SigmaR, w, tol);
    %outImg = bilateral_filter(Img, SigmaS, SigmaR);
    
     BF_top=outImg;
     axes(handles.axessegtop);
     imshow(outImg,[]);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
global Input_top;
global adj_top;

x=get(handles.popupmenu1,'value');
if x==1
    adj_top=Input_top;
    axes(handles.axestop);
    imshow(Input_top);
end
if x==2
    y=imread('sample.jpg');
    adj_top=imhistmatch(Input_top,y);
    axes(handles.axestop);
    imshow(adj_top);
end
    



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
global Input_middle;
global adj_middle

x=get(handles.popupmenu2,'value');
if x==1
    adj_middle=Input_middle;
    axes(handles.axesmiddle);
    imshow(Input_middle);
end
if x==2
    y=imread('sample.jpg');
    adj_middle=imhistmatch(Input_middle,y);
    axes(handles.axesmiddle);
    imshow(adj_middle);
end


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
global Input_bottom;
global adj_bottom;

x=get(handles.popupmenu3,'value');
if x==1
    adj_bottom=Input_bottom;
    axes(handles.axesbottom);
    imshow(Input_bottom);
end
if x==2
    y=imread('sample.jpg');
    adj_bottom=imhistmatch(Input_bottom,y);
    axes(handles.axesbottom);
    imshow(adj_bottom);
end



% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonapplymaskmiddle.
function pushbuttonapplymaskmiddle_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonapplymaskmiddle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonapplymaskbottom.
function pushbuttonapplymaskbottom_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonapplymaskbottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonapplymasktop.
function pushbuttonapplymasktop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonapplymasktop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editaver1_Callback(hObject, eventdata, handles)
% hObject    handle to editaver1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editaver1 as text
%        str2double(get(hObject,'String')) returns contents of editaver1 as a double



% --- Executes during object creation, after setting all properties.
function editaver1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editaver1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editaver2_Callback(hObject, eventdata, handles)
% hObject    handle to editaver2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editaver2 as text
%        str2double(get(hObject,'String')) returns contents of editaver2 as a double


% --- Executes during object creation, after setting all properties.
function editaver2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editaver2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editaver3_Callback(hObject, eventdata, handles)
% hObject    handle to editaver3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editaver3 as text
%        str2double(get(hObject,'String')) returns contents of editaver3 as a double


% --- Executes during object creation, after setting all properties.
function editaver3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editaver3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on selection change in popupconvex.
function popupconvex_Callback(hObject, eventdata, handles)
% hObject    handle to popupconvex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupconvex contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupconvex
global Convex_Layer;
Convex_Layer=get(handles.popupconvex,'value');


% --- Executes during object creation, after setting all properties.
function popupconvex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupconvex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupbf.
function popupbf_Callback(hObject, eventdata, handles)
% hObject    handle to popupbf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupbf contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupbf
global popupbf;

popupbf=get(handles.popupbf,'value');


% --- Executes during object creation, after setting all properties.
function popupbf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupbf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonibfim.
function pushbuttonibfim_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonibfim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ballSize;
global mrockSize;
global m;
global n;
global r1t;
global r2t;
global r3t;
global r1m;
global r2m;
global r3m;
global r1b;
global r2b;
global r3b;

a_lower=0.8;
a_upper=6.25;
b_lower=0.6;
b_upper=9;

ms = mrockSize(mrockSize(:,4)== 1,:);
r = ms(:,1)/ballSize/ballSize; % !!! make sure ballSize is updated for each input image

[x,y] = find(r >= a_lower & r <= a_upper);

p1 = sum(ms(x,1))/m/n;
[x1,y] = find(r < a_lower & r >= b_lower);
[x2,y] = find(r < b_upper & r >= a_upper);
p2 = (sum(ms(x1,1))+sum(ms(x2,1)))/m/n;

if r1t==0 && r1b==0
    x=1;
end
if (r1t==0 && r1b~=0)||(r1t~=0 && r1b==0)
    x=2;
end
if r1t~=0 && r1b~=0
    x=3;
end
r1m=p1*100;
r2m=p2*100;
r3m=p1*100+p2*100;
average1=(r1t+r1m+r1b)/x;
average2=(r2t+r2m+r2b)/x;
average3=(r3t+r3m+r3b)/x;
set(handles.editr1middle,'String',num2str(r1m));
set(handles.editr2middle,'String',num2str(r2m));
set(handles.editr3middle,'String',num2str(r3m));
set(handles.editaver1,'String',num2str(average1));
set(handles.editaver2,'String',num2str(average2));
set(handles.editaver3,'String',num2str(average3));
set(handles.editIBFI,'String',num2str(100-average3));


% --- Executes on button press in pushbuttonibfib.
function pushbuttonibfib_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonibfib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ballSize;
global brockSize;
global m;
global n;
global r1t;
global r2t;
global r3t;
global r1m;
global r2m;
global r3m;
global r1b;
global r2b;
global r3b;

a_lower=0.8;
a_upper=6.25;
b_lower=0.6;
b_upper=9;

ms = brockSize(brockSize(:,4)== 1,:);
r = ms(:,1)/ballSize/ballSize; % !!! make sure ballSize is updated for each input image

[x,y] = find(r >= a_lower & r <= a_upper);

p1 = sum(ms(x,1))/m/n;
[x1,y] = find(r < a_lower & r >= b_lower);
[x2,y] = find(r < b_upper & r >= a_upper);
p2 = (sum(ms(x1,1))+sum(ms(x2,1)))/m/n;

r1b=p1*100;
r2b=p2*100;
r3b=p1*100+p2*100;

if r1m==0 && r1t==0
    x=1;
end
if (r1m==0 && r1t~=0)||(r1m~=0 && r1t==0)
    x=2;
end
if r1m~=0 && r1t~=0
    x=3;
end
average1=(r1t+r1m+r1b)/x;
average2=(r2t+r2m+r2b)/x;
average3=(r3t+r3m+r3b)/x;
set(handles.editr1bottom,'String',num2str(r1b));
set(handles.editr2bottom,'String',num2str(r2b));
set(handles.editr3bottom,'String',num2str(r3b));
set(handles.editaver1,'String',num2str(average1));
set(handles.editaver2,'String',num2str(average2));
set(handles.editaver3,'String',num2str(average3));
set(handles.editIBFI,'String',num2str(100-average3));



function editr3top_Callback(hObject, eventdata, handles)
% hObject    handle to editr3top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editr3top as text
%        str2double(get(hObject,'String')) returns contents of editr3top as a double


% --- Executes during object creation, after setting all properties.
function editr3top_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editr3top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editr3middle_Callback(hObject, eventdata, handles)
% hObject    handle to editr3middle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editr3middle as text
%        str2double(get(hObject,'String')) returns contents of editr3middle as a double


% --- Executes during object creation, after setting all properties.
function editr3middle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editr3middle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editr3bottom_Callback(hObject, eventdata, handles)
% hObject    handle to editr3bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editr3bottom as text
%        str2double(get(hObject,'String')) returns contents of editr3bottom as a double


% --- Executes during object creation, after setting all properties.
function editr3bottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editr3bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editIBFI_Callback(hObject, eventdata, handles)
% hObject    handle to editIBFI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editIBFI as text
%        str2double(get(hObject,'String')) returns contents of editIBFI as a double


% --- Executes during object creation, after setting all properties.
function editIBFI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editIBFI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonbfm.
function pushbuttonbfm_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonbfm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global SigmaS;
    global SigmaR;
    global BF_middle;
    global adj_middle;
    
    tol = 0.01;
    Img=adj_middle;
    
    % make odd
    if (mod(SigmaS,2) == 0)
      w  = SigmaS + 1;
    else
      w  = SigmaS;
    end
    [outImg, param] =  shiftableBF(double(Img), SigmaS, SigmaR, w, tol);
    %outImg = bilateral_filter(Img, SigmaS, SigmaR);
    
        BF_middle=outImg;
        axes(handles.axessegmiddle);
        imshow(outImg,[]);



% --- Executes on button press in pushbuttonbfb.
function pushbuttonbfb_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonbfb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global SigmaS;
    global SigmaR;
    global BF_bottom;
    global adj_bottom;
    
    tol = 0.01;
    Img=adj_bottom;
    
    % make odd
    if (mod(SigmaS,2) == 0)
      w  = SigmaS + 1;
    else
      w  = SigmaS;
    end
    [outImg, param] =  shiftableBF(double(Img), SigmaS, SigmaR, w, tol);
    %outImg = bilateral_filter(Img, SigmaS, SigmaR);

        BF_bottom=outImg;
        axes(handles.axessegbottom);
        imshow(outImg,[]);


% --- Executes on button press in pushbuttonbfa.
function pushbuttonbfa_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonbfa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global SigmaS;
    global SigmaR;
    global BF_top;
    global BF_middle;
    global BF_bottom;
    global adj_top;
    global adj_middle;
    global adj_bottom;
    
    tol = 0.01;
    Img=adj_top;
    
    % make odd
    if (mod(SigmaS,2) == 0)
      w  = SigmaS + 1;
    else
      w  = SigmaS;
    end
    [outImg, param] =  shiftableBF(double(Img), SigmaS, SigmaR, w, tol);
    %outImg = bilateral_filter(Img, SigmaS, SigmaR);
    
     BF_top=outImg;
     axes(handles.axessegtop);
     imshow(outImg,[]);
     
     %middle layer
     tol = 0.01;
    Img=adj_middle;
    
    % make odd
    if (mod(SigmaS,2) == 0)
      w  = SigmaS + 1;
    else
      w  = SigmaS;
    end
    [outImg, param] =  shiftableBF(double(Img), SigmaS, SigmaR, w, tol);
    %outImg = bilateral_filter(Img, SigmaS, SigmaR);
    
        BF_middle=outImg;
        axes(handles.axessegmiddle);
        imshow(outImg,[]);
        
        %bottom layer

    tol = 0.01;
    Img=adj_bottom;
    
    % make odd
    if (mod(SigmaS,2) == 0)
      w  = SigmaS + 1;
    else
      w  = SigmaS;
    end
    [outImg, param] =  shiftableBF(double(Img), SigmaS, SigmaR, w, tol);
    %outImg = bilateral_filter(Img, SigmaS, SigmaR);

        BF_bottom=outImg;
        axes(handles.axessegbottom);
        imshow(outImg,[]);
    
    

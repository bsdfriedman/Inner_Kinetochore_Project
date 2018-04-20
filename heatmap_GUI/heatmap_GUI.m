function varargout = heatmap_GUI(varargin)
% HEATMAP_GUI MATLAB code for heatmap_GUI.fig
%      HEATMAP_GUI, by itself, creates a new HEATMAP_GUI or raises the existing
%      singleton*.
%
%      H = HEATMAP_GUI returns the handle to a new HEATMAP_GUI or the handle to
%      the existing singleton*.
%
%      HEATMAP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HEATMAP_GUI.M with the given input arguments.
%
%      HEATMAP_GUI('Property','Value',...) creates a new HEATMAP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before heatmap_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to heatmap_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help heatmap_GUI

% Last Modified by GUIDE v2.5 15-Aug-2017 09:57:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @heatmap_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @heatmap_GUI_OutputFcn, ...
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


% --- Executes just before heatmap_GUI is made visible.
function heatmap_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to heatmap_GUI (see VARARGIN)

% Choose default command line output for heatmap_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes heatmap_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = heatmap_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function plane_slider_Callback(hObject, eventdata, handles)
% hObject    handle to plane_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% set slider_value to the plane_slider 'Value' GFP channel

%load in the existing lacO points
gfp_points = handles.gfp_points;
spb_points = handles.spb_points;

% update axes3 first so that GFP channel is "selected" for thresholding
slider_value = ceil(get(handles.plane_slider,'Value'));
%if you have SPB and GFP spots plot them when slider updates
if isempty(gfp_points) == 0 && isempty(spb_points) == 0
    subplot(handles.axes3);
    hold on;
    imshow(handles.im_mat_2(:,:,slider_value),...
        [get(handles.axes3_min_slider, 'Value'),...
        get(handles.axes3_max_slider, 'Value')]);
    plot(gfp_points(:,1), gfp_points(:,2), 'go');
    plot(spb_points(:,1), spb_points(:,2), 'rx');
    hold off;
        subplot(handles.axes1);
    hold on;
    imshow(handles.im_mat(:,:,slider_value),...
        [get(handles.axes1_min_slider, 'Value'),...
        get(handles.axes1_max_slider, 'Value')]);
    plot(gfp_points(:,1), gfp_points(:,2), 'go');
    plot(spb_points(:,1), spb_points(:,2), 'rx');
    hold off;
%if you have GFP spots only, plot them when slider updates
elseif isempty(gfp_points) == 0 && isempty(spb_points) == 1
    subplot(handles.axes3);
    hold on;
    imshow(handles.im_mat_2(:,:,slider_value),...
        [get(handles.axes3_min_slider, 'Value'),...
        get(handles.axes3_max_slider, 'Value')]);
    plot(gfp_points(:,1), gfp_points(:,2), 'go');
    hold off;
        subplot(handles.axes1);
    hold on;
    imshow(handles.im_mat(:,:,slider_value),...
        [get(handles.axes1_min_slider, 'Value'),...
        get(handles.axes1_max_slider, 'Value')]);
    plot(gfp_points(:,1), gfp_points(:,2), 'go');
    hold off;
    %if you have SPB spots only, plot them when slider updates
elseif isempty(gfp_points) == 1 && isempty(spb_points) == 0
    subplot(handles.axes3);
    hold on;
    imshow(handles.im_mat_2(:,:,slider_value),...
        [get(handles.axes3_min_slider, 'Value'),...
        get(handles.axes3_max_slider, 'Value')]);
    plot(spb_points(:,1), spb_points(:,2), 'rx');
    hold off;
        subplot(handles.axes1);
    hold on;
    imshow(handles.im_mat(:,:,slider_value),...
        [get(handles.axes1_min_slider, 'Value'),...
        get(handles.axes1_max_slider, 'Value')]);
    plot(spb_points(:,1), spb_points(:,2), 'rx');
    hold off;
%if you have no spots yet,do not plot when slider updates
 else
    subplot(handles.axes3);
    imshow(handles.im_mat_2(:,:,slider_value),...
        [get(handles.axes3_min_slider, 'Value'),...
        get(handles.axes3_max_slider, 'Value')]);
    subplot(handles.axes1);
    imshow(handles.im_mat(:,:,slider_value),...
        [get(handles.axes1_min_slider, 'Value'),...
        get(handles.axes1_max_slider, 'Value')]);
end


%set slider text
plane_num = handles.plane_num;
set(handles.slider_text, 'String', strcat(num2str(...
    ceil(get(handles.plane_slider,'Value'))),...
    '/', num2str(plane_num)));
%set Stack text
stack_total = plane_num/str2double(get(handles.step_number,'String'));
current_stack = ceil(slider_value/...
    str2double(get(handles.step_number,'String')));
set(handles.stack_text,'String',strcat(num2str(current_stack),...
    num2str('/'), num2str(stack_total)));

% --- Executes during object creation, after setting all properties.
function plane_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plane_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in gfp_1_button.
function gfp_1_button_Callback(hObject, eventdata, handles)
% hObject    handle to gfp_1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Switch to axes1
subplot(handles.axes1);
%wait for user button press on SPB
k = waitforbuttonpress;
click_point = round(get(gca,'CurrentPoint'));
hold on
plot(click_point(1,1),click_point(1,2), '-go');
hold off
%get the current plane to find Z position
current_plane=ceil(get(handles.plane_slider, 'Value'));
%figure out where we are in the z-stack
%get the step number information from the edit text box
step_number = str2double(get(handles.step_number,'String'));
%modulo of current step by step number
stack_mod = mod(ceil(get(handles.plane_slider,'Value')), step_number);
if stack_mod == 0
    stack_min = (1-step_number) + current_plane;
    stack_max = current_plane;
else
    stack_min = (1-stack_mod) + current_plane;
    stack_max = abs(stack_mod - step_number) + current_plane;
end

% let me display the search area and the resulting pixel
% on the command line to check
% current_plane
% stack_min
% stack_max
% select a portion of the stack based on the  up a a larger XYZ search area based on the click_point
if click_point(1,1)<5
    click_point(1,1)=5;
end
if click_point(1,2)<5
    click_point(1,2)=5;
end

%correcting if click region is close to maximum
im_mat = handles.im_mat;
[y,x,~] = size(im_mat);
if click_point(1,1) > (x-4)
    click_point(1,1) = x-4;
end
if click_point(1,2) > (y-4)
    click_point(1,2) = y-4;
end

% select a portion of the stack based on the  up a a larger XYZ search area based on the click_point
crop_im = handles.im_mat((click_point(1,2)-4):(click_point(1,2)+4),...
    (click_point(1,1)-4):(click_point(1,1)+4),stack_min:stack_max);
max_pix = max(crop_im(:));
[max_y, max_x, max_z] = ind2sub(size(crop_im),find(crop_im==max_pix));
if length(max_x)~= 1
    max_x=max_x(1);
    max_y=max_y(1);
    max_z=max_z(1);
end
%convert from crop_im coordinates back to the im_mat_2 coordinates
max_x_im = max_x + click_point(1,1)-4 - 1;
max_y_im = max_y + click_point(1,2)-4 - 1;
max_z_im = max_z + stack_min - 1;
test_max = handles.im_mat(max_y_im, max_x_im, max_z_im);
% %check that the intensity values match by displaying the crop and
% %the original
% test_max
% max_pix
lac1 = [max_x_im, max_y_im, max_z_im, max_pix];
gfp_points = [handles.gfp_points;[max_x_im,max_y_im]];
handles.gfp_points = gfp_points;
handles.lac1 = lac1;
guidata(hObject,handles);

% --- Executes on button press in gfp_2_button.
function gfp_2_button_Callback(hObject, eventdata, handles)
% hObject    handle to gfp_2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Switch to axes1
subplot(handles.axes1);
%wait for user button press on SPB
k = waitforbuttonpress;
click_point = round(get(gca,'CurrentPoint'));
hold on
plot(click_point(1,1),click_point(1,2), '-go');
hold off
%get the current plane to find Z position
current_plane=ceil(get(handles.plane_slider, 'Value'));
%figure out where we are in the z-stack
%get the step number information from the edit text box
step_number = str2double(get(handles.step_number,'String'));
%modulo of current step by step number
stack_mod = mod(ceil(get(handles.plane_slider,'Value')), step_number);
if stack_mod == 0
    stack_min = (1-step_number) + current_plane;
    stack_max = current_plane;
else
    stack_min = (1-stack_mod) + current_plane;
    stack_max = abs(stack_mod - step_number) + current_plane;
end

% let me display the search area and the resulting pixel
% on the command line to check
% current_plane
% stack_min
% stack_max
% select a portion of the stack based on the  up a a larger XYZ search area based on the click_point
if click_point(1,1)<5
    click_point(1,1)=5;
end
if click_point(1,2)<5
    click_point(1,2)=5;
end

%correcting if click region is close to maximum
im_mat = handles.im_mat;
[y,x,~] = size(im_mat);
if click_point(1,1) > (x-4)
    click_point(1,1) = x-4;
end
if click_point(1,2) > (y-4)
    click_point(1,2) = y-4;
end

% select a portion of the stack based on the  up a a larger XYZ search area based on the click_point
crop_im = handles.im_mat((click_point(1,2)-4):(click_point(1,2)+4),...
    (click_point(1,1)-4):(click_point(1,1)+4),stack_min:stack_max);
max_pix = max(crop_im(:));
[max_y, max_x, max_z] = ind2sub(size(crop_im),find(crop_im==max_pix));
if length(max_x)~= 1
    max_x=max_x(1);
    max_y=max_y(1);
    max_z=max_z(1);
end
%convert from crop_im coordinates back to the im_mat_2 coordinates
max_x_im = max_x + click_point(1,1)-4 - 1;
max_y_im = max_y + click_point(1,2)-4 - 1;
max_z_im = max_z + stack_min - 1;
test_max = handles.im_mat(max_y_im, max_x_im, max_z_im);
% %check that the intensity values match by displaying the crop and
% %the original
% test_max
% max_pix
lac2 = [max_x_im, max_y_im, max_z_im, max_pix];
gfp_points = [handles.gfp_points;[max_x_im,max_y_im]];
handles.gfp_points = gfp_points;
handles.lac2 = lac2;
guidata(hObject,handles);

% --- Executes on button press in spb_1_button.
function spb_1_button_Callback(hObject, eventdata, handles)
% hObject    handle to spb_1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Switch to axes3
subplot(handles.axes3);
%wait for user button press on SPB
k = waitforbuttonpress;
click_point = round(get(gca,'CurrentPoint'));
hold on
plot(click_point(1,1),click_point(1,2), '-rx');
hold off
%get the current plane to find Z position
current_plane=ceil(get(handles.plane_slider, 'Value'));
%figure out where we are in the z-stack
%get the step number information from the edit text box
step_number = str2double(get(handles.step_number,'String'));
%modulo of current step by step number
stack_mod = mod(ceil(get(handles.plane_slider,'Value')), step_number);
if stack_mod == 0
    stack_min = (1-step_number) + current_plane;
    stack_max = current_plane;
else
    stack_min = (1-stack_mod) + current_plane;
    stack_max = abs(stack_mod - step_number) + current_plane;
end

% let me display the search area and the resulting pixel
% on the command line to check
% current_plane
% stack_min
% stack_max

% select a portion of the stack based on the  up a a larger XYZ search area based on the click_point
if click_point(1,1)<5
    click_point(1,1)=5;
end
if click_point(1,2)<5
    click_point(1,2)=5;
end

%correcting if click region is close to maximum
im_mat = handles.im_mat;
[y,x,~] = size(im_mat);
if click_point(1,1) > (x-4)
    click_point(1,1) = x-4;
end
if click_point(1,2) > (y-4)
    click_point(1,2) = y-4;
end
crop_im = handles.im_mat_2((click_point(1,2)-4):(click_point(1,2)+4),...
    (click_point(1,1)-4):(click_point(1,1)+4),stack_min:stack_max);
max_pix = max(crop_im(:));
[max_y, max_x, max_z] = ind2sub(size(crop_im),find(crop_im==max_pix));
if length(max_x)~= 1
    max_x=max_x(1);
    max_y=max_y(1);
    max_z=max_z(1);
end
%convert from crop_im coordinates back to the im_mat_2 coordinates
max_x_im = max_x + click_point(1,1)-4 - 1;
max_y_im = max_y + click_point(1,2)-4 - 1;
max_z_im = max_z + stack_min - 1;
test_max = handles.im_mat_2(max_y_im, max_x_im, max_z_im);
% %check that the intensity values match by displaying the crop and
% %the original
% % test_max
% % max_pix
% save the X Y and Plane data for SPB1 to the handles structure
spb1 = [max_x_im, max_y_im, max_z_im, max_pix];
spb_points = [handles.spb_points;[max_x_im,max_y_im]];
handles.spb_points = spb_points;
handles.spb1 = spb1;
guidata(hObject,handles);

% --- Executes on button press in spb_2_button.
function spb_2_button_Callback(hObject, eventdata, handles)
% hObject    handle to spb_2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Switch to axes3
subplot(handles.axes3);
%wait for user button press on SPB
k = waitforbuttonpress;
click_point = round(get(gca,'CurrentPoint'));
hold on
plot(click_point(1,1),click_point(1,2), '-rx');
hold off
%get the current plane to find Z position
current_plane=ceil(get(handles.plane_slider, 'Value'));
%figure out where we are in the z-stack
%get the step number information from the edit text box
step_number = str2double(get(handles.step_number,'String'));
%modulo of current step by step number
stack_mod = mod(ceil(get(handles.plane_slider,'Value')), step_number);
if stack_mod == 0
    stack_min = (1-step_number) + current_plane;
    stack_max = current_plane;
else
    stack_min = (1-stack_mod) + current_plane;
    stack_max = abs(stack_mod - step_number) + current_plane;
end

% let me display the search area and the resulting pixel
% on the command line to check
% current_plane
% stack_min
% stack_max
% select a portion of the stack based on the  up a a larger XYZ search area based on the click_point
if click_point(1,1)<5
    click_point(1,1)=5;
end
if click_point(1,2)<5
    click_point(1,2)=5;
end

%correcting if click region is close to maximum
im_mat = handles.im_mat;
[y,x,~] = size(im_mat);
if click_point(1,1) > (x-4)
    click_point(1,1) = x-4;
end
if click_point(1,2) > (y-4)
    click_point(1,2) = y-4;
end

% select a portion of the stack based on the  up a a larger XYZ search area based on the click_point
crop_im = handles.im_mat_2((click_point(1,2)-4):(click_point(1,2)+4),...
    (click_point(1,1)-4):(click_point(1,1)+4),stack_min:stack_max);
max_pix = max(crop_im(:));
[max_y, max_x, max_z] = ind2sub(size(crop_im),find(crop_im==max_pix));
if length(max_x)~= 1
    max_x=max_x(1);
    max_y=max_y(1);
    max_z=max_z(1);
end
%convert from crop_im coordinates back to the im_mat_2 coordinates
max_x_im = max_x + click_point(1,1)-4 - 1;
max_y_im = max_y + click_point(1,2)-4 - 1;
max_z_im = max_z + stack_min - 1;
test_max = handles.im_mat_2(max_y_im, max_x_im, max_z_im);
% %check that the intensity values match by displaying the crop and
% %the original
% test_max
% max_pix
spb2 = [max_x_im, max_y_im, max_z_im, max_pix];
spb_points = [handles.spb_points;[max_x_im,max_y_im]];
handles.spb_points = spb_points;
handles.spb2 = spb2;
guidata(hObject,handles);


function step_number_Callback(hObject, eventdata, handles)
% hObject    handle to step_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step_number as text
%        str2double(get(hObject,'String')) returns contents of step_number as a double


% --- Executes during object creation, after setting all properties.
function step_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function step_size_Callback(hObject, eventdata, handles)
% hObject    handle to step_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step_size as text
%        str2double(get(hObject,'String')) returns contents of step_size as a double


% --- Executes during object creation, after setting all properties.
function step_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in open_image.
function open_image_Callback(hObject, eventdata, handles)
% hObject    handle to open_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set working directory to current directory
wdir = pwd;

%select the file
[filename pathname] = uigetfile('*.*', 'Choose the GFP image stack');
%read in the file using tif3Dread.m
im_cell = bfopen(strcat(pathname,filename));
im_mat = bf2mat(im_cell);
%change image to double
im_mat = double(im_mat);
%set contrast sliders for Axes1
set(handles.axes1_min_slider, 'SliderStep', [1/(max(im_mat(:))-min(im_mat(:))), 0.1],...
    'Min', min(im_mat(:)), 'Max', max(im_mat(:)), 'Value', min(im_mat(:)));
set(handles.axes1_max_slider, 'SliderStep', [1/(max(im_mat(:))-min(im_mat(:))), 0.1],...
    'Min', min(im_mat(:)), 'Max', max(im_mat(:)), 'Value', max(im_mat(:)));
%set up data to axes1 plot
subplot(handles.axes1);
%show first frame of image stack
imshow(im_mat(:,:,1),...
        [get(handles.axes1_min_slider, 'Value'),...
        get(handles.axes1_max_slider, 'Value')]);
%read in the number of frames
[~,~,plane_num] = size(im_mat);
handles.plane_num = plane_num;
%set slider
slider_step = 1/plane_num;
set(handles.plane_slider, 'SliderStep', [slider_step slider_step],...
    'Min', 1, 'Max', plane_num, 'Value', 1);
handles.im_mat = im_mat;
%set slider text
set(handles.slider_text, 'String', strcat(num2str(...
    get(handles.plane_slider,'Value')),...
    '/', num2str(plane_num)));
%set Stack text
stack_total = plane_num/str2double(get(handles.step_number,'String'));
current_stack = 1;
set(handles.stack_text,'String',strcat(num2str(current_stack),...
    num2str('/'), num2str(stack_total)));

%open the RFP stack in axes3
[filename pathname] = uigetfile('*.*', 'Choose the RFP image stack');
im_cell_2 = bfopen(strcat(pathname,filename));
im_mat_2 = bf2mat(im_cell_2);
im_mat_2 = double(im_mat_2);
%set contrast sliders for Axes3
set(handles.axes3_min_slider, 'SliderStep', [1/(max(im_mat_2(:))-min(im_mat_2(:))), 0.1],...
    'Min', min(im_mat_2(:)), 'Max', max(im_mat_2(:)), 'Value', min(im_mat_2(:)));
set(handles.axes3_max_slider, 'SliderStep', [1/(max(im_mat_2(:))-min(im_mat_2(:))), 0.1],...
    'Min', min(im_mat_2(:)), 'Max', max(im_mat_2(:)), 'Value', max(im_mat_2(:)));
handles.im_mat_2 = im_mat_2;
subplot(handles.axes3);
imshow(im_mat_2(:,:, ceil(get(handles.plane_slider,'Value'))),...
        [get(handles.axes3_min_slider, 'Value'),...
        get(handles.axes3_max_slider, 'Value')]);
%set up the cell array that will hold the image stack data
data_cell = {'lacO 1','lacO 1 Stretch','lacO 2','lacO 2 Stretch',...
    'SPB 1', 'SPB 2','Step Size (nm)',...
    'Pixel Size (nm)'};
handles.data_cell = data_cell;
record_num = 1;
handles.record_num = record_num;
%instantiate the handles.gfp_points and handles.spb_points
handles.gfp_points = [];
handles.spb_points = [];
% update the handles data structure
guidata(hObject,handles);


% --- Executes on button press in record_button.
function record_button_Callback(hObject, eventdata, handles)
% hObject    handle to record_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This button will record the data in some a strucuture that
%can be saved later
spb1 = handles.spb1;
spb2 = handles.spb2;
lac1 = handles.lac1;
lac1_stretch = get(handles.lac1_stretch,'Value');
lac2 = handles.lac2;
lac2_stretch = get(handles.lac2_stretch,'Value');
step_size = str2double(get(handles.step_size,'String'));
pixel_size = str2double(get(handles.pixel_size,'String'));
record_num = handles.record_num;

% put it all in the cell array that we created when we opened the image
handles.data_cell(record_num + 1,:) = {lac1,lac1_stretch,lac2,lac2_stretch,...
    spb1, spb2, step_size, ...
    pixel_size};
%increase the record number
record_num = record_num + 1;
handles.record_num = record_num;
%update the record display text
set(handles.record_text,'String', num2str(record_num));

%clear the cell-specific data
handles.lac1 = [];
handles.lac2 = [];
handles.spb1 = [];
handles.spb2 = [];


%revert text boxes to original states


%update the handles data structure
guidata(hObject,handles);


% --- Executes on button press in save_data_button.
function save_data_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_data_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%load the cell array from the handles structure
data_cell = handles.data_cell;
uisave('data_cell');



function pixel_size_Callback(hObject, eventdata, handles)
% hObject    handle to pixel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixel_size as text
%        str2double(get(hObject,'String')) returns contents of pixel_size as a double


% --- Executes during object creation, after setting all properties.
function pixel_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clear_lac.
function clear_lac_Callback(hObject, eventdata, handles)
% hObject    handle to clear_lac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.gfp_points = [];
guidata(hObject,handles);


% --- Executes on button press in lac1_stretch.
function lac1_stretch_Callback(hObject, eventdata, handles)
% hObject    handle to lac1_stretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lac1_stretch


% --- Executes on button press in lac2_stretch.
function lac2_stretch_Callback(hObject, eventdata, handles)
% hObject    handle to lac2_stretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lac2_stretch


% --- Executes on button press in clear_spbs.
function clear_spbs_Callback(hObject, eventdata, handles)
% hObject    handle to clear_spbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spb_points = [];
guidata(hObject,handles);


% --- Executes on slider movement.
function axes1_max_slider_Callback(hObject, eventdata, handles)
% hObject    handle to axes1_max_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function axes1_max_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1_max_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function axes3_max_slider_Callback(hObject, eventdata, handles)
% hObject    handle to axes3_max_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function axes3_max_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3_max_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function axes1_min_slider_Callback(hObject, eventdata, handles)
% hObject    handle to axes1_min_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function axes1_min_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1_min_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function axes3_min_slider_Callback(hObject, eventdata, handles)
% hObject    handle to axes3_min_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function axes3_min_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3_min_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

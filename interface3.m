function varargout = interface3(varargin)
% INTERFACE3 MATLAB code for interface3.fig
%      INTERFACE3, by itself, creates a new INTERFACE3 or raises the existing
%      singleton*.
%
%      H = INTERFACE3 returns the handle to a new INTERFACE3 or the handle to
%      the existing singleton*.
%
%      INTERFACE3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERFACE3.M with the given input arguments.
%
%      INTERFACE3('Property','Value',...) creates a new INTERFACE3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before interface3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to interface3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help interface3

% Last Modified by GUIDE v2.5 06-Sep-2017 19:49:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @interface3_OpeningFcn, ...
                   'gui_OutputFcn',  @interface3_OutputFcn, ...
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


% --- Executes just before interface3 is made visible.
function interface3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to interface3 (see VARARGIN)

% Choose default command line output for interface3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes interface3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = interface3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

[filename pathname] = uigetfile({'*.mat';'*.bmp'},'File Selector');
imageName = strcat(pathname, filename);
input_data = load(imageName);
% rgb = imread(input_data2.cjdata.image);
I=input_data.cjdata.image
axes(handles.axes1);
imshow(I,[]), title('Original Image');


hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
 figure
imshow(gradmag,[]), title('Step-1: Computing Gradient Magnitude')

L = watershed(gradmag);
Lrgb = label2rgb(L);
figure
imshow(Lrgb), title('Step-2: Applying Watershed Transform on Gradient Magnitude')

se = strel('disk', 20);
Io = imopen(I, se);
figure
imshow(Io), title('Step-3: Computing Opening')

Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
 figure
imshow(Iobr), title('Step-4: Computing Opening by Reconstruction')

Ioc = imclose(Io, se);
figure
imshow(Ioc), title('Step-5: Computing Closing Followed By Opening for Removing Dark Spots')

Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
 figure
imshow(Iobrcbr), title('Step-6: Computing Opening-Closing by Reconstruction')

fgm = imregionalmax(Iobrcbr);
figure
imshow(fgm), title('Step-7: Computing Regional Maxima of Opening-Closing by Reconstruction')

I2 = I;
I2(fgm) = 255;
figure
imshow(I2), title('Step-8: Computing Regional Maxima Superimposed on Original Image')

se2 = strel(ones(5,5));
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);

fgm4 = bwareaopen(fgm3, 20);
I3 = I;
I3(fgm4) = 255;
figure
imshow(I3)
title('Step-9: Computing Modified Regional Maxima Superimposed on Original Image')

gradmag2 = imimposemin(gradmag, fgm4);

%%
%%Finally we are ready to compute the watershed-based segmentation.

L = watershed(gradmag2);

I4 = I;
I4(imdilate(L == 0, ones(3, 3)) | fgm4) = 255;
figure
imshow(I4)
title('Markers and object boundaries superimposed on original image')
bw = edge(fgm,0.01);

axes(handles.axes2);
J2=showmask(double(I),imdilate(bw,ones(2,2))); imagesc(J2);axis off
title('Final Segmented Image');


%making the mask for region extraction.
m = uint8(fgm4); %converting the fgm matrix from logical values to uint8 
m = m*255; %scaling the values for visual representation
res = bitand(m,I); %applying the mask on the image

imgmat = graycomatrix(res); %extracting graylevel-cooccurance matrix from the extracted tumour region
stats = graycoprops(imgmat);%extracting 4 glcm features.. Contrast, Correlation, Energy, Homogeneity

cellarr = struct2cell(stats);%converting feature structure to matix. 
matarr = cell2mat(cellarr);

matarr = matarr';%converting the matrix from column to row

matarr = horzcat(matarr,1);%concatenating the label for the region as 1 as it is the tumour region. 

 dlmwrite('tfeat.txt',matarr, '-append',  'delimiter', ' ' );%writing the feature vector to file tfeat.txt
 
%now for the rest of the image and its features to be calculated and labeled as 0. 

%now to extract the image other than tumour part
%inverting the mask calculated in the previous step. 

m2 = 255 - m ; %subtracting each pixel value from the 255. so white becomes black and black becomes white. 

res2 = bitand(m2,I);%anding the mask for region other than tumour with the original image

%finding the feature vector from the extracted part of image other than the tumour part 
imgmat = graycomatrix(res2); %extracting graylevel-cooccurance matrix
stats = graycoprops(imgmat);%extracting 4 glcm features as above. 

cellarr = struct2cell(stats);%converting feature structure to matix. 
matarr2 = cell2mat(cellarr);

matarr2 = matarr2';%from column matrix to row. 

matarr2 = horzcat(matarr2,0);%concatenating the label for rest of region as 0. 
 dlmwrite('tfeat.txt',matarr2, '-append',  'delimiter', ' ' );

% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

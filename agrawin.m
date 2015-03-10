function varargout = agrawin(varargin)
% AGRAWIN MATLAB code for agrawin.fig
%      AGRAWIN, by itself, creates a new AGRAWIN or raises the existing
%      singleton*.
%
%      H = AGRAWIN returns the handle to a new AGRAWIN or the handle to
%      the existing singleton*.
%
%      AGRAWIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AGRAWIN.M with the given input arguments.
%
%      AGRAWIN('Property','Value',...) creates a new AGRAWIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before agrawin_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to agrawin_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help agrawin

% Last Modified by GUIDE v2.5 08-Mar-2015 22:32:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @agrawin_OpeningFcn, ...
                   'gui_OutputFcn',  @agrawin_OutputFcn, ...
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


% --- Executes just before agrawin is made visible.
function agrawin_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to agrawin (see VARARGIN)

% Choose default command line output for agrawin
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes agrawin wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global agra
limits=[agra.x_lim;agra.y_lim];
set(handles.space,'XGrid','on','YGrid','on','Box','on','XLim',limits(1,:),'YLim',limits(2,:));
set(handles.limits,'Data',limits);
set(handles.V,'Value',1);set(handles.dV,'Value',1);set(handles.ROA,'Value',0);

% --- Outputs from this function are returned to the command line.
function varargout = agrawin_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Draw.
function Draw_Callback(hObject, eventdata, handles)
% hObject    handle to Draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global agra;
cla;box on;grid on;hold on;
flag_V=get(handles.V,'Value');
flag_dV=get(handles.dV,'Value');
flag_ROA=get(handles.ROA,'Value');
if flag_ROA, agra.show_ROA(); end
if flag_dV,  agra.show_dVdt('r'); end
if flag_V, agra.show_V('b'); end
if flag_dV, agra.plot_critical('g'); end


% --- Executes when entered data in editable cell(s) in limits.
function limits_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to limits (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

limits=get(handles.limits,'Data');
set(handles.space,'XGrid','on','YGrid','on','Box','on','XLim',limits(1,:),'YLim',limits(2,:));
global agra;
agra.x_lim=limits(1,:);agra.y_lim=limits(2,:);

% --- Executes on button press in V.
function V_Callback(hObject, eventdata, handles)
% hObject    handle to V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of V


% --- Executes on button press in ROA.
function ROA_Callback(hObject, eventdata, handles)
% hObject    handle to ROA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ROA


% --- Executes on button press in dV.
function dV_Callback(hObject, eventdata, handles)
% hObject    handle to dV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dV


% --- Executes on button press in Trajectory.
function Trajectory_Callback(hObject, eventdata, handles)
% hObject    handle to Trajectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global agra;
agra.show_trajectory();

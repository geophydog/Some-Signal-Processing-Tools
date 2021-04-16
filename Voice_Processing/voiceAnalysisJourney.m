function varargout = voiceAnalysisJourney(varargin)
% VOICEANALYSISJOURNEY MATLAB code for voiceAnalysisJourney.fig
%      VOICEANALYSISJOURNEY, by itself, creates a new VOICEANALYSISJOURNEY or raises the existing
%      singleton*.
%
%      H = VOICEANALYSISJOURNEY returns the handle to a new VOICEANALYSISJOURNEY or the handle to
%      the existing singleton*.
%
%      VOICEANALYSISJOURNEY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOICEANALYSISJOURNEY.M with the given input arguments.
%
%      VOICEANALYSISJOURNEY('Property','Value',...) creates a new VOICEANALYSISJOURNEY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before voiceAnalysisJourney_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to voiceAnalysisJourney_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help voiceAnalysisJourney

% Last Modified by GUIDE v2.5 28-Apr-2019 20:06:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @voiceAnalysisJourney_OpeningFcn, ...
                   'gui_OutputFcn',  @voiceAnalysisJourney_OutputFcn, ...
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


% --- Executes just before voiceAnalysisJourney is made visible.
function voiceAnalysisJourney_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to voiceAnalysisJourney (see VARARGIN)

% Choose default command line output for voiceAnalysisJourney
handles.output = hObject;

logo = imread( 'AudioLogo.bmp' );
javaImage = im2java( logo );
newIcon = javax.swing.ImageIcon( javaImage );
figFrame = get( handles.figure1, 'JavaFrame' );
figFrame.setFigureIcon( newIcon );

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes voiceAnalysisJourney wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = voiceAnalysisJourney_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%                           LoadVoiceDataFile.                            %
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
function LoadVoiceDataFile_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to LoadVoiceDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName, pathName, filterIndex] = uigetfile({'*.wav; *.mp3'}, ...
    'Load audio file (s)', 'MultiSelect','on');%#ok

ftype = class(fileName);
if iscell(ftype)
    audioFile = strcat( [pathName, '\', cell2mat(fileName(1))]);
elseif ischar(ftype)
    audioFile = strcat( [pathName, '\', fileName]);
else
    errordlg('File type error!', 'File error!!!', 'modal');
end

% for i = 1: length(fileName)
%     fprintf('File: %s ...\n', cell2mat(fileName(i)));
% end

set( handles.LoadVoiceDataFile, 'UserData', audioFile);



% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%                           PlotAudioWaveform.                            %
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
function PlotAudioWaveform_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to PlotAudioWaveform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
audioFile = get(handles.LoadVoiceDataFile, 'UserData');
[s, fs] = audioread( audioFile );
nd = length(s);
t = linspace(0, nd-1, nd) / fs;

cla(handles.AudioWaveformsAxes);
axes(handles.AudioWaveformsAxes);
plot(t, s, 'w', 'linewidth', 0.5);
xlabel('Time [s]');
ylabel('Amplitude');
xlim([t(1), t(end)]);
set(gca, 'fontsize', 9, 'fontweight', 'bold', 'color', [0, 0, 0], ...
    'xcolor', 'w', 'ycolor', 'w');

data.fs = fs;
data.t = t;
data.s = s;
set( handles.PlotAudioWaveform, 'UserData', data);



% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%            Color map choice for time-frequency analysis                 %
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
function ColorMapChoice_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to ColorMapChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ColorMapChoice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ColorMapChoice
list = get( hObject, 'String' );
val1 = get( hObject, 'Value' );
cmap = list{ val1 };
if isempty( cmap )
    cmap = 'jet';
end
set( handles.ColorMapChoice, 'UserData', cmap );


% --- Executes during object creation, after setting all properties.
function ColorMapChoice_CreateFcn(hObject, eventdata, handles)%#ok
% hObject    handle to ColorMapChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%                          TimeFrequencyAnalysisPlot.                     %
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
function TimeFrequencyAnalysisPlot_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to TimeFrequencyAnalysisPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Voice data.
% Table data.
tableData = get( handles.uitable2, 'Data');
tableData = str2double(tableData);
tRange = tableData(1, :);
TFAFRange = tableData(2, :);

% Waveform data.
data = get( handles.PlotAudioWaveform, 'UserData');
t = data.t;
s = data.s;
dt = t(2) - t(1);

% Time range.
if ~isnan(tRange(1))
    t1 = tRange(1);
else
    t1 = t(1);
end
if ~isnan(tRange(2))
    t2 = tRange(2);
else
    t2 = t(end);
end
tn1 = floor((t1-t(1))/dt);
tn2 = floor((t2-t(1))/dt);
if tn1 < 1
    tn1 = 1;
end
if tn2 > length(t)
    tn2 = length(t);
end

% Get color map.
cmap = get(handles.ColorMapChoice, 'UserData');
if isempty(cmap)
    cmap = 'jet';
end

% Do time-fre. analysis.
s = s(tn1: tn2);
nd = length(s);
win = floor(nd/250);
overlap = floor(nd/300);
nfft = 2 ^ (nextpow2(win) + 1);
[p, F, T] = spectrogram(s, win, overlap, nfft, data.fs);
p = power(abs(p), 0.35);
T = T + t1;


% TFA frequency range.
if ~isnan(TFAFRange(1))
    tf1 = TFAFRange(1);
else
    tf1 = F(1);
end
if ~isnan(TFAFRange(2))
    tf2 = TFAFRange(2);
else
    tf2 = F(end);
end
fn1 = floor((tf1-F(1))/(F(2)-F(1)));
fn2 = floor((tf2-F(1))/(F(2)-F(1)));
if fn1 < 1
    fn1 = 1;
end
if fn2 > length(F)
    fn2 = length(F);
end

% Plot TFA results.
cla(handles.AudioTimeFrequencyAnalysisAxes);
axes(handles.AudioTimeFrequencyAnalysisAxes);
colormap(cmap);
imagesc(T, F(fn1: fn2), p(fn1: fn2, :));
xlabel('Time [s]');
ylabel('Frequency [Hz]');
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'color', [0, 0, 0], ...
    'xcolor', 'w', 'ycolor', 'w', 'ydir', 'normal');
rotate3d(handles.AudioTimeFrequencyAnalysisAxes, 'on');

% Update waveforms plot.
cla(handles.AudioWaveformsAxes);
axes(handles.AudioWaveformsAxes);
plot(t, data.s, 'w', 'linewidth', 0.5);
xlabel('Time [s]');
ylabel('Amplitude');
xlim([t1, t2]);
set(gca, 'fontsize', 9, 'fontweight', 'bold', 'color', [0, 0, 0], ...
    'xcolor', 'w', 'ycolor', 'w');





% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%                             FFTAmplitudeSpectrum.                       %
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
function FFTAmplitudeSpectrum_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to FFTAmplitudeSpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get( handles.PlotAudioWaveform, 'UserData');

% FFT: Amplitude spectrum.
nd = length(data.s);
fs = data.fs;
nfft = 2 ^ nextpow2(nd);
f = linspace(0, nfft-1, nfft) * fs / nfft;
fft_s = fft(data.s, nfft);
amp = abs(fft_s);
h_nfft = floor(nfft/2);

% Table data.
tableData = get( handles.uitable2, 'Data');
tableData = str2double(tableData);
AmpSpeFRange = tableData(3, :);
% Fre. range of amplitude spectrum.
if ~isnan(AmpSpeFRange(1))
    fn1 = floor((AmpSpeFRange(1)-f(1))/(f(2)-f(1)));
    if fn1 < 1
        fn1 = 1;
    end
else
    fn1 = 1;
end
if ~isnan(AmpSpeFRange(2))
    fn2 = floor((AmpSpeFRange(2)-f(1))/(f(2)-f(1)));
    if fn2 > h_nfft
        fn2 = h_nfft;
    end
else
    fn2 = h_nfft;
end

if fn1 > fn2
    errordlg('Fre. range error!', 'Data error!!!', 'modal');
    return;
end


% Plot amplitude spectrum.
cla(handles.FFTAmplitudeAxes);
axes(handles.FFTAmplitudeAxes);
plot(f(fn1: fn2), amp(fn1: fn2), 'y', 'linewidth', 0.75);
xlabel('Frequency [Hz]');
ylabel('Amplitude');
xlim([f(fn1), f(fn2)]);
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'color', [0, 0, 0], ...
    'xcolor', 'w', 'ycolor', 'w');


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%                               FFTPhaseSpectrum.                         %
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
function FFTPhaseSpectrum_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to FFTPhaseSpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get( handles.PlotAudioWaveform, 'UserData');

nd = length(data.s);
fs = data.fs;
nfft = 2 ^ nextpow2(nd);
f = linspace(0, nfft-1, nfft) * fs / nfft;
fft_s = fft(data.s, nfft);
pha = angle(fft_s) / pi * 180;
h_nfft = floor(nfft/2);

% Table data.
tableData = get( handles.uitable2, 'Data');
tableData = str2double(tableData);
PhaSpeFRange = tableData(4, :);
% Frequency range.
% Fre. range of phase spectrum.
if ~isnan(PhaSpeFRange(1))
    fn1 = floor((PhaSpeFRange(1)-f(1))/(f(2)-f(1)));
    if fn1 < 1
        fn1 = 1;
    end
else
    fn1 = 1;
end
if ~isnan(PhaSpeFRange(2))
    fn2 = floor((PhaSpeFRange(2)-f(1))/(f(2)-f(1)));
    if fn2 < 1
        fn2 = 1;
    end
else
    fn2 = h_nfft;
end
if fn1 > fn2
    errordlg('Fre. range error!', 'Data error!!!', 'modal');
    return;
end

cla(handles.FFTPhaseAxes);
axes(handles.FFTPhaseAxes);
plot(f(fn1: fn2), pha(fn1: fn2), 'y', 'linewidth', 0.75);
xlabel('Frequency [Hz]');
ylabel('Phase [^o]');
xlim([f(fn1), f(fn2)]);
set(gca, 'fontsize', 10, 'fontweight', 'bold', 'color', [0, 0, 0], ...
    'xcolor', 'w', 'ycolor', 'w');



% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%                            Listen to this Audio.                        %
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
function ListenAudio_Callback(hObject, eventdata, handles)%#ok
% hObject    handle to ListenAudio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get( handles.PlotAudioWaveform, 'UserData');
tableData = get( handles.uitable2, 'Data');
tableData = str2double(tableData);
tRange = tableData(1, :);

s = data.s;
t = data.t;

% Time range.
if ~isnan(tRange(1))
    tn1 = floor((tRange(1)-t(1))/(t(2)-t(1)));
    if tn1 < 1
        tn1 = 1;
    end
else
    tn1 = 1;
end
if ~isnan(tRange(2))
    tn2 = floor((tRange(2)-t(1))/(t(2)-t(1)));
    if tn2 > length(t)
        tn2 = length(t);
    end
else
    tn2 = length(t);
end

% Listen to this.
sound(s(tn1: tn2), data.fs);

clc;clear all;close all;
%{
This script will go on a process that takes out tidal components from raw
sea level data recorded at each tidal gauge. Those tidal components consist of 4 major tides and
background tides. 
%}

addpath(genpath('S:/home/user_006/08_MATLIB'))
spath = 'C:/Users/admin/Desktop/OBS_TIDE/figures/';

%% Select target tidal gauge
tg_tidal_gauge = '울릉도';

switch tg_tidal_gauge
    case '마산'
        lon_tidal_gauge = 128.5838; % 128.5769 (TMD doesn't represent MASAN obs sites owing to low resolution)
        lat_tidal_gauge = 35.1883; % 35.1928
        tg_tidal_gauge_Eng = 'MASAN';
    case '부산'
        lon_tidal_gauge = 129.0275;
        lat_tidal_gauge = 35.0942;
        tg_tidal_gauge_Eng = 'BUSAN';
    case '동해항'
        lon_tidal_gauge = 129.1938; % 129.1438
        lat_tidal_gauge = 37.4947;
        tg_tidal_gauge_Eng = 'PORT_DONGHAE';
    case '서귀포'
        lon_tidal_gauge = 126.5689;
        lat_tidal_gauge = 33.2306;
        tg_tidal_gauge_Eng = 'SEOGUIPO';
    case '성산포'
        lon_tidal_gauge = 126.9275;
        lat_tidal_gauge = 33.4806;
        tg_tidal_gauge_Eng = 'SEONGSANPO';
    case '울산'
        lon_tidal_gauge = 129.4025;
        lat_tidal_gauge = 35.4994;
        tg_tidal_gauge_Eng = 'ULSAN';
    case '묵호'
        lon_tidal_gauge = 129.1163;
        lat_tidal_gauge = 37.5502;
        tg_tidal_gauge_Eng = 'MUKHO';
    case '제주'
        lon_tidal_gauge = 126.5494;
        lat_tidal_gauge = 33.5364;
        tg_tidal_gauge_Eng = 'JEJU';
    case '속초'
        lon_tidal_gauge = 128.6741; % 128.5941
        lat_tidal_gauge = 38.2272; % 38.2072
        tg_tidal_gauge_Eng = 'SOCKCHO';
    case '포항'
        lon_tidal_gauge = 129.4231; % 129.3831
        lat_tidal_gauge = 36.0403; % 36.0403
        tg_tidal_gauge_Eng = 'POHANG';
    case '덕적도'
        lon_tidal_gauge = 126.1644;
        lat_tidal_gauge = 37.2311;
        tg_tidal_gauge_Eng = 'DEOKJEOKDO';
    case '안흥'
        lon_tidal_gauge = 126.0950;
        lat_tidal_gauge = 36.6622;
        tg_tidal_gauge_Eng = 'ANHEUNG';
    case '어청도'
        lon_tidal_gauge = 125.9825;
        lat_tidal_gauge = 36.1117;
        tg_tidal_gauge_Eng = 'EOCHEONGDO';
    case '위도'
        lon_tidal_gauge = 126.3078;
        lat_tidal_gauge = 35.6303;
        tg_tidal_gauge_Eng = 'WIDO';
    case '울릉도'
        lon_tidal_gauge = 130.9136;
        lat_tidal_gauge = 37.4913;
        tg_tidal_gauge_Eng = 'ULLEONGDO';
    case '흑산도'
        lon_tidal_gauge = 125.4219;
        lat_tidal_gauge = 34.7400;
        tg_tidal_gauge_Eng = 'HEUKSANDO';
    case '후포'
        lon_tidal_gauge = 129.4530;
        lat_tidal_gauge = 36.6775;
        tg_tidal_gauge_Eng = 'HUPO';
end

%% Load Storm Surge Height data
tgt_tc = 'MAEMI';
switch tgt_tc
    case 'MAYSAK'
        tgt_yr = '2020';
        s_mth = 8; s_day = 30;
        e_mth = 9; e_day = 6;
    case 'HAISHEN'
        tgt_yr = '2020'; 
        s_mth = 9; s_day = 1;
        e_mth = 9; e_day = 7;
    case 'MAEMI'
        tgt_yr = '2003';
        s_mth = 9; s_day = 7;
        e_mth = 9; e_day = 15;
end

sdate = datenum(datetime(str2double(tgt_yr),s_mth,s_day,0,0,0));
edate = datenum(datetime(str2double(tgt_yr),e_mth,e_day,0,0,0));

obsTimeSpan = sdate:1/24/60:edate;

%% Load sea level observation data
tpath = 'S:/home/user_006/03_DATA/SeaLevel_Obs';
cd(fullfile(tpath, tg_tidal_gauge))

ds_SL = readtable([tgt_yr '년 ' tg_tidal_gauge ' 조위관측소.txt']);

SL_fnames = fieldnames(ds_SL);
SL_date_1 = datenum(ds_SL.(SL_fnames{1}));
% SL_date_2 = datenum(ds_SL.(SL_fnames{2}));
SL_date = SL_date_1;% + SL_date_2;
SL_tide = ds_SL.(SL_fnames{2});

if all(~isnumeric(SL_tide))
    disp('No recorded sea level data!!')
end

%% Filter out redundant period
raw_SL_date = datenum(SL_date);
idx_SL_tgt_date = and(raw_SL_date > obsTimeSpan(1), raw_SL_date < obsTimeSpan(end));
tgt_SL_date = raw_SL_date(idx_SL_tgt_date);
tgt_SL_tide = SL_tide(idx_SL_tgt_date);
%{
for i = 1:length(tgt_SL_tide)
    tgt_SL_tide_new(i) = str2double(tgt_SL_tide{i});
end
tgt_SL_tide = tgt_SL_tide_new';
%}
clf; set(gcf, 'Position', [-1919 -4 1920 1003]);
plot(tgt_SL_date, tgt_SL_tide, 'k');

% tgt_SL_tide(tgt_SL_tide == 0) = NaN;
% tgt_SL_tide(tgt_SL_tide == 14) = NaN;

%% Get tidal info at the target location using TMD
cd 'C:/MATLIB/tmd_toolbox/'
[TMD_tide,ConList] = tmd_tide_pred('Model_Ind_2016',tgt_SL_date-9/24,lat_tidal_gauge,lon_tidal_gauge,'z');

%% Sizes are not identical from time to time
if ~isequal(size(TMD_tide),size(tgt_SL_tide))
    TMD_tide = TMD_tide';
end

%% Take out tidal signal from raw sea level info
valid_idx = ~isnan(tgt_SL_tide) & ~isnan(TMD_tide);
raw_tide = nan(size(tgt_SL_tide));
raw_tide(valid_idx) = tgt_SL_tide(valid_idx) - TMD_tide(valid_idx) * 100;
% raw_tide = tgt_SL_tide - TMD_tide*100;

clf; set(gcf, 'Position', [-1919 -4 1920 1003]); hold on;
plot(tgt_SL_date,tgt_SL_tide,'k','LineWidth',2);
plot(tgt_SL_date,TMD_tide*100,'b','LineWidth',2);
plot(tgt_SL_date,raw_tide,'r','LineWidth',2);

%% Apply FFT to raw_tide
Fs = 1000;
T = 1/Fs;
L = length(raw_tide);
t = (0:L-1)*T;

if any(isnan(raw_tide))
    raw_tide = fillmissing(raw_tide,'linear');
end

Y = fft(raw_tide);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

%% Check out which freq needs to be removed
clf; set(gcf, 'Position', [-1919 -4 1920 1003]); hold on;
plot(f, P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Take out the target noise
Fs = 1000;               
f_target = 1.36062;

wo = f_target / (Fs / 2);
bw = wo/2;

notchFilter = designfilt('bandstopiir', 'FilterOrder', 2, ...
                         'HalfPowerFrequency1', wo - bw, ...
                         'HalfPowerFrequency2', wo + bw, ...
                         'DesignMethod', 'butter');

filtered_tide = filtfilt(notchFilter, raw_tide);

clf; set(gcf, 'Position', [-1919 -4 1920 1003]); hold on;
plot(tgt_SL_date,raw_tide,'k','LineWidth',2);
plot(tgt_SL_date,filtered_tide,'r','LineWidth',2);

%% Moveing average to take out high frequency signal
prcsd_tide = movmean(filtered_tide,20,'omitnan');
prcsd_tide(isnan(tgt_SL_tide)) = NaN;
plot(tgt_SL_date,prcsd_tide,'b','LineWidth',2);
spath = 'S:/home/user_006/03_DATA/StormSurgeValue';
load(fullfile(spath, 'DATA.mat'))
% DATA.(tgt_tc).('Date') = tgt_SL_date;
DATA.(tgt_tc).(tg_tidal_gauge_Eng) = prcsd_tide;
save(fullfile(spath, 'DATA.mat'), 'DATA')







clc;clear all;close all;
%{
This script will go on a process that takes out tidal components from raw
sea level data recorded at each tidal gauge. Those tidal components consist of 4 major tides and
background tides. 
%}

addpath(genpath('S:/home/user_006/08_MATLIB'))
spath = 'C:/Users/admin/Desktop/OBS_TIDE/figures/';

%% Select target tidal gauge
tg_tidal_gauge = '제주';

switch tg_tidal_gauge
    case '마산'
        lon_tidal_gauge = 128.5838; % 128.5769 (TMD doesn't represent MASAN obs sites owing to low resolution)
        lat_tidal_gauge = 35.1883; % 35.1928 
    case '부산'
        lon_tidal_gauge = 129.0275;
        lat_tidal_gauge = 35.0942;
    case '서귀포'
        lon_tidal_gauge = 126.5689;
        lat_tidal_gauge = 33.2306;
    case '성산포'
        lon_tidal_gauge = 126.9275;
        lat_tidal_gauge = 33.4806;
    case '울산'
        lon_tidal_gauge = 129.4025;
        lat_tidal_gauge = 35.4994;
    case '제주'
        lon_tidal_gauge = 126.5494;
        lat_tidal_gauge = 33.5364;
    case '포항'
        lon_tidal_gauge = 129.3831;
        lat_tidal_gauge = 36.0403;
    case '덕적도'
        lon_tidal_gauge = 126.1644;
        lat_tidal_gauge = 37.2311;
    case '안흥'
        lon_tidal_gauge = 126.0950;
        lat_tidal_gauge = 36.6622;
    case '어청도'
        lon_tidal_gauge = 125.9825;
        lat_tidal_gauge = 36.1117;
    case '위도'
        lon_tidal_gauge = 126.3078;
        lat_tidal_gauge = 35.6303;
    case '흑산도'
        lon_tidal_gauge = 125.4219;
        lat_tidal_gauge = 34.7400;
end

%% Load Storm Surge Height data
tg_tc = 'HINNAMNOR'; % need to change in every iteration
tg_yr = '2022';

cd S:/home/user_006/01_WORK/240902_TY/03_RUN
tc_list = dir(['*' tg_tc]);
cd(tc_list.name)
fgs = grd_to_opnml('fort.14');
load('zeta.mat')
raw_zeta = read_adcirc_fort63('fort.63');
zeta = raw_zeta.zeta;

idx_tgt_fgs_x = and(fgs.x > lon_tidal_gauge - 0.1, fgs.x < lon_tidal_gauge + 0.1);
idx_tgt_fgs_y = and(fgs.y > lat_tidal_gauge - 0.1, fgs.y < lat_tidal_gauge + 0.1);
idx_tgt_fgs = and(idx_tgt_fgs_x,idx_tgt_fgs_y);

[~,num_iter] = size(zeta);
for i = 1:num_iter
    TS(i) = griddata(fgs.x(idx_tgt_fgs),fgs.y(idx_tgt_fgs),zeta(idx_tgt_fgs,i),lon_tidal_gauge,lat_tidal_gauge);
end

%% Take time info during the time TC alive
dlist = dir('*_PRESS');
cd(dlist(1).name)

flist = dir('*.dat');

raw_sctime = flist(1).name;
raw_ectime = flist(end).name;
sctime = datenum(datevec(raw_sctime(1:end-4),'YYYY-mm-dd_HH'));
ectime = datenum(datevec(raw_ectime(1:end-4),'YYYY-mm-dd_HH'));

cal_time = sctime:1/24/6:ectime;

%% Load sea level observation data
tpath = 'C:/Users/admin/Desktop/폭풍해일_관측자료/';
cd([tpath tg_tidal_gauge])

ds_SL = readtable([tg_yr '년 ' tg_tidal_gauge ' 조위관측소.txt']);

SL_fnames = fieldnames(ds_SL);
SL_date_1 = datenum(ds_SL.(SL_fnames{1}));
SL_date_2 = datenum(ds_SL.(SL_fnames{2}));
SL_date = SL_date_1 + SL_date_2;
SL_tide = ds_SL.(SL_fnames{3});

if all(~isnumeric(SL_tide))
    disp('No recorded sea level data!!')
end

%% Filter out redundant period
raw_SL_date = datenum(SL_date);
idx_SL_tgt_date = and(raw_SL_date > cal_time(1), raw_SL_date < cal_time(end));
tgt_SL_date = raw_SL_date(idx_SL_tgt_date);
tgt_SL_tide = SL_tide(idx_SL_tgt_date);

%% 
new_SL_date = tgt_SL_date(1):1/24/60:tgt_SL_date(end);
new_SL_tide = interp1(tgt_SL_date,tgt_SL_tide,new_SL_date,"linear");

if isequal(tg_tidal_gauge,'울산')
    new_SL_tide = new_SL_tide(1:14499);
    new_SL_date = new_SL_date(1:14499);
end

%% Get tidal info at the target location using TMD
cd 'C:/MATLIB/tmd_toolbox/'
[TMD_tide,ConList] = tmd_tide_pred('Model_Ind_2016',new_SL_date-9/24,lat_tidal_gauge,lon_tidal_gauge,'z');

%% Sizes are not identical from time to time
if ~isequal(size(TMD_tide),size(new_SL_tide))
    TMD_tide = TMD_tide';
end

%% Take out tidal signal from raw sea level info
raw_tide = new_SL_tide - TMD_tide*100;

figure;hold on;set(gcf,'Position',[1 1 1920 1003]);
plot(new_SL_date,new_SL_tide,'k','LineWidth',2);
plot(new_SL_date,TMD_tide*100,'b','LineWidth',2);
plot(new_SL_date,raw_tide,'r','LineWidth',2);

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
figure;hold on;set(gcf,'Position',[1 1 1920 1003]);
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

% f_target = 2.72124;
% filtered_tide = filtfilt(notchFilter, filtered_tide);

figure;hold on;set(gcf,'Position',[1 1 1920 1003]);
plot(new_SL_date,raw_tide,'k','LineWidth',2);
plot(new_SL_date,filtered_tide,'r','LineWidth',2);

%% Moveing average to take out high frequency signal
prcsd_tide = movmean(filtered_tide,40,'omitnan');
% for i = 1:2
%     prcsd_tide = movmean(prcsd_tide,40,'omitnan');
% end

%% Background tide
switch tg_tidal_gauge
    case '마산'
        bg_tide = 110;
    case '서귀포'
        bg_tide = 190;
    case '부산'
        bg_tide = 90;
    case '성산포'
        bg_tide = 160;
    case '울산'
        bg_tide = 50;
    case '제주'
        bg_tide = 180;
    case '포항'
        bg_tide = 50;
end

%% Plot
figure;hold on;set(gcf,'position',[1 1 1920 1003]);
plot(new_SL_date-9/24,prcsd_tide-bg_tide,'ko','MarkerFaceColor','k','MarkerSize',8);
plot(cal_time,TS*100,'r','LineWidth',5);
xlim([new_SL_date(1) new_SL_date(end)]);
ylim([-150 250]);
switch tg_tc
    case 'HINNAMNOR'
        sample_sdate = datenum(2022,8,31,0,0,0);
        sample_edate = datenum(2022,9,8,0,0,0);
        xtickVals = linspace(sample_sdate, sample_edate, 5);
        xlim([sample_sdate sample_edate]);
end

xticks(xtickVals);
xticklabels(datestr(xtickVals, 'mmm-dd'));
xlabel('Time (UTC, mm-dd)');
ylabel('Surge height (cm)');
grid on;
grid minor
set(gca, 'MinorGridLineWidth',1.5)
set(gcf, 'Color','w')
set(gca, 'FontSize', 35, 'FontWeight', 'bold', 'LineWidth', 4, 'Box', 'on');
lg = legend([{'OBS.'},{'JMA-MSM'}],'Location','SouthWest','FontSize',35);
leg_pos = get(lg,'position');
set(lg,'position',[leg_pos(1),leg_pos(2),leg_pos(3)*1.1,leg_pos(4)]);

print('-vector',strcat(spath,tg_tc,'/',tg_tc,'_',tg_tidal_gauge,'.png'), '-dpng','-r300');



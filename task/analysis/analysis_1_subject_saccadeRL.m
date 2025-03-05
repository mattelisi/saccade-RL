% anti-saccade task; analysis single subject
clear all

addpath('../functions');
addpath('./analysis_functions');


%% screen settings
scr.subDist = 80;   % subject distance (cm)
scr.width   = 570;  % monitor width (mm)

scr.xres = 1920;
scr.yres = 1080;
scr.xCenter = scr.xres/2;
scr.yCenter = scr.yres/2 ;
ppd = va2pix(1,scr); % pixel per degree

%% data table
dtab = readtable('../data/S3');
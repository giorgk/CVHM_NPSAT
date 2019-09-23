
% Set time span
startTime = [1995 7];
endTime = [1999 12];

opt.prefix = 'cvhm';
opt.simFolder = ['Sim_' num2str(startTime(1)) 'm' num2str(startTime(2)) '_'  num2str(endTime(1)) 'm' num2str(endTime(2))];
opt.timestring = [num2str(startTime(1)) 'm' num2str(startTime(2)) '_'  num2str(endTime(1)) 'm' num2str(endTime(2))];
%opt.cbcf_path = '../CVHM_DATA/ClaudiaRun/';
opt.cbcf_path = '/media/giorgk/FIONA/UBU_BACKUP/Documents/UCDAVIS/CVHM_DATA/ClaudiaRun/';
opt.gis_path = 'gis_data/';
opt.std_Htol = 0.5;
opt.LU = 'LU_2000';
opt.do_plot = false;
%% Do NOT run this at once. Just line by line by selecting and pressing F9
[sumSTM, sumWLL, sumRCH, totDays, opt] = CVHMInputRCH(startTime, endTime, opt);

CVHMInputHEAD(startTime, endTime, opt);

avSTRM = sumSTM/totDays;
CVHMInputStreamA(avSTRM, opt);
% Run the CVHM_HOU.hipnc animation to generate the streams
% Then run the next line to write the input file
CVHMInputStreamB(opt);

avWells = sumWLL/totDays;
Wells_Gen = CVHMgenerateWells(avWells, opt);
save([opt.simFolder filesep opt.timestring '_DATA'], 'sumSTM', 'sumWLL', 'sumRCH', 'totDays', 'Wells_Gen');

botFile = 'CVHM_Bot_elev.npsat';
topFile = [opt.simFolder filesep opt.prefix '_' opt.timestring  '_Top_iter_1.npsat'];
opt.above_thres = 10;
opt.below_thres = 30;
opt.minScreen = 10;
opt.iter = 1;
Wells_Gen_mod = writeWells(Wells_Gen, botFile,topFile, opt);
save([opt.simFolder filesep opt.timestring '_DATA'], 'Wells_Gen_mod', '-append');

%% This snippet is used to adjust manually the water table elevation.
% Ideally we should let the npsat_engine do that. However Central Valley
% is a special case
clear
clc
% set the previous iteration
iter = 2;

% First read the water table
wtcfile = [opt.simFolder filesep 'output' filesep opt.prefix '_' opt.timestring '_H' num2str(opt.std_Htol) '_top_002_'];
WTC = ReadWaterTableCloudPoints(wtcfile,1,true, iter);
% If there are extreme wtc points fix them here
WTC_mod = WTC;
% create an interpolant
FnewElev = scatteredInterpolant(WTC_mod(:,1), WTC_mod(:,2), WTC_mod(:,3));
% Read old top
TopOld = read_Scattered([opt.simFolder filesep opt.prefix '_' opt.timestring  '_Top_iter_1.npsat'], 2);
% interpolate the new top elevation
new_xy_top = [TopOld.p(:,1) TopOld.p(:,2) FnewElev(TopOld.p(:,1), TopOld.p(:,2))];
topnewfile = [opt.simFolder filesep opt.prefix '_' opt.timestring  '_Top_iter_' num2str(iter+1) '.npsat'];
writeScatteredData(topnewfile, struct('PDIM',2,'TYPE','HOR','MODE','SIMPLE'), new_xy_top);
% load initial generated well data
load([opt.simFolder filesep opt.timestring '_DATA'],'Wells_Gen');
botFile = 'CVHM_Bot_elev.npsat';
opt.iter = iter+1;
Wells_Gen_mod = writeWells(Wells_Gen, botFile, topnewfile, opt);


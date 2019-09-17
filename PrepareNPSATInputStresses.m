% Set time span
startTime = [1995 7];
endTime = [1999 12];

opt.prefix = 'cvhm';
opt.simFolder = ['Sim_' num2str(startTime(1)) 'm' num2str(startTime(2)) '_'  num2str(endTime(1)) 'm' num2str(endTime(2))];
opt.timestring = [num2str(startTime(1)) 'm' num2str(startTime(2)) '_'  num2str(endTime(1)) 'm' num2str(endTime(2))];
opt.cbcf_path = '../CVHM_DATA/ClaudiaRun/';
opt.gis_path = 'gis_data/';
opt.std_Htol = 5;
opt.LU = 'LU_2000';
opt.do_plot = false;

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
topFile = [opt.simFolder filesep opt.prefix '_' opt.timestring  '_Top.npsat'];
opt.above_thres = 10;
opt.below_thres = 30;
opt.minScreen = 10;
opt.wellVersion = 1;
Wells_Gen_mod = writeWells(Wells_Gen, botFile,topFile, opt);
save([opt.simFolder filesep opt.timestring '_DATA'], 'Wells_Gen_mod', '-append');
% Define data paths
c2vsim_path = '/home/giorgk/Documents/UCDAVIS/C2VSIM_DATA/';
%% load the groundwater budget file of the C2Vsim simulation
load([c2vsim_path 'BaseRun.mat']);
%% Extract the storage
m = 10; y = 1921;
STRG = zeros(1056,1);
for ii = 1:21
    STRG = STRG + out.gw(ii).Data(:,2);
end
%% prepare the cumulative change in storage from water year 1961 up to 2009
id_oct_1961 = 481;
monthlySTRG = STRG(id_oct_1961:end)-STRG(id_oct_1961);
yr = reshape(repmat(1961:2009,12,1),12*49,1);
yr = yr(10:end-3);
plot(monthlySTRG/10^6)
%%
id_1961 = 41;
yearlySTRG = mean(reshape(STRG, 12, 88))';
yearlySTRG = yearlySTRG(id_1961:end) - yearlySTRG(id_1961);
plot(yearlySTRG/10^6)

%% print as JS variables
writeTimeSeries2JS('Jscripts/Storage', 'C2VsimMonthly', ones(length(monthlySTRG),1), repmat([10:12 1:9]',576/12,1), yr, monthlySTRG/10^6, true);
writeTimeSeries2JS('Jscripts/Storage', 'C2VsimYearly', ones(length(yearlySTRG),1), 2*ones(length(yearlySTRG),1),[1962:2009]', yearlySTRG/10^6, true);
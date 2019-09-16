function [sumSTM, sumWLL, totDays, opt] = CVHMInputRCH(startTime, endTime, opt)
% startTime, endTime are [1x2] matrices [y m]

if ~isfield(opt,'cbcf_path')
    opt.cbcf_path = '../CVHM_DATA/ClaudiaRun/';
end
if ~isfield(opt,'gis_path')
    opt.gis_path = 'gis_data/';
end

if ~isfield(opt, 'simFolder')
    opt.simFolder = ['Sim_' num2str(startTime(1)) 'm' num2str(startTime(2)) '_'  num2str(endTime(1)) 'm' num2str(endTime(2))];
end
if ~exist(opt.simFolder,'dir')
    mkdir(opt.simFolder)
end

% ======== load Data
load([opt.cbcf_path 'CBCdaily.mat'])
load([opt.gis_path 'BAS_active_shp.mat']);
load('BufferPnts.mat')
% make unique list of buffer outline points
buff_pnt = [buff(1,1).X(1) buff(1,1).Y(1)];
for ii = 1:length(buff(1,1).X) - 1
    dst = sqrt((buff_pnt(:,1) - buff(1,1).X(ii)).^2 + (buff_pnt(:,2) - buff(1,1).Y(ii)).^2);
    if min(dst) > 0.1
        buff_pnt = [buff_pnt; buff(1,1).X(ii) buff(1,1).Y(ii)];
    end
end


% =========Calculate water balance for the selected period
% find indices for startTime and endTime
istart = find(CBCdaily.ym(:,1) == startTime(1) & CBCdaily.ym(:,2) == startTime(2));
iend = find(CBCdaily.ym(:,1) == endTime(1) & CBCdaily.ym(:,2) == endTime(2));
timestring = [num2str(startTime(1)) 'm' num2str(startTime(2)) '_'  num2str(endTime(1)) 'm' num2str(endTime(2))];

% Find the sum of stresses
sumRCH = zeros(441,98);
sumSTM = zeros(441,98);
sumWLL = zeros(441,98);
totDays = 0;
for ii = istart:iend
    ndays = eomday(CBCdaily.ym(ii,1), CBCdaily.ym(ii,2));
    sumRCH = sumRCH + CBCdaily.RCH{ii,1}*ndays;
    ind_str = sub2ind([441 98], str{ii,1}(:,1), str{ii,1}(:,2));
    sumSTM(ind_str) = sumSTM(ind_str) + str{ii,1}(:,3)*ndays;
    ind_wll = sub2ind([441 98],mnw{ii,1}(:,1), mnw{ii,1}(:,2));
    sumWLL(ind_wll) = sumWLL(ind_wll) + mnw{ii,1}(:,3)*ndays;
    ind_wll = sub2ind([441 98],frw{ii,1}(:,1), frw{ii,1}(:,2));
    sumWLL(ind_wll) = sumWLL(ind_wll) + frw{ii,1}(:,3)*ndays;
    totDays = totDays + ndays;
end
err = sum(sum(sumRCH)) + sum(sum(sumSTM)) + sum(sum(sumWLL));
posFlows = sum(sumRCH(sumRCH>0)) + sum(sumSTM(sumSTM>0)) + sum(sumWLL(sumWLL>0));
negFlows = sum(sumRCH(sumRCH<0)) + sum(sumSTM(sumSTM<0)) + sum(sumWLL(sumWLL<0));
err_ratio = posFlows/abs(negFlows); 
display(['Inflows: ' num2str(posFlows/totDays) ', Outflows: ' num2str(negFlows/totDays)]);
display(['Ratio: ' num2str(err_ratio)]);

% Decrease the negative values so that the water balance is zero
if err_ratio < 1
    sumRCH(sumRCH<0) = sumRCH(sumRCH<0)*err_ratio;
    sumSTM(sumSTM<0) = sumSTM(sumSTM<0)*err_ratio;
    sumWLL(sumWLL<0) = sumWLL(sumWLL<0)*err_ratio;
else
    sumRCH(sumRCH>0) = sumRCH(sumRCH>0)/err_ratio;
    sumSTM(sumSTM>0) = sumSTM(sumSTM>0)/err_ratio;
    sumWLL(sumWLL>0) = sumWLL(sumWLL>0)/err_ratio;
end

% ================= Write Recharge files
avRCH = sumRCH/totDays;

R = [bas_active.ROW]';
C = [bas_active.COLUMN_]';
xy_rch = [];
for ii = 1:size(avRCH,1)
    for jj = 1:size(avRCH,2)
        id = find(R == ii & C == jj);
        if isempty(id)
            continue;
        elseif length(id) == 1
            xc = mean(bas_active(id,1).X(1:end-1));
            yc = mean(bas_active(id,1).Y(1:end-1));
            xy_rch = [xy_rch; xc yc avRCH(ii,jj)];
        else
            error('Thats so wrong');
        end
    end
end
% convert volume to rate
xy_rch(:,3) = xy_rch(:,3)./(1609.344*1609.344);
% create an interpolant with the existing data
Frch = scatteredInterpolant(xy_rch(:,1), xy_rch(:,2), xy_rch(:,3));
Frch.Method = 'nearest';
Frch.ExtrapolationMethod = 'nearest';
% Extrapolate values on the buffer nodes
buf_val = Frch(buff_pnt(:,1), buff_pnt(:,2));
xy_rch = [xy_rch; buff_pnt buf_val];
% write the rch file
writeScatteredData([opt.simFolder filesep opt.prefix '_' timestring  '_Rch.npsat'], struct('PDIM',2,'TYPE','HOR','MODE','SIMPLE'), xy_rch);


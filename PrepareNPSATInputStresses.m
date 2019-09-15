function [sumSTM, sumWLL, totDays] = PrepareNPSATInputStresses(startTime, endTime, prefix, opt)
% startTime, endTime are [1x2] matrices [y m]

if ~isfield(opt,'cbcf_path')
    opt.cbcf_path = '../CVHM_DATA/ClaudiaRun/';
end
if ~isfield(opt,'gis_path')
    opt.gis_path = 'gis_data/';
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

% find indices for startTime and endTime
istart = find(CBCdaily.ym(:,1) == startTime(1) & CBCdaily.ym(:,2) == startTime(2));
iend = find(CBCdaily.ym(:,1) == endTime(1) & CBCdaily.ym(:,2) == endTime(2));
timestring = [num2str(startTime(1)) 'm' num2str(startTime(2)) '_'  num2str(endTime(1)) 'm' num2str(endTime(2))]

% Find the sum of stresses
sumRCH = zeros(441,98);
sumSTM = zeros(441,98);
sumWLL = zeros(441,98);
totDays = 0;
for ii = istart:iend
    ndays = eomday(CBCdaily.ym(ii,1), CBCdaily.ym(ii,2));
    sumRCH = sumRCH + CBCdaily.RCH{ii,1}*ndays;
    ind_str = sub2ind([441 98],str{1,1}(:,1), str{1,1}(:,2));
    sumSTM(ind_str) = sumSTM(ind_str) + str{1,1}(:,3)*ndays;
    ind_wll = sub2ind([441 98],mnw{1,1}(:,1), mnw{1,1}(:,2));
    sumWLL(ind_wll) = sumWLL(ind_wll) + mnw{1,1}(:,3)*ndays;
    ind_wll = sub2ind([441 98],frw{1,1}(:,1), frw{1,1}(:,2));
    sumWLL(ind_wll) = sumWLL(ind_wll) + frw{1,1}(:,3)*ndays;
    totDays = totDays + ndays;
end
err = sum(sum(sumRCH)) + sum(sum(sumSTM)) + sum(sum(sumWLL));
posFlows = sum(sumRCH(sumRCH>0)) + sum(sumSTM(sumSTM>0)) + sum(sumWLL(sumWLL>0));
negFlows = sum(sumRCH(sumRCH<0)) + sum(sumSTM(sumSTM<0)) + sum(sumWLL(sumWLL<0));
err_ratio = posFlows/abs(negFlows); 

% Decrease the negative values so that the water balance is zero
sumRCH(sumRCH<0) = sumRCH(sumRCH<0)*err_ratio;
sumSTM(sumSTM<0) = sumSTM(sumSTM<0)*err_ratio;
sumWLL(sumWLL<0) = sumWLL(sumWLL<0)*err_ratio;

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
writeScatteredData([prefix '_' timestring  'rch.npsat'], struct('PDIM',2,'TYPE','HOR','MODE','SIMPLE'), xy_rch);

% ================== Write Hydraulic head ============== 
load('HeadAnalysisData.mat')
load([opt.cbcf_path 'HEADS.mat'])
%------ loop through the edges of the bas and set the ij ids for those that 
% match the outline segments
LNS_IJ = nan(size(PNTS,1),2);
for ii = 1:size(mesh,1)
    %[Xs, Ys] = polysplit(mesh(ii,1).X, mesh(ii,1).Y);
    Xs{1,1} = mesh(ii,1).X; Ys{1,1} = mesh(ii,1).Y;
    for jj = 1:size(Xs,1)
        for k = 1:length(Xs{jj,1})-1
            pa = [Xs{jj,1}(k) Ys{jj,1}(k)];
            dst = sqrt((pa(1) - PNTS(:, 1)).^2 + (pa(2) - PNTS(:, 2)).^2);
            ida = find(dst < 0.01);
            if ~isempty(ida)
                pb = [Xs{jj,1}(k+1) Ys{jj,1}(k+1)];
                dst = sqrt((pb(1) - PNTS(:, 1)).^2 + (pb(2) - PNTS(:, 2)).^2);
                idb = find(dst < 0.01);
                if ~isempty(idb)
                    for kk = 1:length(ida)
                        for kkk = 1:length(idb)
                            idmn = min([ida(kk) idb(kkk)]);
                            idmx = max([ida(kk) idb(kkk)]);
                            if idmx - idmn == 1
                                LNS_IJ(idmn,:) = [mesh(ii,1).R mesh(ii,1).C];
                            end
                        end
                    end
                end
            end
        end
    end
end

% For each time step compute the top most head value
clear TOPHEAD toplay
cnt = 1;
for it = istart:iend
    it
    TOPHEAD{cnt,1} = nan(size(HEADS{1,1},1), size(HEADS{1,1},2));
    toplay{cnt,1} = nan(size(HEADS{1,1},1), size(HEADS{1,1},2));
    for ii = 1:size(HEADS{1,1},1)
        for jj = 1:size(HEADS{1,1},2)
            for kk = 1:size(HEADS{1,1},3)%:-1:1
                if HEADS{it,1}(ii,jj,kk) ~= HEADS{1,1}(1,1,1)
                    TOPHEAD{cnt,1}(ii,jj) = HEADS{it,1}(ii,jj,kk);
                    toplay{cnt,1}(ii, jj) = kk;
                    break
                end
            end
        end
    end
end
%% Calculate average head and Standard deviation
head_std = nan(size(HEADS{1,1},1), size(HEADS{1,1},2));
head_Av = nan(size(HEADS{1,1},1), size(HEADS{1,1},2));
timespan = istart2:1:iend2;
for ii = 1:size(HEADS{1,1},1)
    for jj = 1:size(HEADS{1,1},2)
        temp = nan(length(timespan),1);
        for kk = 1:length(timespan)
            temp(kk,1) = TOPHEAD{timespan(kk),1}(ii,jj);
        end
        head_std(ii,jj) = std(temp);
        head_Av(ii,jj) = mean(temp);
    end
end

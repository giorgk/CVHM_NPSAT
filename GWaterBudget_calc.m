% path where the cbc flow data in matlab format are located
% Claudia's Caped Run
%
%% PATHS
% D218 path = /home/UCDAVIS/CVHM_NPSAT/ClaudiaRun
%% Home Linux
cbcf_path = '/home/giorgk/Documents/UCDAVIS/CVHM_DATA/ClaudiaRun/';
%% UCDAVIS Linux
% cbcf_path = '/media/giorgk/DATA/giorgk/Documents/CVHM_DATA/ClaudiaRun/';
%% Get Cell by cell flow data.
%{
m=4;y=1961;
t2=0;
for ii=1:510
    if m==13
        y
        m=1;y=y+1;
    end
    cnst=floor((ii-1)/10);
    if t2 < cnst*10+1
        clear OUT
        load([cbcf_path 'cbcf_' num2str(cnst*10+1) '_' num2str(10*cnst+10) '.mat']);
        t2 = cnst*10+1;
    end
    CBCdaily.ym(ii,:)=[y m];
    CBCdaily.STRG{ii,1} = sum(OUT{ii,1}{1,2},3); %Storage (We loose the vertical distribution but not the spatial)
    CBCdaily.IBST{ii,1} = sum(OUT{ii,1}{7,2},3); %Inst. IB Storage
    if ii == 1
        CBCdaily.STRM.IJ = OUT{ii,1}{8,2}(:,1:3);%Stream Leackage
        CBCdaily.STRM.data = zeros(size(OUT{ii,1}{8,2},1),510);
    end
    CBCdaily.STRM.data(:,ii) = OUT{ii,1}{8,2}(:,4);
    CBCdaily.RCH{ii,1} = reshape(OUT{ii,1}{11,2}(:,2),98,441)';
    CBCdaily.WELLS.MNW{ii,1} = OUT{ii,1}{9,2};
    CBCdaily.WELLS.FRM{ii,1} = OUT{ii,1}{10,2};
    m=m+1;
end
%}
%% save cell by cell values
%save([cbcf_path 'CBCdaily.mat'],'CBCdaily')
%% load instead of running the above snippets
load([cbcf_path 'CBCdaily.mat'])
%% load the cell shapefile
% keep only the usefull fields
bas = shaperead(['gis_data' filesep 'BAS_active.shp']);
bas_rc = [[bas.ROW]' [bas.COLUMN_]'];
bud = rmfield(bas, {'lay1_act','lay2_act','lay3_act','lay45_act','lay6_act','lay7_act','lay8_act','lay9_act','lay10_act',...
        'strt_hd_01','strt_hd_02','strt_hd_03','strt_hd_04','strt_hd_05','strt_hd_06','strt_hd_07','strt_hd_08','strt_hd_09','strt_hd_10'});
%% add Subbasins and farm id
farms = shaperead(['gis_data' filesep 'FARMS_poly.shp']);
XYmean = zeros(length(bud),2);
for ii = 1:length(bud)
   XYmean(ii,:) = [mean(bud(ii,1).X(1:4)) mean(bud(ii,1).Y(1:4))]; 
end
for ii = 1:length(farms)
   id_in = find(inpolygon(XYmean(:,1), XYmean(:,2), farms(ii,1).X, farms(ii,1).Y));
   for jj = 1:length(id_in)
       bud(id_in(jj),1).farm_id = farms(ii,1).dwr_sbrgns;
       bud(id_in(jj),1).bas_id = farms(ii,1).Basins;
   end
end
%% Calculate monthly regional storage time series
% CV(1) Basins(2-4) Farms(5-25)
% To get monthly values we have to multiply the daily values with the number of days per month 
RR = [bud.ROW]';
CC = [bud.COLUMN_]';
B_ID = [bud.bas_id]';
F_ID = [bud.farm_id]';
basin_ids = unique(B_ID)';
farm_ids = unique(F_ID)';
%}
%% preprocess streams and wells so that you have unique rows and columns
%{
clear str mnw frw
for ii = 1:510
    ii
    % STREAMS
    rc = CBCdaily.STRM.IJ(:,2:3);
    str{ii,1} = unique(rc,'rows');
    str{ii,1} = [str{ii,1} zeros(size(str{ii,1},1),1)];
    for jj = 1:size(str{ii,1},1)
        id = find(rc(:,1) == str{ii,1}(jj,1) & rc(:,2) == str{ii,1}(jj,2));
        str{ii,1}(jj,3) = sum(CBCdaily.STRM.data(id,ii));
    end
    
    % MNW
    rc = CBCdaily.WELLS.MNW{ii,1}(:,2:3);
    mnw{ii,1} = unique(rc,'rows');
    mnw{ii,1} = [mnw{ii,1} zeros(size(mnw{ii,1},1),1)];
    for jj = 1:size(mnw{ii,1},1)
        id = find(rc(:,1) == mnw{ii,1}(jj,1) & rc(:,2) == mnw{ii,1}(jj,2));
        mnw{ii,1}(jj,3) = sum(CBCdaily.WELLS.MNW{ii,1}(id,4));
    end
    
    % FRM
    rc = CBCdaily.WELLS.FRM{ii,1}(:,2:3);
    frw{ii,1} = unique(rc,'rows');
    frw{ii,1} = [frw{ii,1} zeros(size(frw{ii,1},1),1)];
    for jj = 1:size(frw{ii,1},1)
        id = find(rc(:,1) == frw{ii,1}(jj,1) & rc(:,2) == frw{ii,1}(jj,2));
        frw{ii,1}(jj,3) = sum(CBCdaily.WELLS.FRM{ii,1}(id,4));
    end
end
%%
save([cbcf_path 'CBCdaily.mat'],'str','mnw','frw','-append')
%}
%%
%
clear StorageTimeSeries DiscrepTimeSeries
StorageTimeSeries = zeros(510, 25);
% DiscrepTimeSeries contains the water balance that will actually used
% during averaging and it is calculated as RCH + STRM - Wells. In practice
% this should be very close to storage
DiscrepTimeSeries = zeros(510, 25);
for ii = 1:510
    ii
    ndays = eomday(CBCdaily.ym(ii,1), CBCdaily.ym(ii,2));
    % CV
    StorageTimeSeries(ii,1) = (sum(sum(CBCdaily.STRG{ii,1}))+sum(sum(CBCdaily.IBST{ii,1})))*ndays;
    DiscrepTimeSeries(ii,1) = (sum(sum(CBCdaily.RCH{ii,1})) + sum(CBCdaily.STRM.data(:,ii)) + ...
        sum(CBCdaily.WELLS.MNW{ii,1}(:,4)) + sum(CBCdaily.WELLS.FRM{ii,1}(:,4)))*ndays;
    
    % Basins
    for jj = basin_ids
        cell_id = find(B_ID == jj);
        r = RR(cell_id,1);
        c = CC(cell_id,1);
        ind = sub2ind([441 98],r,c);
        StorageTimeSeries(ii,2 +jj) = (sum(CBCdaily.STRG{ii,1}(ind))+sum(CBCdaily.IBST{ii,1}(ind)))*ndays;
        % find the streams and wells that belong to the basin
        [~,~,ib_strm] = intersect([r c], str{ii,1}(:,1:2),'rows');
        [~,~,ib_mnw] = intersect([r c], mnw{ii,1}(:,1:2),'rows');
        [~,~,ib_frm] = intersect([r c], frw{ii,1}(:,1:2),'rows');
        
        DiscrepTimeSeries(ii,2 +jj) = (sum(CBCdaily.RCH{ii,1}(ind)) + sum(str{ii,1}(ib_strm,3)) + ...
            sum(mnw{ii,1}(ib_mnw,3)) + sum(frw{ii,1}(ib_frm,3)))*ndays;
    end
    
    % Farms
    for jj = farm_ids
        cell_id = find(F_ID == jj);
        r = RR(cell_id,1);
        c = CC(cell_id,1);
        ind = sub2ind([441 98],r,c);
        StorageTimeSeries(ii,4 +jj) = (sum(CBCdaily.STRG{ii,1}(ind))+sum(CBCdaily.IBST{ii,1}(ind)))*ndays;
        % find the streams and wells that belong to the farm
        [~,~,ib_strm] = intersect([r c], str{ii,1}(:,1:2),'rows');
        [~,~,ib_mnw] = intersect([r c], mnw{ii,1}(:,1:2),'rows');
        [~,~,ib_frm] = intersect([r c], frw{ii,1}(:,1:2),'rows');
        
        DiscrepTimeSeries(ii,4 +jj) = (sum(CBCdaily.RCH{ii,1}(ind)) + sum(str{ii,1}(ib_strm,3)) + ...
            sum(mnw{ii,1}(ib_mnw,3)) + sum(frw{ii,1}(ib_frm,3)))*ndays;
    end
end
%% Do some plots and comparisons
k = 2;
clf
plot(cumsum([0; -StorageTimeSeries(7:end,k)]/10^6/1233.48184),'b')
hold on
plot(cumsum([0; DiscrepTimeSeries(7:end,k)]/10^6/1233.48184),'r')
%% write CV and basins to js files
writeTimeSeries2JS('Jscripts/StorageBasins', 'cvhmCVMonthly', ones(504,1), CBCdaily.ym(7:end,2), CBCdaily.ym(7:end,1), cumsum([0; -StorageTimeSeries(7:end,1)]/10^6/1233.48184), false);
writeTimeSeries2JS('Jscripts/StorageBasins', 'cvhmTLBMonthly', ones(504,1), CBCdaily.ym(7:end,2), CBCdaily.ym(7:end,1), cumsum([0; -StorageTimeSeries(7:end,2)]/10^6/1233.48184), true);
writeTimeSeries2JS('Jscripts/StorageBasins', 'cvhmSJVMonthly', ones(504,1), CBCdaily.ym(7:end,2), CBCdaily.ym(7:end,1), cumsum([0; -StorageTimeSeries(7:end,3)]/10^6/1233.48184), true);
writeTimeSeries2JS('Jscripts/StorageBasins', 'cvhmSACMonthly', ones(504,1), CBCdaily.ym(7:end,2), CBCdaily.ym(7:end,1), cumsum([0; -StorageTimeSeries(7:end,4)]/10^6/1233.48184), true);
%% Write Farms 
append = false;
for ii = 1:21
    if ii > 1; append = true;end
    writeTimeSeries2JS('Jscripts/StorageFarms', ['Farm_' num2str(ii)], ones(504,1), CBCdaily.ym(7:end,2), CBCdaily.ym(7:end,1), cumsum([0; -StorageTimeSeries(7:end,4+ii)]/10^6/1233.48184), append);
end
%}
%% Find periods with minimum storage change
Descrepancy = nan(510*510,4);
Nyears = nan(510*510,1);
TimeSpan = nan(510*510,2);
cnt = 1;
TLB_cells = sub2ind([441 98], RR(B_ID == 0,1), CC(B_ID == 0,1));
SJV_cells = sub2ind([441 98], RR(B_ID == 1,1), CC(B_ID == 1,1));
SAC_cells = sub2ind([441 98], RR(B_ID == 2,1), CC(B_ID == 2,1));

for sy = 7:1:510-11
    sy
    for ey = sy+11:1:510
        temp_rch = zeros(441,98);
        temp_well = zeros(441,98);
        temp_strm = zeros(441,98);
        
        for k = sy:ey
            ndays =eomday(CBCdaily.ym(k,1), CBCdaily.ym(k,2));
            % recharge
            temp_rch = temp_rch + CBCdaily.RCH{k,1}*ndays;
            % streams
            ind_str = sub2ind([441 98],str{1,1}(:,1), str{1,1}(:,2));
            temp_strm(ind_str) = temp_strm(ind_str) + str{1,1}(:,3)*ndays;
            % MN wells
            ind_wll = sub2ind([441 98],mnw{1,1}(:,1), mnw{1,1}(:,2));
            temp_well(ind_wll) = temp_well(ind_wll) + mnw{1,1}(:,3)*ndays;
            % Farm wells
            ind_wll = sub2ind([441 98],frw{1,1}(:,1), frw{1,1}(:,2));
            temp_well(ind_wll) = temp_well(ind_wll) + frw{1,1}(:,3)*ndays;
        end
        Descrepancy(cnt , 1) = (sum(sum(temp_rch)) + sum(sum(temp_strm)) + sum(sum(temp_well)))/10^6/1233.48184;
        Descrepancy(cnt , 2) = (sum(sum(temp_rch(TLB_cells))) + sum(sum(temp_strm(TLB_cells))) + sum(sum(temp_well(TLB_cells))))/10^6/1233.48184;
        Descrepancy(cnt , 3) = (sum(sum(temp_rch(SJV_cells))) + sum(sum(temp_strm(SJV_cells))) + sum(sum(temp_well(SJV_cells))))/10^6/1233.48184;
        Descrepancy(cnt , 4) = (sum(sum(temp_rch(SAC_cells))) + sum(sum(temp_strm(SAC_cells))) + sum(sum(temp_well(SAC_cells))))/10^6/1233.48184;
        Nyears(cnt,1) = length(sy:ey)/12;
        TimeSpan(cnt,:) = [sy ey];
        cnt = cnt + 1;
    end
end
Descrepancy(cnt:end,:) = [];
Nyears(cnt:end,:) = [];
TimeSpan(cnt:end,:) = [];
%% isolate selected periods e.g. the last decade
ind_sy = 391;
temp_inds = TimeSpan(:,1) >= ind_sy;
temp_ts = TimeSpan(temp_inds,:);
temp_Ny = Nyears(temp_inds,:);
temp_dp = Descrepancy(temp_inds,:);
%% find out the time spans for the
plot(temp_Ny, abs(temp_dp(:,1)),'.')
plot(temp_Ny, sum(abs(temp_dp(:,2:4)),2),'.')
%% write the selected points as javascript variable for scatter plot
nmfile = 'Jscripts/TimeSpanDiscrepancyBasins';
varname = 'Basinsdiscrep';
fid = fopen([nmfile '.js'],'w');
fprintf(fid, 'var %s = [\n', varname);
for ii = 1:length(temp_Ny)
    fprintf(fid, '{ x: %0.5f, y: %0.5f, sy: new Date(%d,%d,%d), ey: new Date(%d,%d,%d) },\n', [temp_Ny(ii) sum(abs(temp_dp(ii,2:4)),2) CBCdaily.ym(temp_ts(ii,1),:) 1 CBCdaily.ym(temp_ts(ii,2),:) 1]);
end
fprintf(fid, '];');
fclose(fid);
%%
totDays = 0;
for ii = 1:510
    ndays = eomday(RCH.ym(ii,1), RCH.ym(ii,2));
    totDays = totDays + ndays;
    R(ii,1) = sum(sum(RCH.data{ii,1}.*ndays));
    S(ii,1) = sum(STRLK.data(:,ii).*ndays);
    W(ii,1) = sum(WELLS.data{ii,1}(:,4).*ndays) + sum(WELLS.data{ii,2}(:,4).*ndays);
    Diff(ii,1) = R(ii,1) + S(ii,1) + W(ii,1);
end
%%
plot(cumsum(R),'r')
hold on
plot(cumsum(S),'g')
plot(cumsum(W),'b')
%
% From file https://water.usgs.gov/ogw/modflow/mf6io.pdf pg21:
% The cell-by-cell storage term gives the net flow to or from storage in a variable-head cell. The net storage
% for each cell in the grid is saved in transient simulations if the appropriate flags are set. Withdrawal from storage
% in the cell is considered positive, whereas accumulation in storage is considered negative.
%% CVHM default run
% ds218 : /home/UCDAVIS/CVHM_NPSAT
load('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/CBC.mat');

%%
m=4;y=1961;
basin_ids = unique([bud.bas_id]')';
farm_ids = unique([bud.farm_id]')';
for ii=1:510
    if m==13
        m=1;y=y+1;
    end
    dailySUM(ii,1)=y;
    dailySUM(ii,2)=m;
    
    % All CV
    dailySUM(ii,3)=sum(sum(sum(CBC{ii,1}{1,2}))) + sum(sum(sum(CBC{ii,1}{7,2})));
    
    % group by Basin
    for jj = basin_ids
        cell_id = find([bud.bas_id] == jj);
        r = [bud(cell_id,1).ROW]';
        c = [bud(cell_id,1).COLUMN_]';
        tmp_strg = 0;
        for k = 1:10
            tmp = CBC{ii,1}{1,2}(:,:,k) + CBC{ii,1}{7,2}(:,:,k);
            tmp_strg = tmp_strg + sum(tmp(sub2ind(size(tmp),r,c)));
        end
        dailySUM(ii,4+jj) = tmp_strg;
    end
    
    % group by farm
    for jj = farm_ids
        cell_id = find([bud.farm_id] == jj);
        r = [bud(cell_id,1).ROW]';
        c = [bud(cell_id,1).COLUMN_]';
        tmp_strg = 0;
        for k = 1:10
            tmp = CBC{ii,1}{1,2}(:,:,k) + CBC{ii,1}{7,2}(:,:,k);
            tmp_strg = tmp_strg + sum(tmp(sub2ind(size(tmp),r,c)));
        end
        dailySUM(ii,6+jj) = tmp_strg;
    end
    m=m+1;
end
%%
% Test the above code
% All lines should be identical
plot(-cumsum([0;dailySUM(:,3)]),'b')
hold on
plot(-cumsum([0;sum(dailySUM(:,4:6),2)]),'r')
plot(-cumsum([0;sum(dailySUM(:,7:27),2)]),'g')
%
%% Multiply the values by the number of days and average per year
cnt=1;
yc=1;
clear monthlySUM YearlySum
for ii=7:510
    ndays = eomday(dailySUM(ii,1), dailySUM(ii,2));
    monthlySUM(cnt,:) = [dailySUM(ii,1:2) dailySUM(ii,3:end)*ndays];
    if mod(cnt,12) == 0
        YearlySum(yc,:) = [dailySUM(ii,1) sum(monthlySUM(cnt-11:cnt,3:end),1)];
        yc = yc + 1; 
    end
    cnt=cnt+1;
end
%% Compare the plots
plot(-cumsum([0;YearlySum(:,2)]),'b')
hold on
plot(-cumsum([0;sum(YearlySum(:,3:5),2)]),'r')
plot(-cumsum([0;sum(YearlySum(:,6:26),2)]),'g')
%% write the storage time series to a js file
writeTimeSeries2JS('Jscripts/StorageBasins', 'CVHMMonthly', ones(504,1), monthlySUM(:,2), monthlySUM(:,1), cumsum([0; -monthlySUM(1:end-1,3)]/10^6/1233.48184), false);
writeTimeSeries2JS('Jscripts/StorageBasins', 'TLBMonthly', ones(504,1), monthlySUM(:,2), monthlySUM(:,1), cumsum([0; -monthlySUM(1:end-1,4)]/10^6/1233.48184), true);
writeTimeSeries2JS('Jscripts/StorageBasins', 'SJVMonthly', ones(504,1), monthlySUM(:,2), monthlySUM(:,1), cumsum([0; -monthlySUM(1:end-1,5)]/10^6/1233.48184), true);
writeTimeSeries2JS('Jscripts/StorageBasins', 'SACMonthly', ones(504,1), monthlySUM(:,2), monthlySUM(:,1), cumsum([0; -monthlySUM(1:end-1,6)]/10^6/1233.48184), true);
%%
basin_id = 2;
name_file ={'TLB','SJV','SAC'};
basin_farm_id  = unique([bud(find([bud.bas_id] == basin_id),1).farm_id])
isfirst = true;
for ii = basin_farm_id
    writeTimeSeries2JS(['Jscripts/StorageFarm' name_file{basin_id+1}], ...
        [name_file{basin_id+1} 'Farm_' num2str(ii)], ones(504,1), monthlySUM(:,2), monthlySUM(:,1), ...
        cumsum([0; -monthlySUM(1:end-1,6+ii)]/10^6/1233.48184), ~isfirst);
    isfirst = false;
end

    
%%
clear temp YearlySum
cnt=1;
yc=1;
for i=7:510
    temp(cnt,:)=[monthlySUM(i,3:6) monthlySUM(i,end) sum(monthlySUM(i,7:8))]*eomday(monthlySUM(i,1), monthlySUM(i,2));
    if cnt==12
        YearlySum(yc,:)=[monthlySUM(i,1) sum(temp,1)];
        cnt=0;
        yc=yc+1;
        clear temp
    end
    cnt=cnt+1;
end
%
%%
clear STRG
for i = 1:size(monthlySUM, 1)
    STRG(i,1) = eomday(monthlySUM(i,1), monthlySUM(i,2))*monthlySUM(i,3) + eomday(monthlySUM(i,1), monthlySUM(i,2))*monthlySUM(i,5);
end
%% average storage per year
STRG_year = sum(reshape(STRG, 12, 42),1);


%% write storage as javascript file
writeTimeSeries2JS('Jscripts/Storage', 'CVHMMonthly', ones(504,1), monthlySUM(7:end,2), monthlySUM(7:end,1), cumsum([0; -STRG(1:end-1)]/10^6/1233.48184), false);
writeTimeSeries2JS('Jscripts/Storage', 'CVHMYearly', ones(length(STRG_year),1), 2*ones(length(STRG_year),1),[1962:2003]', [0;-cumsum(STRG_year')]./10^6/1233.48184, true);

%%
plot(cumsum([0; -STRG]/10^6/1233.48184),'.-')
grid on
grid minor

%% Period 1
istart1 = 199; % Oct 1977
iend1 = 455; % Feb 1999
%% Period 2
istart2 = 367; % Noe 1991
iend2 = 509; % Sep 2003
%% find inflows and outflows for each of the selected period for the streams

days_per_month1 = eomday(STRLK.ym(istart1:iend1,1), STRLK.ym(istart1:iend1,2));
Qstream1 = bsxfun(@times, STRLK.data(:,istart1:iend1)', days_per_month1)'; %m^3/month
days_per_month2 = eomday(STRLK.ym(istart2:iend2,1), STRLK.ym(istart2:iend2,2));
Qstream2 = bsxfun(@times, STRLK.data(:,istart2:iend2)', days_per_month2)'; %m^3/month
for ii = 1:size(STRLK.rc,1)
    Q_pos1 = 0; Q_neg1 = 0; Q_pos2 = 0; Q_neg2 = 0;
    Q_pos1 = sum(Qstream1(ii,Qstream1(ii,:)>0));
    Q_neg1 = sum(Qstream1(ii,Qstream1(ii,:)<0));
    Q_pos2 = sum(Qstream2(ii,Qstream2(ii,:)>0));
    Q_neg2 = sum(Qstream2(ii,Qstream2(ii,:)<0));
    if Q_pos1 ~=0 || Q_neg1 ~= 0 || Q_pos2 ~=0 || Q_neg2 ~= 0
        ind = find(bas_rc(:,1) == STRLK.rc(ii,2) & bas_rc(:,2) == STRLK.rc(ii,3));
        bud(ind,1).QstrmPos_77_99 = Q_pos1;
        bud(ind,1).QstrmNeg_77_99 = Q_neg1;
        bud(ind,1).QstrmPos_91_03 = Q_pos2;
        bud(ind,1).QstrmNeg_91_03 = Q_neg2;
    end
end
% fill the empty values with zeros
for ii = 1:size(bud,1)
    if isempty(bud(ii,1).QstrmPos_77_99); bud(ii,1).QstrmPos_77_99 = 0;end
    if isempty(bud(ii,1).QstrmNeg_77_99); bud(ii,1).QstrmNeg_77_99 = 0;end
    if isempty(bud(ii,1).QstrmPos_91_03); bud(ii,1).QstrmPos_91_03 = 0;end
    if isempty(bud(ii,1).QstrmNeg_91_03); bud(ii,1).QstrmNeg_91_03 = 0;end
end
%% Repeat for the recharge
rch_pos1 = zeros(441,98);
rch_neg1 = zeros(441,98);
rch_pos2 = zeros(441,98);
rch_neg2 = zeros(441,98);
for ii = istart1:iend1
    days_per_month = eomday(RCH.ym(ii,1), RCH.ym(ii,2));
    id_pos = find(RCH.data{ii,1} > 0);
    id_neg = find(RCH.data{ii,1} < 0);
    rch_pos1(id_pos) = rch_pos1(id_pos) + RCH.data{ii,1}(id_pos)*days_per_month;
    rch_neg1(id_neg) = rch_neg1(id_neg) + RCH.data{ii,1}(id_neg)*days_per_month;
end
for ii = istart2:iend2
    days_per_month = eomday(RCH.ym(ii,1), RCH.ym(ii,2));
    id_pos = find(RCH.data{ii,1} > 0);
    id_neg = find(RCH.data{ii,1} < 0);
    rch_pos2(id_pos) = rch_pos2(id_pos) + RCH.data{ii,1}(id_pos)*days_per_month;
    rch_neg2(id_neg) = rch_neg2(id_neg) + RCH.data{ii,1}(id_neg)*days_per_month;
end

II = [bud.ROW]';
JJ = [bud.COLUMN_]';
for ii = 1:441
    for jj = 1:98
        id = find(II == ii & JJ == jj);
        if isempty(id)
            continue;
        end
        bud(id,1).QrchPos_77_99 = rch_pos1(ii,jj);
        bud(id,1).QrchNeg_77_99 = rch_neg1(ii,jj);
        bud(id,1).QrchPos_91_03 = rch_pos2(ii,jj);
        bud(id,1).QrchNeg_91_03 = rch_neg2(ii,jj);
    end
end
%
%% Repeat for the wells
II = [bud.ROW]';
JJ = [bud.COLUMN_]';
%bud = rmfield(bud,{'QwellPos_77_99','QwellNeg_77_99','QwellPos_91_03','QwellNeg_91_03'});
bud(1,1).QwellPos_77_99 = [];
bud(1,1).QwellNeg_77_99 = [];
bud(1,1).QwellPos_91_03 = [];
bud(1,1).QwellNeg_91_03 = [];
for ii = istart1:iend1
    ii
    days_per_month = eomday(WELLS.ym(ii,1), WELLS.ym(ii,2));
    for k = 1:2
        for jj = 1:size(WELLS.data{ii,k},1)
            id = find(II == WELLS.data{ii,k}(jj,2) & JJ == WELLS.data{ii,k}(jj,3));
            if isempty(id)
                continue;
            end
            if isempty(bud(id,1).QwellPos_77_99)
                bud(id,1).QwellPos_77_99 = 0;
            end
            if isempty(bud(id,1).QwellNeg_77_99)
                bud(id,1).QwellNeg_77_99 = 0;
            end
            if WELLS.data{ii,k}(jj,4) >= 0
                bud(id,1).QwellPos_77_99 = bud(id,1).QwellPos_77_99 + WELLS.data{ii,k}(jj,4)*days_per_month;
            else
                bud(id,1).QwellNeg_77_99 = bud(id,1).QwellNeg_77_99 + WELLS.data{ii,k}(jj,4)*days_per_month;
            end
        end
    end
    
end

for ii = istart2:iend2
    ii
    days_per_month = eomday(WELLS.ym(ii,1), WELLS.ym(ii,2));
    for k = 1:2
        for jj = 1:size(WELLS.data{ii,k},1)
            id = find(II == WELLS.data{ii,k}(jj,2) & JJ == WELLS.data{ii,k}(jj,3));
            if isempty(id)
                continue;
            end
            if isempty(bud(id,1).QwellPos_91_03)
                bud(id,1).QwellPos_91_03 = 0;
            end
            if isempty(bud(id,1).QwellNeg_91_03)
                bud(id,1).QwellNeg_91_03 = 0;
            end
            if WELLS.data{ii,k}(jj,4) >= 0
                bud(id,1).QwellPos_91_03 = bud(id,1).QwellPos_91_03 + WELLS.data{ii,k}(jj,4)*days_per_month;
            else
                bud(id,1).QwellNeg_91_03 = bud(id,1).QwellNeg_91_03 + WELLS.data{ii,k}(jj,4)*days_per_month;
            end
        end
    end
end
% fill the empty values with zeros
for ii = 1:size(bud,1)
    if isempty(bud(ii,1).QwellPos_77_99); bud(ii,1).QwellPos_77_99 = 0;end
    if isempty(bud(ii,1).QwellNeg_77_99); bud(ii,1).QwellNeg_77_99 = 0;end
    if isempty(bud(ii,1).QwellPos_91_03); bud(ii,1).QwellPos_91_03 = 0;end
    if isempty(bud(ii,1).QwellNeg_91_03); bud(ii,1).QwellNeg_91_03 = 0;end
end
%% write to shapefile
shapewrite(bud, 'gis_data/CVHM_cellBudget');
%% Caclulate average ranges over the two periods
TotDays1 = sum(eomday(STRLK.ym(istart1:iend1,1), STRLK.ym(istart1:iend1,2)));
TotDays2 = sum(eomday(STRLK.ym(istart2:iend2,1), STRLK.ym(istart2:iend2,2)));
StrmAvPos_77_99 = sum([bud.QstrmPos_77_99]'./TotDays1)/10^6;
StrmAvNeg_77_99 = sum([bud.QstrmNeg_77_99]'./TotDays1)/10^6;
StrmAvPos_91_03 = sum([bud.QstrmPos_91_03]'./TotDays2)/10^6;
StrmAvNeg_91_03 = sum([bud.QstrmNeg_91_03]'./TotDays2)/10^6;

%% Average Wells
totdays = 0;
AVwells = zeros(size(RCH.data{1,1},1), size(RCH.data{1,1},2));
for ii = istart2:iend2
    days = eomday(WELLS.ym(ii,1), WELLS.ym(ii,2));
    totdays = totdays + days;
    for jj = 1:size(WELLS.data{ii,1},1)
        AVwells(WELLS.data{ii,1}(jj,2), WELLS.data{ii,1}(jj,3)) = ...
            AVwells(WELLS.data{ii,1}(jj,2), WELLS.data{ii,1}(jj,3)) + ...
            WELLS.data{ii,1}(jj,4)*days;
    end
    for jj = 1:size(WELLS.data{ii,2},1)
        AVwells(WELLS.data{ii,2}(jj,2), WELLS.data{ii,2}(jj,3)) = ...
            AVwells(WELLS.data{ii,2}(jj,2), WELLS.data{ii,2}(jj,3)) + ... 
            WELLS.data{ii,2}(jj,4)*days;
    end
end
AVwells = AVwells/totdays;
fprintf('%f\n',sum(sum(AVwells)))
%% Average Streams
AVstrm = zeros(size(STRLK.data,1),1); %m^3/day
ym = STRLK.ym(istart2:iend2,:);
totdays = 0;
for ii = istart2:iend2
    totdays = totdays + eomday(STRLK.ym(ii,1), STRLK.ym(ii,2));
    AVstrm = AVstrm + STRLK.data(:,ii)*eomday(STRLK.ym(ii,1), STRLK.ym(ii,2)); %m^3/month
end
AVstrm = AVstrm/totdays; %m^3/day Average over the period of interest
STRMS = [STRLK.rc AVstrm];
fprintf('%f\n',sum(STRMS(:,4)))
%%
save('AvStresses','STRMS')
%% Average Recharge
AVrch = zeros(size(RCH.data{1,1},1), size(RCH.data{1,1},2));
%ym = RCH.ym(istart1:iend1,:);
totdays = 0;
for ii = istart2:iend2
    totdays = totdays + eomday(RCH.ym(ii,1), RCH.ym(ii,2));
    AVrch = AVrch + RCH.data{ii,1}*eomday(RCH.ym(ii,1), RCH.ym(ii,2)); %m^3/month;
end
AVrch = AVrch/totdays;
fprintf('%f\n',sum(sum(AVrch)))
%%
save('AvStresses','AVrch','-append')
%% Prepare recharge for writing to Scattered format
% From the gridded format convert to x,y,value format based on the actuall 
% coordinate system.
% Run the first 3 sections of the PrepareGeometryData.m to import and
% process the bas variable
% Loop through the modflow grid cells and find the center coordinate of the
% cell in the real system
load('gis_data/BAS_active_shp.mat');
bas = bas_active;
R = [bas.ROW]';
C = [bas.COLUMN_]';
xy_rch = [];
for ii = 1:size(AVrch,1)
    for jj = 1:size(AVrch,2)
        id = find(R == ii & C == jj);
        if isempty(id)
            continue;
        elseif length(id) == 1
            xc = mean(bas(id,1).X(1:end-1));
            yc = mean(bas(id,1).Y(1:end-1));
            xy_rch = [xy_rch; xc yc AVrch(ii,jj)];
        else
            error('Thats so wrong');
        end
    end
end
%% convert volume to rate
xy_rch(:,3) = xy_rch(:,3)./(1609.344*1609.344);
%% 
% add the buffer nodes.
% run the section from PrepareGeometryData.m which creates the buff_pnt
% variable
% create an interpolant with the existing data
Frch = scatteredInterpolant(xy_rch(:,1), xy_rch(:,2), xy_rch(:,3));
Frch.Method = 'nearest';
Frch.ExtrapolationMethod = 'nearest';
buf_val = Frch(buff_pnt(:,1), buff_pnt(:,2));
xy_rch = [xy_rch; buff_pnt buf_val];
%% write the rch file
writeScatteredData('CVHM_RCH_per2.npsat', struct('PDIM',2,'TYPE','HOR','MODE','SIMPLE'), xy_rch);
%}
%% Average Well data
AVwell = zeros(size(RCH.data{1,1},1), size(RCH.data{1,1},2));
totdays = 0;
for ii = istart2:iend2
    ii
    monthDays = eomday(WELLS.ym(ii,1), WELLS.ym(ii,2));
    totdays = totdays + monthDays;
    rc = WELLS.data{ii,1}(:,2:3);
    rc = unique(rc, 'rows');
    for jj = 1:size(rc,1)
        id = find(WELLS.data{ii,1}(:,2) == rc(jj,1) & WELLS.data{ii,1}(:,3) == rc(jj,2));
        AVwell(rc(jj,1),rc(jj,2)) = AVwell(rc(jj,1),rc(jj,2)) + sum(WELLS.data{ii,1}(id,4))*monthDays;
    end
    
    rc = WELLS.data{ii,2}(:,2:3);
    rc = unique(rc, 'rows');
    for jj = 1:size(rc,1)
        id = find(WELLS.data{ii,2}(:,2) == rc(jj,1) & WELLS.data{ii,2}(:,3) == rc(jj,2));
        AVwell(rc(jj,1),rc(jj,2)) = AVwell(rc(jj,1),rc(jj,2)) + sum(WELLS.data{ii,2}(id,4))*monthDays;
    end
end
AVwell = AVwell/totdays;

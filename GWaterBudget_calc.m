% path where the cbc flow data in matlab format are located
% Claudia's Caped Run
%
cbcf_path = '/home/giorgk/Documents/UCDAVIS/CVHM_DATA/ClaudiaRun/';
m=4;y=1961;
t2=0;
for i=1:510
    if m==13
        y
        m=1;y=y+1;
    end
    cnst=floor((i-1)/10);
    if t2<cnst*10+1
        clear OUT
        load([cbcf_path 'cbcf_' num2str(cnst*10+1) '_' num2str(10*cnst+10) '.mat']);
        t2=cnst*10+1;
    end
    monthlySUM(i,1)=y;
    monthlySUM(i,2)=m;
    monthlySUM(i,3)=sum(sum(sum(OUT{i,1}{1,2})));%Storage
    monthlySUM(i,4)=sum(OUT{i,1}{6,2}(:,4));%Head Dep bounds
    monthlySUM(i,5)=sum(sum(sum(OUT{i,1}{7,2})));%Inst. IB Storage
    monthlySUM(i,6)=sum(OUT{i,1}{8,2}(:,4));%Stream Leackage
    monthlySUM(i,7)=sum(OUT{i,1}{9,2}(:,4));%MNW
    monthlySUM(i,8)=sum(OUT{i,1}{10,2}(:,4));%Farm wells
    monthlySUM(i,9)=sum(OUT{i,1}{11,2}(:,2));%Recharge
    if i == 1
        STRLK.rc = OUT{i,1}{8,2}(:,1:3);
        STRLK.data = zeros(size(OUT{i,1}{8,2}(:,1:3),1),510);
        STRLK.ym = zeros(510,2);
    end
    STRLK.data(:,i) = OUT{i,1}{8,2}(:,4);
    STRLK.ym(i,:) = [y m];
    RCH.data{i,1} = reshape(OUT{i,1}{11,2}(:,2),98,441)';
    RCH.ym(i,:) = [y m];
    WELLS.data{i,1} = OUT{i,1}{9,2};%MNW
    WELLS.data{i,2} = OUT{i,1}{10,2};%FARM wells
    WELLS.ym(i,:) = [y m];
    m=m+1;
end
%
% From file https://water.usgs.gov/ogw/modflow/mf6io.pdf pg21:
% The cell-by-cell storage term gives the net flow to or from storage in a variable-head cell. The net storage
% for each cell in the grid is saved in transient simulations if the appropriate flags are set. Withdrawal from storage
% in the cell is considered positive, whereas accumulation in storage is considered negative.
%% CVHM default run
load('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/CBC.mat');
%%
m=4;y=1961;
for ii=1:510
    if m==13
        m=1;y=y+1;
    end
    monthlySUM(ii,1)=y;
    monthlySUM(ii,2)=m;
    monthlySUM(ii,3)=sum(sum(sum(CBC{ii,1}{1,2}))) + sum(sum(sum(CBC{ii,1}{7,2})));
    m=m+1;
end
%%
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
%% load the cell shapefile
% and write the averaged info
bas = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/BAS_active.shp');
bas_rc = [[bas.ROW]' [bas.COLUMN_]'];
bud = rmfield(bas, {'lay1_act','lay2_act','lay3_act','lay45_act','lay6_act','lay7_act','lay8_act','lay9_act','lay10_act',...
        'strt_hd_01','strt_hd_02','strt_hd_03','strt_hd_04','strt_hd_05','strt_hd_06','strt_hd_07','strt_hd_08','strt_hd_09','strt_hd_10'});
%% add Subbasins and farm id
farms = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/FARMS_poly.shp');
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
%}
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


%% Average Streams
AVstrm = STRLK.data(:,istart2:iend2); %m^3/day
ym = STRLK.ym(istart2:iend2,:);
totdays = 0;
for ii = 1:size(AVstrm,2)
    totdays = totdays + eomday(ym(ii,1), ym(ii,2));
    AVstrm(:,ii) = AVstrm(:,ii)*eomday(ym(ii,1), ym(ii,2)); %m^3/month
end
AVstrm = sum(AVstrm,2)/totdays; %m^3/day Average over the period of interest
STRMS = [STRLK.rc AVstrm];
save('AvStresses','STRMS')
%% Average Recharge
AVrch = zeros(size(RCH.data{1,1},1), size(RCH.data{1,1},2));
ym = RCH.ym(istart2:iend2,:);
totdays = 0;
for ii = 1:size(ym,1)
    totdays = totdays + eomday(ym(ii,1), ym(ii,2));
    AVrch = AVrch + RCH.data{ii,1}*eomday(ym(ii,1), ym(ii,2)); %m^3/month;
end
AVrch = AVrch/totdays;
%%
save('AvStresses','AVrch','-append')
%% Prepare recharge for writing to Scattered format
% From the gridded format convert to x,y,value format based on the actuall 
% coordinate system.
% Run the first 3 sections of the PrepareGeometryData.m to import and
% process the bas variable
% Loop through the modflow grid cells and find the center coordinate of the
% cell in the real system
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



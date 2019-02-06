% path where the cbc flow data in matlab format are located
% Claudia's Caped Run
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
for i = 7:size(monthlySUM, 1)
    STRG(i-6,1) = eomday(monthlySUM(i,1), monthlySUM(i,2))*monthlySUM(i,3) + eomday(monthlySUM(i,1), monthlySUM(i,2))*monthlySUM(i,5);
end
%% average storage per year
STRG_year = sum(reshape(STRG, 12, 42),1);


%% write storage as javascript file
writeTimeSeries2JS('Jscripts/MonthlyStorage', 'StorageMonthly', ones(504,1), monthlySUM(7:end,2), monthlySUM(7:end,1), cumsum([0; -STRG(1:end-1)]/10^6/1233.48184), false);
writeTimeSeries2JS('Jscripts/MonthlyStorage', 'StorageYearly', ones(length(STRG_year),1), 2*ones(length(STRG_year),1),[1962:2003]', [0;-cumsum(STRG_year')]./10^6/1233.48184, true);

%%
plot(cumsum([0; -STRG]/10^6/1233.48184),'.-')
grid on
grid minor
%% Period 1
istart = 200; % Noe 1977
iend = 456; % Mar 1999
%% Period 2
istart = 368; % Noe 1991
iend = 503; % Sep 2003
%% Average Streams
AVstrm = STRLK.data(:,istart:iend); %m^3/day
ym = STRLK.ym(istart:iend,:);
totdays = 0;
for ii = 1:size(AVstrm,2)
    totdays = totdays + eomday(ym(ii,1), ym(ii,2));
    AVstrm(:,ii) = AVstrm(:,ii)*eomday(ym(ii,1), ym(ii,2)); %m^3/month
end
AVstrm = sum(AVstrm,2)/totdays; %m^3/day Average over the period of interest
STRMS = [STRLK.rc AVstrm];
save('AvStresses','STRMS','-append')
%% Average Recharge
AVrch = zeros(size(RCH.data{1,1},1), size(RCH.data{1,1},2));
ym = RCH.ym(istart:iend,:);
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
writeScatteredData('CVHM_RCH.npsat', struct('PDIM',2,'TYPE','HOR','MODE','SIMPLE'), xy_rch);



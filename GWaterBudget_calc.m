% path where the cbc flow data in matlab format are located
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
    m=m+1;
end
%
% From file https://water.usgs.gov/ogw/modflow/mf6io.pdf pg21:
% The cell-by-cell storage term gives the net flow to or from storage in a variable-head cell. The net storage
% for each cell in the grid is saved in transient simulations if the appropriate flags are set. Withdrawal from storage
% in the cell is considered positive, whereas accumulation in storage is considered negative.
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
for i = 1:size(monthlySUM, 1)
    STRG(i,1) = eomday(monthlySUM(i,1), monthlySUM(i,2))*monthlySUM(i,3);
end
%%
plot(cumsum([0; -STRG]/10^6/1233.48184),'.-')
grid on
grid minor
%% 
istart = 199; % Oct 1977
iend = 503; % Feb 2003
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

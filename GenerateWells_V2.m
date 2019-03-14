%% WELL GENERATION ALGORITHM V2 
%
% In this version we are going to take into account the pumping 
% distribution of CVHM and combined information from the well dataset and
% CVHM model
%% Load Well data set
% bot = depth from land surface to the bottom of the perforated interval [ft]
% top = depth from land surface to the top of the perforated interval [ft]
% bot - top is the length of the screened interval
% Q pumping in gpm
% year of construction
% S = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CV_only_PubAg.shp');
%S = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CVpbAG_mth.shp');
load('/media/giorgk/DATA/giorgk/Documents/CVHM_DATA/WellData/CVpbAG_mth_shp.mat')
%% Process well data set
WellTable = readtable('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/giorgos_with_method_ll.csv');
%% Create shapefile from the table
for ii = 1:size(WellTable,1)
    S(ii,1).Geometry = 'Point';
    S(ii,1).X = WellTable.lon(ii,1);
    S(ii,1).Y = WellTable.lat(ii,1);
    S(ii,1).WCRNumber = cell2mat(WellTable.WCRNumber(ii,1));
    S(ii,1).top = str2double(cell2mat(WellTable.top(ii,1)));
    S(ii,1).bot = str2double(cell2mat(WellTable.bot(ii,1)));
    S(ii,1).depth = str2double(cell2mat(WellTable.depth(ii,1)));
    S(ii,1).year = str2double(cell2mat(WellTable.year(ii,1)));
    S(ii,1).Q = str2double(cell2mat(WellTable.Q(ii,1)));
    S(ii,1).type = str2double(cell2mat(WellTable.type(ii,1)));
    if strcmp(cell2mat(WellTable.type(ii,1)), 'agriculture')
        S(ii,1).type = 0;
    elseif strcmp(cell2mat(WellTable.type(ii,1)), 'public')
        S(ii,1).type = 1;
    else
        S(ii,1).type = -9;
    end
    
    if strcmp(cell2mat(WellTable.MethodofDeterminationLL(ii,1)), 'Derived from Address') || ...
           strcmp(cell2mat(WellTable.MethodofDeterminationLL(ii,1)), 'GPS')
        S(ii,1).method = 0;
    else
        S(ii,1).method = 1;
    end
    
end
shapewrite(S, '/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CVpbAG_mth');
%% Well Age analysis
% discard the wells without year of construction.
% discard prehistoric and other old wells prior to 1900
% discard future wells
id_dlt = isnan([S.year]');
S(id_dlt,:) = [];
id_dlt = [S.year]' < 1900;
S(id_dlt,:) = [];
id_dlt = [S.year]' > 2018;
S(id_dlt,:) = [];
%% plot histogram of age
H = histogram([S.year]','BinWidth',5);
xlabel('Year of construction');
% print as js data set
for ii = 1:length(H.BinEdges)-1
    fprintf('{ x: new Date( %d , 0), y: %d },\n', [H.BinEdges(ii) H.Values(ii)]);
end
fprintf('{ x: new Date( %d , 0), y: %d }\n', [H.BinEdges(end) H.Values(end)]);
%% Find the number of wells for a given threshold
age_thres = 25;
for age_thres = 1:81
    Nwells(age_thres,1) = length(find([S.year]' > 2018 - age_thres));
end
%%
yrNw = [];
for yr = 1900:2018
    yrNw = [yrNw; yr length(find([S.year]' > yr))];
    fprintf('{ x: new Date( %d , 0), y: %.3f },\n', [yrNw(end,1) yrNw(end,2)/1000]);
end
%%
H = histogram([S.year]','BinWidth',1,'Normalization','cumcount');
plot(1:81, Nwells/1000,'linewidth',2)
xlabel('Well Age [year]')
ylabel('# of wells [10^3]')
grid on
grid minor
%% spacial distribution for a given age threshold
age_thres = 30;
id_dlt = [S.year]' < 2018 - age_thres;
Sact = S;
Sact(id_dlt,:) = [];
plot([Sact.X]',[Sact.Y]','.')
axis equal
axis off
%% save to local shapefile.
shapewrite(Sact, '/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CVwells_30y')
%% CVHM WELL PUMPING
% starting with the Bas active shapefile we will add all the well info on
% each cell.
Bas = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/BAS_active.shp');
%%
load('gis_data/BAS_active_shp')
Bas = bas_active;
for ii = [1 2 3 45 6 7 8 9 10]
    Bas = rmfield(Bas,['lay' num2str(ii) '_act']);
end
for ii = 1:10
    Bas = rmfield(Bas,['strt_hd_' num2str(ii,'%02.0f')]);
end
Bas = rmfield(Bas, 'grid_active');

CVHM_wells = Bas;
%% load Well Budgets
cbcf_path = '/home/giorgk/Documents/UCDAVIS/CVHM_DATA/ClaudiaRun/';
m=4;y=1961;
t2=0;
for ii=1:510
    if m==13
        y
        m = 1;
        y = y+1;
    end
    cnst=floor((ii-1)/10);
    if t2<cnst*10+1
        clear OUT
        load([cbcf_path 'cbcf_' num2str(cnst*10+1) '_' num2str(10*cnst+10) '.mat']);
        t2=cnst*10+1;
    end
    Welldates(ii,1)=y;
    Welldates(ii,2)=m;
    FARMWELLS{ii,1} = OUT{ii,1}{10,2};
    MNWELLS{ii,1} = OUT{ii,1}{9,2};
    m = m + 1;
end
%% OR
%load('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/Well_Budget.mat');
load('/media/giorgk/DATA/giorgk/Documents/CVHM_DATA/WellData/Well_Budget.mat');
%% 
Ndays = 0;
for ii = 1:size(Welldates,1)
    Ndays(ii,1) = eomday(Welldates(ii,1),Welldates(ii,2));
end
%
%% Average the stresses for the simulation period
% get istart and iend from GWaterBudget_calc
MNW = zeros(441, 98);
FRW = zeros(441, 98);
totdays = 0;
for ii = istart2:iend2
    totdays = totdays + Ndays(ii,1);
    for jj = 1:size(FARMWELLS{ii,1},1)
        FRW(FARMWELLS{ii,1}(jj,2), FARMWELLS{ii,1}(jj,3)) = ...
            FRW(FARMWELLS{ii,1}(jj,2), FARMWELLS{ii,1}(jj,3)) + ...
            FARMWELLS{ii,1}(jj,4)*Ndays(ii,1);
    end
    for jj = 1:size(MNWELLS{ii,1},1)
        MNW(MNWELLS{ii,1}(jj,2), MNWELLS{ii,1}(jj,3)) = ...
            MNW(MNWELLS{ii,1}(jj,2), MNWELLS{ii,1}(jj,3)) + ...
            MNWELLS{ii,1}(jj,4)*Ndays(ii,1);
    end
end
FRW = FRW./totdays;
MNW = MNW./totdays;
%% add the well rates to CVHM_wells shapefiles
for ii = 1:size(CVHM_wells,1)
    CVHM_wells(ii,1).farm = FRW(CVHM_wells(ii,1).ROW, CVHM_wells(ii,1).COLUMN_);
    CVHM_wells(ii,1).mnw = MNW(CVHM_wells(ii,1).ROW, CVHM_wells(ii,1).COLUMN_);
end
%%
shapewrite(CVHM_wells, 'gis_data/CVHM_wellBud');
%
%% Find the volume of water that is used for Urban and Ag for each farm
%CVHM_wells = shaperead('gis_data/CVHM_wellBud');
load('gis_data/FMP_shp'); %Farms = shaperead('gis_data/FMP.shp');
load('gis_data/FARMS_poly_shp'); %Farms_poly = shaperead('gis_data/FARMS_poly.shp');
for ii = 1:length(Farms_poly)
    Farms_poly(ii,1).UrbanPump = 0;
    Farms_poly(ii,1).AgPump = 0;
end
ROWfrm = [Farms.ROW]';
COLfrm = [Farms.COLUMN_]';
%%
for ii = 1:size(CVHM_wells,1)
    r = CVHM_wells(ii,1).ROW;
    c = CVHM_wells(ii,1).COLUMN_;
    ir = find(ROWfrm == r & COLfrm == c);
    ifarm = find([Farms_poly.dwr_sbrgns]' == Farms(ir,1).dwr_sbrgns);
    if Farms(ir,1).LU_2000 == 2
        Farms_poly(ifarm,1).UrbanPump = Farms_poly(ifarm,1).UrbanPump + ...
            CVHM_wells(ii,1).farm + CVHM_wells(ii,1).mnw;
    else
        Farms_poly(ifarm,1).AgPump = Farms_poly(ifarm,1).AgPump + ...
            CVHM_wells(ii,1).farm + CVHM_wells(ii,1).mnw;
    end
end
%%
shapewrite(Farms_poly, 'gis_data/FARMS_polyPump.shp');
%
%% Generate Pumping wells
Farm_basin = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/FARMS_poly.shp');
Farms_poly = shaperead('gis_data/FARMS_polyPump.shp');
CVHMwells = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/WellData/CVwells_30y_proj.shp');
cvhmDomain = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/BAS_active_outline.shp');
%
%%
%
clear
clc
do_plot = false;
load('gis_data/WellGenerationSHP','Farm_basin','cvhmDomain','CVHMwells');
Farm_basin = Farm_basin';
CVHMwells = CVHMwells';
load('gis_data/FARMS_polyPump_91_03.mat')
load('gis_data/CVHM_Mesh_outline_modif_shp.mat')

% RCH : 25715545.162923 m^3/day
% STRM: 899912.211251 m^3/day
Total_pump_vol = 25715545.162923 + 899912.211251; %m^3/day
Act_pump_vol = abs(sum([Farms_poly.AgPump]')+sum([Farms_poly.UrbanPump]'));
corr_ratio = Total_pump_vol/Act_pump_vol;

% prepare the streams for vectorized distance calculations
load('gis_data/CVHM_streams_shp');
stream_L = [];
for ii = 1:length(CVHM_streams)
    for jj = 1:length(CVHM_streams(ii,1).X)-1
        l = [CVHM_streams(ii,1).X(jj) CVHM_streams(ii,1).Y(jj) CVHM_streams(ii,1).X(jj+1) CVHM_streams(ii,1).Y(jj+1)];
        if ~any(isnan(l))
            stream_L = [stream_L;l];
        end
    end
end

%
% Generate random pumping
Ql = 0.15;
Qh = 0.925;
WX = nan(20000,1); WY = nan(20000,1); cnt = 1;
clear Wells_Gen;
for ii = 1:length(Farms_poly)
    if Farms_poly(ii,1).dwr_sbrgns == 0
        continue;
    end
    disp(['Farm: ' num2str(Farms_poly(ii,1).dwr_sbrgns) ' [' num2str(ii) ']'])
    in_wells = inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
                         Farms_poly(ii,1).X, Farms_poly(ii,1).Y) & ... 
               inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
                         cvhmDomain.X, cvhmDomain.Y);
    tempW = CVHMwells(in_wells,1);
    ipb = find([tempW.type]'==1);
    iag = find([tempW.type]'==0);

    tempWpb = tempW(ipb,1);
    tempWag = tempW(iag,1);
    
    figure(1); clf
    figure(1); plot([tempWpb.X]',[tempWpb.Y]','or')
    figure(1); hold on
    figure(1); plot([tempWag.X]',[tempWag.Y]','.b')
    figure(1); axis equal
    AG_pdf = Calc_2D_PDF([tempWag.X]',[tempWag.Y]', 15);
    PB_pdf = Calc_2D_PDF([tempWpb.X]',[tempWpb.Y]', 15);
    
    figure(2); clf
    figure(2); plot([tempWag.X]',[tempWag.Y]','+r')
    figure(2);hold on
    figure(2);surf(AG_pdf.X,AG_pdf.Y,AG_pdf.V,'edgecolor','none')
    figure(2); view(0,90)
    figure(2); alpha(0.5)
    figure(2); axis equal
    drawnow
    
    figure(3); clf
    figure(3); plot([tempWpb.X]',[tempWpb.Y]','+r')
    figure(3);hold on
    figure(3);surf(PB_pdf.X,PB_pdf.Y,PB_pdf.V,'edgecolor','none')
    figure(3); view(0,90)
    figure(3); alpha(0.5)
    figure(3); axis equal
    drawnow
    
    
    % Pumping distribution (different for Ag-Urban wells unless there are
    % very few urban wells in the farm.
    Qag = [tempW(iag,1).Q]';
    Qpb = [tempW(ipb,1).Q]';
    Qag(isnan(Qag), :) = [];
    Qpb(isnan(Qpb), :) = [];
    Qag(Qag <= 0, :) = [];
    Qpb(Qpb <= 0, :) = [];
    [xag, fag] = ecdf(log10(Qag));
    if length(Qpb) < 5
        Qpb = [Qpb;Qag];
    end
    [xpb, fpb] = ecdf(log10(Qpb));
    figure(4); clf
    figure(4); plot(fag, xag, 'b')
    figure(4); hold on
    figure(4); plot(fpb, xpb, 'r')
    drawnow
    
    %Pumping - depth distribution
    Qtemp = [tempW.Q]';
    Dtemp = [tempW.depth]';
    id = find(~isnan(Qtemp) & ~isnan(Dtemp));
    Qtemp = Qtemp(id);
    Dtemp = Dtemp(id);
    id = find(Qtemp > 0 & Dtemp > 0);
    Qtemp = Qtemp(id);
    Dtemp = Dtemp(id);
    QD_pdf = Calc_2D_PDF(log10(Qtemp), log10(Dtemp), 15);
    figure(5); clf
    figure(5); plot(log10(Qtemp), log10(Dtemp),'+r')
    figure(5);hold on
    figure(5);surf(QD_pdf.X, QD_pdf.Y, QD_pdf.V, 'edgecolor','none')
    figure(5); view(0,90)
    figure(5); alpha(0.5)
    figure(5); axis equal
    drawnow
    
    % Depth - screen length
    Dtemp = [tempW.depth]';
    SLtemp = [tempW.bot]' - [tempW.top]';
    id = find(~isnan(Dtemp) & ~isnan(SLtemp));
    Dtemp = Dtemp(id);
    SLtemp = SLtemp(id);
    id = find(Dtemp > 0 & SLtemp > 0);
    Dtemp = Dtemp(id);
    SLtemp = SLtemp(id);
    DS_pdf = Calc_2D_PDF(log10(Dtemp), log10(SLtemp), 15);
    figure(6); clf
    figure(6); plot(log10(Dtemp), log10(SLtemp),'+r')
    figure(6);hold on
    figure(6);surf(DS_pdf.X, DS_pdf.Y, DS_pdf.V, 'edgecolor','none')
    figure(6); view(0,90)
    figure(6); alpha(0.5)
    figure(6); axis equal
    drawnow
    

    
    Qag_farm = 0;
    Qpb_farm = 0;
    
    % The data in Rich's well data set are in gpm while the units of CVHM are
    % in m^3/day. To convert the pumping to m^/day Q = Q*5.451.
    % However this correpond to the design pumping. The argument is that
    % each year only six months operate with this rate therefor we divide
    % the pumping by 6 months Q = Q*5.451/6. 
    % As we dont want rates lower than 400 m^3/day this would be
    % log10(400*(6/5.451)) = 2.6437 gpm.
    % 10^2.6437*5.451/6
    
    if isempty(find(fag < 2.6437, 1, 'last'))
        Ql = xag(1);
    else
        Ql = xag(find(fag < 2.6437, 1, 'last')+1);
    end
    while Qag_farm < abs(Farms_poly(ii,1).AgPump)*corr_ratio
        % generate a well location
        while 1 
            xr = Farms_poly(ii,1).BoundingBox(1,1) + (Farms_poly(ii,1).BoundingBox(2,1) - Farms_poly(ii,1).BoundingBox(1,1))*rand;
            yr = Farms_poly(ii,1).BoundingBox(1,2) + (Farms_poly(ii,1).BoundingBox(2,2) - Farms_poly(ii,1).BoundingBox(1,2))*rand;
            in = inpolygon(xr, yr, Farms_poly(ii,1).X, Farms_poly(ii,1).Y);
            if ~in; continue;end
            in = inpolygon(xr, yr, CVHM_Mesh_outline_modif.X, CVHM_Mesh_outline_modif.Y);
            if ~in; continue;end
            if nanmin(sqrt((WX - xr).^2 + (WY - yr).^2)) < 400; continue;end
            if  min(Dist_Point_LineSegment(xr, yr, stream_L)) < 400; continue;end
            
            r_accept = rand;
            if r_accept < AG_pdf.F(xr, yr)
                break;
            end
        end
        
        

            
        % Generate pumping Depth screen length
        rq = Ql + (Qh - Ql)*rand;
        Qr = interp1(xag, fag, rq);
        if abs(Farms_poly(ii,1).AgPump)*corr_ratio - (Qag_farm + 10^Qr*5.451/6) < 400
            Qr = log10((abs(Farms_poly(ii,1).AgPump)*corr_ratio - Qag_farm)*6/5.451);
        end
        
        [ Qw, Dw, Sw ] = AssignQDS_v2(QD_pdf, DS_pdf, Qr, nan, nan);
        Qag_farm = Qag_farm + 10^Qw*5.451/6; % m^3/day for six months
        
        
        WX(cnt,1) = xr; WY(cnt,1) = yr;
        
        if do_plot
            figure(2);plot(xr,yr,'.b')
            figure(5);plot(Qw,Dw,'.b')
            figure(6);plot(Dw,Sw,'.b')
        end
        figure(4);title([num2str(abs(Farms_poly(ii,1).AgPump)*corr_ratio - Qag_farm) ' ' num2str(cnt)]);
        Wells_Gen(cnt,1).X = xr;
        Wells_Gen(cnt,1).Y = yr;
        Wells_Gen(cnt,1).Q = 10^Qw*5.451/6;
        Wells_Gen(cnt,1).D = 10^Dw;
        Wells_Gen(cnt,1).SL = 10^Sw;
        Wells_Gen(cnt,1).type = 0;
        cnt = cnt + 1;
        drawnow
    end
    
    if isempty(find(fpb < 2.6437, 1, 'last'))
        Ql = xpb(1);
    else
        Ql = xpb(find(fpb < 2.6437, 1, 'last')+1);
    end
    while Qpb_farm < abs(Farms_poly(ii,1).UrbanPump)*corr_ratio
        while 1 
            xr = Farms_poly(ii,1).BoundingBox(1,1) + (Farms_poly(ii,1).BoundingBox(2,1) - Farms_poly(ii,1).BoundingBox(1,1))*rand;
            yr = Farms_poly(ii,1).BoundingBox(1,2) + (Farms_poly(ii,1).BoundingBox(2,2) - Farms_poly(ii,1).BoundingBox(1,2))*rand;
            in = inpolygon(xr, yr, Farms_poly(ii,1).X, Farms_poly(ii,1).Y);
            if ~in; continue;end
            in = inpolygon(xr, yr, CVHM_Mesh_outline_modif.X, CVHM_Mesh_outline_modif.Y);
            if ~in; continue;end
            if nanmin(sqrt((WX - xr).^2 + (WY - yr).^2)) < 400; continue;end
            if  min(Dist_Point_LineSegment(xr, yr, stream_L)) < 400; continue;end
            
            r_accept = rand;
            if r_accept < PB_pdf.F(xr, yr)
                break;
            end
        end

            
        % Generate pumping Depth screen length
        rq = Ql + (Qh - Ql)*rand;
        Qr = interp1(xpb, fpb, rq);
        if abs(Farms_poly(ii,1).UrbanPump)*corr_ratio - Qpb_farm - 10^Qr*5.451/6 < 400
            Qr = log10((abs(Farms_poly(ii,1).UrbanPump)*corr_ratio - Qpb_farm)*6/5.451);
        end
        
        [ Qw, Dw, Sw ] = AssignQDS_v2(QD_pdf, DS_pdf, Qr, nan, nan);
        Qpb_farm = Qpb_farm + 10^Qw*5.451/6;
        WX(cnt,1) = xr; WY(cnt,1) = yr;
        
        if do_plot
            figure(3);plot(xr,yr,'.b')
            figure(5);plot(Qw,Dw,'.b')
            figure(6);plot(Dw,Sw,'.b')
        end
        figure(4);title([num2str(abs(Farms_poly(ii,1).UrbanPump)*corr_ratio - Qpb_farm) ' ' num2str(cnt)] );
        Wells_Gen(cnt,1).X = xr;
        Wells_Gen(cnt,1).Y = yr;
        Wells_Gen(cnt,1).Q = 10^Qw*5.451/6;
        Wells_Gen(cnt,1).D = 10^Dw;
        Wells_Gen(cnt,1).SL = 10^Sw;
        Wells_Gen(cnt,1).type = 1;
        cnt = cnt + 1;
        drawnow 
    end
end
%
%
%% 
do_plot = true;
clear WELLS4CVHM
for ii = 1:length(Farms_poly)
    if Farms_poly(ii,1).dwr_sbrgns == 0
        continue;
    end
        
    disp(['Farm: ' num2str(Farms_poly(ii,1).dwr_sbrgns) ' [' num2str(ii) ']'])
    in_wells = inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
                         Farms_poly(ii,1).X, Farms_poly(ii,1).Y) & ... 
               inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
                         cvhmDomain.X, cvhmDomain.Y);
                     
    tempW = CVHMwells(in_wells,1);
    
    negQ = find([tempW.Q]' <= 0);
    for jj = 1:length(negQ)
        tempW(negQ(jj),1).Q = nan;
    end
    negD = find([tempW.depth]' <= 0);
    for jj = 1:length(negD)
        tempW(negD(jj),1).depth = nan;
    end
    negSL = find([tempW.bot]' - [tempW.top]' <= 0);
    for jj = 1:length(negSL)
        tempW(negSL(jj),1).top = nan;
        tempW(negSL(jj),1).bot = nan;
    end
    
    
    
    
    QD_PDF = Calc_2D_PDF(log10([tempW.Q]'), log10([tempW.depth]'));
    DS_PDF = Calc_2D_PDF(log10([tempW.depth]'), log10([tempW.bot]' - [tempW.top]'));
    if do_plot
        figure(1); clf
        figure(1); hold on
        figure(1);surf(QD_PDF.X,QD_PDF.Y,QD_PDF.V,'edgecolor','none')
        figure(1);plot(log10([tempW.Q]'), log10([tempW.depth]'),'.')
        figure(1); view(0,90)
        figure(1); alpha(0.5)
        drawnow;
        
        figure(2); clf
        figure(2); hold on
        figure(2);surf(DS_PDF.X,DS_PDF.Y,DS_PDF.V,'edgecolor','none')
        figure(2);plot(log10([tempW.depth]'), log10([tempW.bot]' - [tempW.top]'),'.')
        figure(2); view(0,90)
        figure(2); alpha(0.5)
        drawnow;
    end
    
    tempW_mod = tempW;
    for jj = 1:length(tempW)
        jj
        [ Qr, Dr, Sr ] = AssignQDS_v2(QD_PDF, DS_PDF, ...
                                    log10(tempW(jj,1).Q), ...
                                    log10(tempW(jj,1).depth), ...
                                    log10(tempW(jj,1).bot - tempW(jj,1).top));
         tempW_mod(jj,1).Q = 10^Qr;
         tempW_mod(jj,1).depth = 10^Dr;
         tempW_mod(jj,1).SL = 10^Sr;
    end
    
    if do_plot
        figure(1);plot(log10([tempW_mod.Q]'), log10([tempW_mod.depth]'),'.r')
        figure(2);plot(log10([tempW_mod.depth]'), log10([tempW_mod.SL]'),'.r')
        drawnow
    end
        
    
    
     
%     ipub = find([tempW.type]'==1);
%     iag = find([tempW.type]'==0);
%     tempWpb = tempW(ipub,1);
%     tempWag = tempW(iag,1);
%     figure(1);clf
%     figure(1); plot([tempWag.Q]',[tempWag.depth]','.b')
%     figure(1); hold on
%     figure(1); plot([tempWpb.Q]',[tempWpb.depth]','.r')
%     
%     figure(2);clf
%     figure(2); plot([tempWag.depth]',[tempWag.bot]' - [tempWag.top]','.b')
%     figure(2); hold on
%     figure(2); plot([tempWpb.depth]',[tempWpb.bot]' - [tempWpb.top]','.r')
%     continue;
    
    
%     tempWag_mod = tempWag;
%     tempWpb_mod = tempWpb;
%     if length(tempWpb_mod) < 15
%         % use all public records that belong to basin
%         ibas = Farm_basin(ii,1).Basins;
%         in_basin = false(length(CVHMwells),1);
%         for kk = 1:length(Farm_basin)
%             if Farm_basin(kk,1).Basins ~= ibas
%                 continue;
%             end
%             in_basin(inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
%                          Farm_basin(kk,1).X, Farm_basin(kk,1).Y)) = true;
%         end
%         tempWbas = CVHMwells(in_basin,1);
%         ipub = find([tempWbas.type]'==1);
%         tempWpb = tempWbas(ipub,1);
%     end
% 
%     % ----------AG wells --------
%     %%% Pumping ECDF
%     idQ = ~isnan([tempWag.Q]') & [tempWag.Q]' > 0; % isolate the records that have pumping
%     display(sum(idQ))
%     logQ = log10([tempWag(idQ,1).Q]');
%     [fQ, xQ]=ecdf(logQ); % create an ECDF based on the data
%     
%     %%% Depth ECDF
%     idD = ~isnan([tempWag.depth]') & [tempWag.depth]' > 0;
%     display(sum(idD))
%     logD = log10([tempWag(idD,1).depth]');
%     [fD, xD]=ecdf(logD);
%     
%     %%% Screen length ECDF
%     tempSL = [tempWag.bot]' - [tempWag.top]';
%     idS = ~isnan(tempSL) & tempSL < 2000  & tempSL > 0;
%     display(sum(idS))
%     logS = log10([tempWag(idS,1).bot]' - [tempWag(idS,1).top]');
%     [fS, xS]=ecdf(logS);
%     
%     clear Stats
%     Stats.fQ = fQ; Stats.lowQ = 0.15;
%     Stats.xQ = xQ; Stats.uppQ = 0.85;
%     Stats.fD = fD; Stats.lowD = 0;
%     Stats.xD = xD; Stats.uppD = 1;
%     Stats.fS = fS; Stats.lowS = 0;
%     Stats.xS = xS; Stats.uppS = 1;
%     
%     for jj = 1:length(tempWag_mod)
%         [ tempWag_mod(jj,1).Q, ...
%           tempWag_mod(jj,1).depth, ...
%           tempWag_mod(jj,1).SL ] = AssignQDS(Stats, ...
%                                              tempWag_mod(jj,1).Q, ...
%                                              tempWag_mod(jj,1).depth, ...
%                                              tempWag_mod(jj,1).bot - tempWag_mod(jj,1).top);
%     end
% 
%     
%     % --------------- Urban wells-------------
%     %%% Pumping ECDF
%     idQ = ~isnan([tempWpb.Q]') & [tempWpb.Q]' > 0; % isolate the records that have pumping
%     display(sum(idQ))
%     if sum(idQ) < 15
%         % use all public records that belong to basin
%         ibas = Farm_basin(ii,1).Basins;
%         in_basin = false(length(CVHMwells),1);
%         for kk = 1:length(Farm_basin)
%             if Farm_basin(kk,1).Basins ~= ibas
%                 continue;
%             end
%             in_basin(inpolygon([CVHMwells.X]',[CVHMwells.Y]', ...
%                          Farm_basin(kk,1).X, Farm_basin(kk,1).Y)) = true;
%         end
%         tempWbas = CVHMwells(in_basin,1);
%         ipub = find([tempWbas.type]'==1);
%         tempWpbQ = tempWbas(ipub,1);
%         idQ = ~isnan([tempWpbQ.Q]') & [tempWpbQ.Q]' > 0;
%         display(sum(idQ))
%         logQ = log10([tempWpbQ(idQ,1).Q]');
%         [fQ, xQ]=ecdf(logQ);
%     else
%         logQ = log10([tempWpb(idQ,1).Q]');
%         [fQ, xQ]=ecdf(logQ);
%     end
%     
%     
%     
%     %%% Depth ECDF
%     idD = ~isnan([tempWpb.depth]') & [tempWpb.depth]' > 0;
%     display(sum(idD))
%     logD = log10([tempWpb(idD,1).depth]');
%     [fD, xD]=ecdf(logD);
%     
%     %%% Screen length ECDF
%     tempSL = [tempWpb.bot]' - [tempWpb.top]';
%     idS = ~isnan(tempSL) & tempSL < 2000  & tempSL > 0;
%     display(sum(idS))
%     logS = log10([tempWpb(idS,1).bot]' - [tempWpb(idS,1).top]');
%     [fS, xS]=ecdf(logS);
%     
%     clear Stats
%     Stats.fQ = fQ; Stats.lowQ = 0.15;
%     Stats.xQ = xQ; Stats.uppQ = 0.85;
%     Stats.fD = fD; Stats.lowD = 0;
%     Stats.xD = xD; Stats.uppD = 1;
%     Stats.fS = fS; Stats.lowS = 0;
%     Stats.xS = xS; Stats.uppS = 1;
%     
%     for jj = 1:length(tempWpb_mod)
%         [ tempWpb_mod(jj,1).Q, ...
%           tempWpb_mod(jj,1).depth, ...
%           tempWpb_mod(jj,1).SL ] = AssignQDS(Stats, ...
%                                              tempWpb_mod(jj,1).Q, ...
%                                              tempWpb_mod(jj,1).depth, ...
%                                              tempWpb_mod(jj,1).bot - tempWpb_mod(jj,1).top);
%     end
    
    ipub = find([tempW_mod.type]'==1);
    iag = find([tempW_mod.type]'==0);
    tempWpb_mod = tempW_mod(ipub,1);
    tempWag_mod = tempW_mod(iag,1);
    
    
    WELLS4CVHM(ii,1).farm_id = Farms_poly(ii,1).dwr_sbrgns;
    WELLS4CVHM(ii,1).UrbanPump = Farms_poly(ii,1).UrbanPump;
    WELLS4CVHM(ii,1).AgPump = Farms_poly(ii,1).AgPump;
    WELLS4CVHM(ii,1).Pubwells = tempWpb_mod;
    WELLS4CVHM(ii,1).Agwells = tempWag_mod;
    WELLS4CVHM(ii,1).Npb = length(tempWpb_mod);
    WELLS4CVHM(ii,1).Nag = length(tempWag_mod);
end
WELLS4CVHM(17,:) = [];
%
%% Scale the pumping rates to match the desired pumping
% RCH : 25715545.162923 m^3/day
% STRM: 899912.211251 m^3/day
Total_pump_vol = 25715545.162923 + 899912.211251; %m^3/day
Min_pumping = 400;
init_vol = abs(sum([WELLS4CVHM.UrbanPump]')) + abs(sum([WELLS4CVHM.AgPump]'));
for ii = 1:length(WELLS4CVHM)
    
    urb_perc = abs(WELLS4CVHM(ii,1).UrbanPump) / init_vol;
    ag_perc = abs(WELLS4CVHM(ii,1).AgPump) / init_vol;
    
    urb_act = Total_pump_vol*urb_perc;
    ag_act = Total_pump_vol*ag_perc;
    
    %%% Public wells
    Q_init = [WELLS4CVHM(ii,1).Pubwells.Q]';
    Q_perc = Q_init/sum(Q_init);
    Q_temp = urb_act*Q_perc;
    c = min(Q_temp);
    if c < Min_pumping
        cnt = 1;
        while 1
            Q_init = [WELLS4CVHM(ii,1).Pubwells.Q]';
            [c, d] = sort(Q_init);
            Q_init(d(1:cnt,1),:) = [];
            Q_perc = Q_init/sum(Q_init);
            Q_temp = urb_act*Q_perc;
            c = min(Q_temp);
            if c > Min_pumping
                break;
            else
                cnt = cnt + 1;
            end
        end
        WELLS4CVHM(ii,1).Pubwells(d(1:cnt,1),:) = [];
    end
   
    Q_init = [WELLS4CVHM(ii,1).Pubwells.Q]';
    Q_perc = Q_init/sum(Q_init);
    for jj = 1:length(Q_perc)
        WELLS4CVHM(ii,1).Pubwells(jj,1).Q_act = urb_act*Q_perc(jj);
    end
    
    %%% Ag wells
    Q_init = [WELLS4CVHM(ii,1).Agwells.Q]';
    Q_perc = Q_init/sum(Q_init);
    Q_temp = ag_act*Q_perc;
    c = min(Q_temp);
    if c < Min_pumping
        cnt = 1;
        while 1
            Q_init = [WELLS4CVHM(ii,1).Agwells.Q]';
            [c, d] = sort(Q_init);
            Q_init(d(1:cnt,1),:) = [];
            Q_perc = Q_init/sum(Q_init);
            Q_temp = ag_act*Q_perc;
            c = min(Q_temp);
            if c > Min_pumping
                break;
            else
                cnt = cnt + 1;
            end
        end
        WELLS4CVHM(ii,1).Agwells(d(1:cnt,1),:) = [];
    end
    
    Q_init = [WELLS4CVHM(ii,1).Agwells.Q]';
    Q_perc = Q_init/sum(Q_init);
    for jj = 1:length(Q_perc)
        WELLS4CVHM(ii,1).Agwells(jj,1).Q_act = ag_act*Q_perc(jj);
    end
end
%% Check total volume
Q_ag = 0;
Q_pb = 0;
for ii = 1:length(WELLS4CVHM)
    Q_ag = Q_ag + sum([WELLS4CVHM(ii,1).Agwells.Q_act]);
    Q_pb = Q_pb + sum([WELLS4CVHM(ii,1).Pubwells.Q_act]);
end
%% Put it all in one variable table
WW = nan(20000,7);
cnt = 1;
for ii = 1:size(WELLS4CVHM,1)
    for jj = 1:length(WELLS4CVHM(ii,1).Agwells)
        WW(cnt,:) = [WELLS4CVHM(ii,1).Agwells(jj,1).X ...
                     WELLS4CVHM(ii,1).Agwells(jj,1).Y ...
                     WELLS4CVHM(ii,1).Agwells(jj,1).depth ...
                     WELLS4CVHM(ii,1).Agwells(jj,1).SL ...
                     WELLS4CVHM(ii,1).Agwells(jj,1).Q_act ...
                     ii jj];
        cnt = cnt + 1;
    end
    
    for jj = 1:length(WELLS4CVHM(ii,1).Pubwells)
        WW(cnt,:) = [WELLS4CVHM(ii,1).Pubwells(jj,1).X ...
                     WELLS4CVHM(ii,1).Pubwells(jj,1).Y ...
                     WELLS4CVHM(ii,1).Pubwells(jj,1).depth ...
                     WELLS4CVHM(ii,1).Pubwells(jj,1).SL ...
                     WELLS4CVHM(ii,1).Pubwells(jj,1).Q_act ...
                     ii jj];
        cnt = cnt + 1;
    end
end
WW(cnt:end,:) = [];
WW_bck = WW; % make a backup table
%% Randomize well locations
WW = WW_bck;
for ii = 1:size(WW,1)
    ii
    c = sqrt((WW(:,1) - WW(ii,1)).^2 + (WW(:,2) - WW(ii,2)).^2);
    id = find(c < 400);
    if length(id) > 1
        xm = WW(ii,1);
        ym = WW(ii,2);
        for jj = 1:length(id)
            cnt = 0;
            thres = 400;
            ifarm = WELLS4CVHM(WW(id(jj),6),1).farm_id;% This is the id of the farm
            ifarm = find([Farms_poly.dwr_sbrgns]' == ifarm); % This is the index of the farm in Farm_poly
            
            while 1
                while 1 % make sure the generated well is inside the farm
                    xr = xm - 600 + (xm + 600 - (xm - 600))*rand;
                    yr = ym - 600 + (ym + 600 - (ym - 600))*rand;
                    in = inpolygon(xr, yr, Farms_poly(ifarm,1).X, Farms_poly(ifarm,1).Y) & ...
                         inpolygon(xr, yr, cvhmDomain.X, cvhmDomain.Y);
                    if in
                        break;
                    end
                end
                if jj == 1
                    WW(id(jj),1) = xr;
                    WW(id(jj),2) = yr;
                    break;
                else
                    dst = sqrt((WW(id(1:jj-1),1) - xr).^2 + (WW(id(1:jj-1),2)- yr).^2);
                    if min(dst) > thres
                        WW(id(jj),1) = xr;
                        WW(id(jj),2) = yr;
                        break;
                    end
                end
                cnt = cnt + 1;
                if cnt > 10
                    thres = thres - 10;
                    thres
                    cnt = 0;
                end
            end
        end
    end
end
%% check minimum distances
for ii = 1:size(WW,1)
    tii = setdiff([1:size(WW,1)]',ii);
    mindst(ii,1) = min(sqrt((WW(tii,1) - WW(ii,1)).^2 + (WW(tii,2) - WW(ii,2)).^2));
end
%% relocate wells that are closer than 400 m
[c d] = sort(mindst);
for ii = 1:length(c)
    ii
    if c(ii) > 400
        break;
    end
    tii = setdiff([1:size(WW,1)]',d(ii));
    xc = WW(d(ii),1);
    yc = WW(d(ii),2);
    ifarm = WELLS4CVHM(WW(d(ii),6),1).farm_id;% This is the id of the farm
    ifarm = find([Farms_poly.dwr_sbrgns]' == ifarm); % This is the index of the farm in Farm_poly
    radius = 10;
    cnt = 0;
    while 1
        % pick a random direction and a random radius
        while 1
            r_rad = radius*rand;
            tr = 2*pi*rand;
            xr = xc + r_rad*cos(tr);
            yr = yc + r_rad*sin(tr);
            in = inpolygon(xr, yr, Farms_poly(ifarm,1).X, Farms_poly(ifarm,1).Y) & ...
                 inpolygon(xr, yr, cvhmDomain.X, cvhmDomain.Y);
            if in
                break;
            end
        end
        
        dst = min(sqrt((WW(tii,1) - xr).^2 + (WW(tii,2) - yr).^2));
        if dst > 400
            WW(d(ii),1) = xr;
            WW(d(ii),2) = yr;
            break;
        end
        cnt = cnt + 1;
        if cnt > 100
            radius = radius + 10;
            cnt = 0;
        end
    end
end
%% Correct those that their Screen Lengths are greater or equal to depth
id = find(WW(:,3) <= WW(:,4));
WW(id,3) = WW(id,4) + 30;
%% Assign top and bottom according to the artificial elevation and CVHM bottom
topelev = read_Scattered('CVHM_top_elev_per2.npsat',2);
Ftop = scatteredInterpolant(topelev.p(:,1), topelev.p(:,2), topelev.v, 'natural');

% abd the bottom of the aquifer
botelev = read_Scattered('CVHM_Bot_elev.npsat',2);
Fbot = scatteredInterpolant(botelev.p(:,1), botelev.p(:,2), botelev.v, 'natural');
%
%%
load('Random_wells.mat', 'Wells_Gen')
WW = [[Wells_Gen.X]' [Wells_Gen.Y]' [Wells_Gen.D]' [Wells_Gen.SL]' [Wells_Gen.Q]'];

%}
%%
cnvrs = 0.3048;
for ii = 1:size(WW,1)
    Wt = Ftop(WW(ii,1), WW(ii,2));
    Bs = Fbot(WW(ii,1), WW(ii,2));
    aquif_height = Wt - Bs;
    depth = WW(ii,3)*cnvrs;
    screen = WW(ii,4)*cnvrs;
    if aquif_height < depth + screen
        if depth > screen
            botw = Bs + 30;
            topw = botw + screen;
            if topw > Wt
                topw = Wt-30;
            end
        else
            topw = Wt - depth;
            botw = Bs + 30;
            if topw < botw
                topw = Wt-30;
            end
        end
    else
        topw = Wt - depth;
        botw = topw - screen;
    end
    
    
%     if aquif_height < WW(ii,3)*cnvrs
%         if aquif_height > 30
%             botw = Bs + 10;
%         else
%             wtf = true;
%         end
%     else
%         botw = Wt - WW(ii,3)*cnvrs;
%     end
%     
%     if botw + WW(ii,4)*cnvrs > Wt
%         if aquif_height > 30
%             topw = Wt - 10;
%         else
%             wtf = true;
%         end
%     else
%         topw = botw+WW(ii,4)*cnvrs;
%     end
    
    WW(ii,8) = topw;
    WW(ii,9) = botw;
    WW(ii,10) = Wt;
    WW(ii,11) = Bs;
end

%% Correct the wells with 
id1 = find(WW(:,9) < WW(:,11))
id2 = find(WW(:,8) < WW(:,9))
id3 = find(WW(:,10) < WW(:,8))
%% Print to file

fid = fopen('CVHM_wells_per2.npsat','w');
fprintf(fid, '%d\n', size(WW, 1));
fprintf(fid, '%f %f %f %f -%f\n', [WW(:,1) WW(:,2) WW(:,8) WW(:,9) WW(:,5)]');
fclose(fid);

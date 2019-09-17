function Wells_Gen = CVHMgenerateWells(avWells, opt)

if ~isfield(opt, 'LU')
    opt.LU = 'LU_2000';
end

if ~isfield(opt, 'do_plot')
    opt.do_plot = false;
end


% Find the volume of water that is used for Urban and Ag for each farm
% This contains the farms outline. 22 polygons in total inlcudeing the out
% of study area polygon
load('gis_data/FARMS_poly_shp'); %Farms_poly = shaperead('gis_data/FARMS_poly.shp');
% This has all grid cells a features with attributes the row columns and land use 
load('gis_data/FMP_shp'); %Farms = shaperead('gis_data/FMP.shp');
load('../CVHM_DATA/WellData/CVwells_30y_proj_shp.mat')
load('gis_data/WellGenerationSHP','cvhmDomain');
load('gis_data/CVHM_Mesh_outline_modif_shp.mat')

cell_farm_ids = [Farms.dwr_sbrgns]';
ROWfrm = [Farms.ROW]';
COLfrm = [Farms.COLUMN_]';
LUfarm = [Farms.(opt.LU)]';
for ii = 1:length(Farms_poly)
    id_farm = Farms_poly(ii,1).dwr_sbrgns;
    farm_i_cells = find(cell_farm_ids == id_farm);
    farm_ind = sub2ind([441 98], ROWfrm(farm_i_cells), COLfrm(farm_i_cells));
    farm_i_lu = LUfarm(farm_i_cells);
    wRates = avWells(farm_ind);
    Farms_poly(ii,1).UrbanPump = sum(wRates(farm_i_lu == 2));
    Farms_poly(ii,1).AgPump = sum(wRates(farm_i_lu ~= 2));
end

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
    figure(1);title(['Farm: ' num2str(Farms_poly(ii,1).dwr_sbrgns) ' (' num2str(ii) '/' num2str(length(Farms_poly)) ')']);
    
    % Caclulate spatial distribution of wells
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
    
    % Pumping distribution (This is different for Ag and Urban wells unless there are
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
    
    % The data in Rich well data set are in gpm while the units of CVHM are
    % in m^3/day. To convert the pumping to m^/day Q = Q*5.451.
    % However this correpond to the design pumping. The argument is that
    % each year only six months operate with this rate therefore we divide
    % the pumping by 6 months Q = Q*5.451/6. 
    % As we dont want rates lower than 400 m^3/day this would be
    % log10(400*(6/5.451)) = 2.6437 gpm.
    % 10^2.6437*5.451/6
    
    if isempty(find(fag < 2.6437, 1, 'last'))
        Ql = xag(1);
    else
        Ql = xag(find(fag < 2.6437, 1, 'last')+1);
    end
    while Qag_farm < abs(Farms_poly(ii,1).AgPump)
        % generate a valid well location
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
        if abs(Farms_poly(ii,1).AgPump) - (Qag_farm + 10^Qr*5.451/6) < 400
           Qr = log10((abs(Farms_poly(ii,1).AgPump) - Qag_farm)*6/5.451); 
        end
        
        [ Qw, Dw, Sw ] = AssignQDS_v2(QD_pdf, DS_pdf, Qr, nan, nan);
        Qag_farm = Qag_farm + 10^Qw*5.451/6; % m^3/day for six months
        
        WX(cnt,1) = xr; WY(cnt,1) = yr;
        if opt.do_plot
            figure(2);plot(xr,yr,'.b')
            drawnow
            figure(5);plot(Qw,Dw,'.b')
            drawnow
            figure(6);plot(Dw,Sw,'.b')
            drawnow
        end
        figure(4);title([num2str(abs(Farms_poly(ii,1).AgPump) - Qag_farm) ' ' num2str(cnt)]);
        drawnow
        
        Wells_Gen(cnt,1).X = xr;
        Wells_Gen(cnt,1).Y = yr;
        Wells_Gen(cnt,1).Q = 10^Qw*5.451/6;
        Wells_Gen(cnt,1).D = 10^Dw;
        Wells_Gen(cnt,1).SL = 10^Sw;
        Wells_Gen(cnt,1).type = 0;
        cnt = cnt + 1;
    end
    
    if isempty(find(fpb < 2.6437, 1, 'last'))
        Ql = xpb(1);
    else
        Ql = xpb(find(fpb < 2.6437, 1, 'last')+1);
    end
    
    while Qpb_farm < abs(Farms_poly(ii,1).UrbanPump)
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
        if abs(Farms_poly(ii,1).UrbanPump) - (Qpb_farm + 10^Qr*5.451/6) < 400
            Qr = log10((abs(Farms_poly(ii,1).UrbanPump) - Qpb_farm)*6/5.451);
        end
        
        [ Qw, Dw, Sw ] = AssignQDS_v2(QD_pdf, DS_pdf, Qr, nan, nan);
        Qpb_farm = Qpb_farm + 10^Qw*5.451/6;
        WX(cnt,1) = xr; WY(cnt,1) = yr;
        if opt.do_plot
            figure(3);plot(xr,yr,'.b')
            drawnow
            figure(5);plot(Qw,Dw,'.b')
            drawnow
            figure(6);plot(Dw,Sw,'.b')
            drawnow
        end
        figure(4);title([num2str(abs(Farms_poly(ii,1).UrbanPump) - Qpb_farm) ' ' num2str(cnt)] );
        drawnow
        Wells_Gen(cnt,1).X = xr;
        Wells_Gen(cnt,1).Y = yr;
        Wells_Gen(cnt,1).Q = 10^Qw*5.451/6;
        Wells_Gen(cnt,1).D = 10^Dw;
        Wells_Gen(cnt,1).SL = 10^Sw;
        Wells_Gen(cnt,1).type = 1;
        cnt = cnt + 1;
    end
end

%{
WW = [[Wells_Gen.X]' [Wells_Gen.Y]' [Wells_Gen.D]' [Wells_Gen.SL]' [Wells_Gen.Q]'];

topelev = read_Scattered([opt.simFolder filesep opt.prefix '_' timestring  '_Top.npsat'], 2);
Ftop = scatteredInterpolant(topelev.p(:,1), topelev.p(:,2), topelev.v, 'natural');

% abd the bottom of the aquifer
botelev = read_Scattered('CVHM_Bot_elev.npsat',2);
Fbot = scatteredInterpolant(botelev.p(:,1), botelev.p(:,2), botelev.v, 'natural');



% Print wells to file
fid = fopen([opt.simFolder filesep opt.prefix '_' timestring  'InitWells.npsat'],'w');

fprintf(fid, '%d\n', size(WW, 1));
fprintf(fid, '%f %f %f %f -%f\n', [WW(:,1) WW(:,2) WW(:,8) WW(:,9) WW(:,5)]');
fclose(fid);
%}
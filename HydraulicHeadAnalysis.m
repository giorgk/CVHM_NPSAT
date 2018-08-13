%% Create a list of lines of the domain outline.
% For each line find in which cell they correspond
%{
mesh_outline = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/BAS_active_outline.shp');
bas = shaperead('/home/giorgk/Documents/UCDAVIS/CVHM_NPSAT/gis_data/BAS_active.shp');
%
%% make a unique list of outline points
%
PNTS = nan(5000,2);
cnt = 1;
for ii = 1:size(mesh_outline,1)
    [Xs, Ys] = polysplit(mesh_outline(ii,1).X, mesh_outline(ii,1).Y);
    for jj = 1:size(Xs,1)
        for k = 1:length(Xs{jj,1})
            if cnt == 1
                PNTS(cnt,:) = [Xs{jj,1}(k) Ys{jj,1}(k)];
                cnt = cnt + 1;
            else
                dst = sqrt((Xs{jj,1}(k) - PNTS(1:cnt-1, 1)).^2 + (Ys{jj,1}(k) - PNTS(1:cnt-1, 2)).^2);
                if min(dst) >= 0.01
                    PNTS(cnt,:) = [Xs{jj,1}(k) Ys{jj,1}(k)];
                    PNTS(cnt,:) = [Xs{jj,1}(k) Ys{jj,1}(k)];
                    cnt = cnt + 1;
                end
            end
        end
    end
end
PNTS(cnt:end,:) = [];
%
%% make a list of unique outline segments
%
LNS = [];
for ii = 1:size(mesh_outline,1)
    [Xs, Ys] = polysplit(mesh_outline(ii,1).X, mesh_outline(ii,1).Y);
    for jj = 1:size(Xs,1)
        for k = 1:length(Xs{jj,1})-1
            pa = [Xs{jj,1}(k) Ys{jj,1}(k)];
            pb = [Xs{jj,1}(k+1) Ys{jj,1}(k+1)];
            dst = sqrt((pa(1) - PNTS(:, 1)).^2 + (pa(2) - PNTS(:, 2)).^2);
            ida = find(dst < 0.01);
            dst = sqrt((pb(1) - PNTS(:, 1)).^2 + (pb(2) - PNTS(:, 2)).^2);
            idb = find(dst < 0.01);
            if isempty(ida) || isempty(idb)
               error('wtf');
            else
                % check if we have already add this line
                idmn = min([ida idb]);
                idmx = max([ida idb]);
                if isempty(LNS)
                    LNS = [LNS; idmn idmx];
                else
                    idln = find(LNS(:,1) == idmn & LNS(:,2) == idmx);
                    if isempty(idln)
                        LNS = [LNS; idmn idmx];
                    end
                end
            end
        end
    end
end
%
%% loop through the edges of the bas and set the ij ids for those that 
% match the outline segments
LNS_IJ = nan(size(LNS,1),2);
for ii = 1:size(bas,1)
    [Xs, Ys] = polysplit(bas(ii,1).X, bas(ii,1).Y);
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
                    idmn = min([ida idb]);
                    idmx = max([ida idb]);
                    idln = find(LNS(:,1) == idmn & LNS(:,2) == idmx);
                    if ~isempty(idln)
                        LNS_IJ(idln,:) = [bas(ii,1).ROW bas(ii,1).COLUMN_];
                    end
                end
            end
        end
    end
end

%% Assign head values and standard deviations
% For the analysis load the file HEADS.mat
% This is not included in the repository because the size is about 224 MB
%
load('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/ClaudiaRun/HEADS.mat')
clear TOPHEAD
% For each time step compute the top most head value
for it = 1:size(HEADS,1)
    it
    TOPHEAD{it,1} = nan(size(HEADS{1,1},1), size(HEADS{1,1},2));
    for ii = 1:size(HEADS{1,1},1)
        for jj = 1:size(HEADS{1,1},2)
            for kk = 1:size(HEADS{1,1},3)%:-1:1
                if HEADS{it,1}(ii,jj,kk) ~= HEADS{1,1}(1,1,1)
                    TOPHEAD{it,1}(ii,jj) = HEADS{it,1}(ii,jj,kk);
                end
            end
        end
    end
end
%
%% Calculate Standard deviation
% get the istart, iend variables from the GWaterBudget_calc script
head_std = nan(size(HEADS{1,1},1), size(HEADS{1,1},2));
head_Av = nan(size(HEADS{1,1},1), size(HEADS{1,1},2));
timespan = istart:1:iend;
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
%%
temp = head_std;
%temp(head_std>2) = nan;
surf(temp,'edgecolor','none')
view(360,-90)
%
%% assign values to segments and create BC shapefile
%
clear S
timespan = istart:1:iend;
for ii = 1:size(LNS,1)
    S(ii,1).Geometry = 'PolyLine';
    S(ii,1).X = [PNTS(LNS(ii,:),1)' nan];
    S(ii,1).Y = [PNTS(LNS(ii,:),2)' nan];
    S(ii,1).ROW = LNS_IJ(ii,1);
    S(ii,1).COL = LNS_IJ(ii,2);
    temp = nan(length(timespan),1);
    for kk = 1:length(timespan)
        temp(kk,1) = TOPHEAD{timespan(kk),1}(S(ii,1).ROW, S(ii,1).COL);
    end
    S(ii,1).Hav = mean(temp);
    S(ii,1).Hstd = std(temp);
end
%%
shapewrite(S,'gis_data/CVHM_outline_head')
%% Add the bottom info
bot = read_Scattered('CVHM_Bot_elev.npsat', 2);
Fbot = scatteredInterpolant(bot.p(:,1), bot.p(:,2), bot.v(:,1));
%% interpolate the head points
for ii = 1:size(S, 1)
    xc = mean(S(ii,1).X(1:end-1));
    yc = mean(S(ii,1).Y(1:end-1));
    S(ii,1).bot = Fbot(xc, yc);
end
%find([S.bot]' > [S.Hav]')
%
%% Interpolate head values on nodes
Pnt_hd = nan(size(PNTS, 1),1);
for ii = 1:size(PNTS, 1)
    % find how many segments have this node in common
    id = find(LNS(:,1) == ii | LNS(:,2) == ii);
    hd_nd = nan(length(id), 1);
    for jj = 1:length(id)
        hd_nd(jj) = head_Av(LNS_IJ(id(jj),1),LNS_IJ(id(jj),2));
    end
    Pnt_hd(ii,1) = mean(hd_nd);
end
%}
%% write the file
% I have first to implement the line boundary
fid = fopen('CVHM_BC.npsat','w');
fprintf(fid, '%d\n', size(LNS,1));
for ii = 1:size(LNS,1)
    fprintf(fid, 'EDGETOP 2 BC_files/%s\n', ['cvhm_linebc_' num2str(ii) '.npsat'] );
    fprintf(fid, '%f %f\n', PNTS(LNS(ii,1),:)); %Pnt_hd(LNS(ii,1),1)
    fprintf(fid, '%f %f\n', PNTS(LNS(ii,2),:)); %Pnt_hd(LNS(ii,2),1)
    
    fid1 = fopen(['BC_files/cvhm_linebc_' num2str(ii) '.npsat'], 'w');
    fprintf(fid1, 'SCATTERED\n');
    fprintf(fid1, 'VERT\n');
    fprintf(fid1, 'STRATIFIED\n');
    fprintf(fid1, '%d %d\n', [2 1]);
    fprintf(fid1, '%f %f\n', [0 Pnt_hd(LNS(ii,1),1)]);
    fprintf(fid1, '%f %f\n', [pdist(PNTS(LNS(ii,:),:)) Pnt_hd(LNS(ii,2),1)]);
    fclose(fid1);
    
end
fclose(fid);

%% Use the average of Head for the simulation period as top function.
% Run the 2 sections which read and average the hydraulic head
% Run the first 3 sections of the PrepareGeometryData.m to import and
% process the bas variable
% Loop through the modflow grid cells and find the center coordinate of the
% cell in the real system
R = [bas.ROW]';
C = [bas.COLUMN_]';
xy_top = [];
for ii = 1:size(head_Av, 1)
    for jj = 1:size(head_Av, 2)
        id = find(R == ii & C == jj);
        if isempty(id)
            continue;
        elseif length(id) == 1
            xc = mean(bas(id,1).X(1:end-1));
            yc = mean(bas(id,1).Y(1:end-1));
            xy_top = [xy_top; xc yc head_Av(ii,jj)];
        else
            error('Thats so wrong');
        end
    end
end
%% 
% add the buffer nodes.
% run the section from PrepareGeometryData.m which creates the buff_pnt
% variable
% create an interpolant with the existing data
Ftop = scatteredInterpolant(xy_top(:,1), xy_top(:,2), xy_top(:,3));
Ftop.Method = 'nearest';
Ftop.ExtrapolationMethod = 'nearest';
buf_val = Ftop(buff_pnt(:,1), buff_pnt(:,2));
xy_top = [xy_top; buff_pnt buf_val];
%% write the rch file
writeScatteredData('CVHM_top_elev.npsat', struct('PDIM',2,'TYPE','HOR','MODE','SIMPLE'), xy_top);

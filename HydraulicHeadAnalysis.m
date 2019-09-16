%% Create a list of lines of the domain outline.
% For each line find in which cell they correspond
%
mesh_outline = shaperead('gis_data/CVHM_Mesh_outline_modif.shp');
mesh = shaperead('gis_data/CVHM_mesh.shp');
%% make a list of nodes that describe the boundary of the domain
% mesh_outline has only one polygon and we split it bellow
[Xs, Ys] = polysplit(mesh_outline(1,1).X, mesh_outline(1,1).Y);
% only the first polygon described the mesh outline
PNTS = [];
for ii = 1:length(Xs{1,1})
    PNTS = [PNTS; Xs{1,1}(ii) Ys{1,1}(ii)];
end
%% loop through the edges of the bas and set the ij ids for those that 
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
%% or load this
load('HeadAnalysisData.mat')
%% Assign head values and standard deviations
% For the analysis load the file HEADS.mat
% This is not included in the repository because the size is about 224 MB
%
%load('/home/giorgk/Documents/UCDAVIS/CVHM_DATA/ClaudiaRun/HEADS.mat')
load('/media/giorgk/DATA/giorgk/Documents/CVHM_DATA/ClaudiaRun/HEADS.mat')
clear TOPHEAD toplay
% For each time step compute the top most head value
for it = 1:size(HEADS,1)
    it
    TOPHEAD{it,1} = nan(size(HEADS{1,1},1), size(HEADS{1,1},2));
    toplay{it,1} = nan(size(HEADS{1,1},1), size(HEADS{1,1},2));
    for ii = 1:size(HEADS{1,1},1)
        for jj = 1:size(HEADS{1,1},2)
            for kk = 1:size(HEADS{1,1},3)%:-1:1
                if HEADS{it,1}(ii,jj,kk) ~= HEADS{1,1}(1,1,1)
                    TOPHEAD{it,1}(ii,jj) = HEADS{it,1}(ii,jj,kk);
                    toplay{it,1}(ii, jj) = kk;
                    break
                end
            end
        end
    end
end
%
%% Calculate Standard deviation
% get the istart2, iend2 variables from the GWaterBudget_calc script
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
%%
temp = head_std;
%temp(head_std>2) = nan;
surf(temp,'edgecolor','none')
view(360,-90)
%
%% Assign head values to segments
timespan = istart2:1:iend2;
LNS_IJ = [LNS_IJ nan(size(LNS_IJ,1),2)];
for ii = 1:size(LNS_IJ, 1) - 1
    r = LNS_IJ(ii,1);
    c = LNS_IJ(ii,2);
    temp = nan(length(timespan),1);
    for kk = 1:length(timespan)
        temp(kk,1) = TOPHEAD{timespan(kk),1}(r, c);
    end
    LNS_IJ(ii,3:4) = [mean(temp) std(temp)];
end

% % % %% assign values to segments and create BC shapefile
% % % %
% % % clear S
% % % timespan = istart:1:iend;
% % % for ii = 1:size(LNS,1)
% % %     S(ii,1).Geometry = 'PolyLine';
% % %     S(ii,1).X = [PNTS(LNS(ii,:),1)' nan];
% % %     S(ii,1).Y = [PNTS(LNS(ii,:),2)' nan];
% % %     S(ii,1).ROW = LNS_IJ(ii,1);
% % %     S(ii,1).COL = LNS_IJ(ii,2);
% % %     temp = nan(length(timespan),1);
% % %     for kk = 1:length(timespan)
% % %         temp(kk,1) = TOPHEAD{timespan(kk),1}(S(ii,1).ROW, S(ii,1).COL);
% % %     end
% % %     S(ii,1).Hav = mean(temp);
% % %     S(ii,1).Hstd = std(temp);
% % % end
% % % %%
% % % shapewrite(S,'gis_data/CVHM_outline_head')
%% Add the bottom info
bot = read_Scattered('CVHM_Bot_elev.npsat', 2);
Fbot = scatteredInterpolant(bot.p(:,1), bot.p(:,2), bot.v(:,1));
% % % %% interpolate the head points
% % % for ii = 1:size(S, 1)
% % %     xc = mean(S(ii,1).X(1:end-1));
% % %     yc = mean(S(ii,1).Y(1:end-1));
% % %     S(ii,1).bot = Fbot(xc, yc);
% % % end
% % % %find([S.bot]' > [S.Hav]')
% % % %
%% Interpolate the head and standard deviations values on the nodes
PNTS = [PNTS nan(size(PNTS,1),2)];
for ii = 1:size(PNTS,1)-1
    if ii == 1
        % if this is the first we average the first and last segments
        PNTS(ii,3:4) = [(LNS_IJ(1,3) + LNS_IJ(end-1,3) )/2 (LNS_IJ(1,4) + LNS_IJ(end-1,4) )/2];
        PNTS(end,3:4) = [(LNS_IJ(1,3) + LNS_IJ(end-1,3) )/2 (LNS_IJ(1,4) + LNS_IJ(end-1,4) )/2];
    else
        PNTS(ii,3:4) = [(LNS_IJ(ii,3) + LNS_IJ(ii-1,3) )/2 (LNS_IJ(ii,4) + LNS_IJ(ii-1,4) )/2];
    end
end
%% Define the number of boundary functions based on the head standard deviation
std_Htol = 2.5;
cnt_bnd = 0;
clear BND_LINES
next_line = 1;
for ii = 1:size(PNTS,1)
    if PNTS(ii,4) < std_Htol
        if next_line == 1
            cnt_bnd = cnt_bnd + 1;
            BND_LINES{cnt_bnd,1} = [];
            next_line = 0;
        end
        BND_LINES{cnt_bnd,1} = [BND_LINES{cnt_bnd,1}; ii];
    else
        next_line = 1;
    end
end
% remove the segments with length 1
dlt = [];
for ii = 1:length(BND_LINES)
    if length(BND_LINES{ii,1}) == 1
        dlt = [dlt; ii];
    end
end
BND_LINES(dlt,:) = [];
%% write the file using boundary function
fid = fopen('CVHM_BC.npsat','w');
fprintf(fid, '%d\n', length(BND_LINES));
for ii = 1:length(BND_LINES)
    fprintf(fid, 'EDGETOP 0 BC_files/%s\n', ['bndfnc_' num2str(ii) '.npsat']);
    
    fid1 = fopen(['BC_files/bndfnc_' num2str(ii) '.npsat'], 'w');
    fprintf(fid1, 'BOUNDARY_LINE\n');
    fprintf(fid1, '%d %d %f\n', [length(BND_LINES{ii,1}) 1 1]); % Npnts Ndata tolerance
    fprintf(fid1, '%f %f %f\n', PNTS(BND_LINES{ii,1},1:3)');
    fclose(fid1);
end
fclose(fid);

%% Use the average of Head for the simulation period as top function.
% Run the 2 sections which read and average the hydraulic head
% Run the first 3 sections of the PrepareGeometryData.m to import and
% process the bas variable
% Loop through the modflow grid cells and find the center coordinate of the
% cell in the real system
R = [bas_active.ROW]';
C = [bas_active.COLUMN_]';
xy_top = [];
for ii = 1:size(head_Av, 1)
    for jj = 1:size(head_Av, 2)
        id = find(R == ii & C == jj);
        if isempty(id)
            continue;
        elseif length(id) == 1
            xc = mean(bas_active(id,1).X(1:end-1));
            yc = mean(bas_active(id,1).Y(1:end-1));
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
%% write the top elevation file
writeScatteredData('CVHM_top_elev_per2.npsat', struct('PDIM',2,'TYPE','HOR','MODE','SIMPLE'), xy_top);
%}
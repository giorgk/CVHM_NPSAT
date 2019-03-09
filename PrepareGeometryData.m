%% =================== MESH FILE based on CVHM grid =======================
% From the complete CVHM database
% https://ca.water.usgs.gov/projects/central-valley/cvhm-database.html
% extract and load the BAS file 
%
bas = shaperead('gis_data/BAS6');
%% Isolate the active grid.
for ii = 1:size(bas,1)
    active = 0;
    for jj = [1 2 3 45 6 7 8 9 10]
        active = active + bas(ii,1).(['lay' num2str(jj) '_act']);
    end
    if active > 0
        bas(ii,1).grid_active = 1;
    else
        bas(ii,1).grid_active = 0;
    end
end
%% delete non active cells
id = find([bas.grid_active]'== 0);
bas(id,:) = [];
%%
shapewrite(bas,'gis_data/BAS_active');
%
%% Make a list of unique nodes and elements
msh_nd = [];
msh_el = [];
msh_basIJ = [];
for ii = 1:size(bas,1)
    ii
    el = [];
    for jj = 1:4
       px = bas(ii,1).X(jj);
       py = bas(ii,1).Y(jj);
       % test if this node already exists
       if isempty(msh_nd)
           msh_nd = [px py];
           el = 1;
       else
           dst = sqrt((px - msh_nd(:,1)).^2 + (py - msh_nd(:,2)).^2);
           id = find(dst < 0.1);
           if isempty(id)
               msh_nd = [msh_nd;px py];
               el = [el size(msh_nd,1)];
           else
               if length(id) == 1
                   el = [el id];
               else
                   wtf = true;
               end
           end
       end
    end
    msh_el = [msh_el; el];
    msh_basIJ = [msh_basIJ; bas(ii,1).ROW bas(ii,1).COLUMN_];
end
%% Create a shapefile for the mesh
for ii = 1:size(msh_el,1)
    xx = msh_nd(msh_el(ii,[1 2 3 4 1]),1)';
    yy = msh_nd(msh_el(ii,[1 2 3 4 1]),2)';
    S(ii,1).Geometry = 'Polygon';
    S(ii,1).BoundingBox = [min(xx) min(yy); max(xx) max(yy)];
    S(ii,1).X = [xx nan];
    S(ii,1).Y = [yy nan];
    S(ii,1).ID = ii;
    S(ii,1).R = msh_basIJ(ii,1);
    S(ii,1).C = msh_basIJ(ii,2);
    %S(ii,1).ND1 = msh_el(ii,1);
    %S(ii,1).ND2 = msh_el(ii,2);
    %S(ii,1).ND3 = msh_el(ii,3);
    %S(ii,1).ND4 = msh_el(ii,4);
end
shapewrite(S,'gis_data/CVHM_mesh');
%% load the modified mesh shapefile
S = shaperead('gis_data/CVHM_mesh');
%% Create unique nodes and elements
msh_nd = [];
msh_el = [];
for ii = 1:length(S)
    el_ids = [];
    for jj = 1:4
        Xnd = S(ii,1).X(jj);
        Ynd = S(ii,1).Y(jj);
        if isempty(msh_nd)
            msh_nd = [msh_nd; Xnd Ynd];
            el_ids = [el_ids 1];
        else
            dst = sqrt((msh_nd(:,1) - Xnd).^2 + (msh_nd(:,2) - Ynd).^2);
            id = find(dst < 0.01);
            if isempty(id)
                msh_nd = [msh_nd; Xnd Ynd];
                el_ids = [el_ids size(msh_nd,1)];
            else
                if length(id) > 1
                    error('More than one nodes are found in the same location');
                end
                el_ids = [el_ids id];
            end 
        end
    end
    msh_el = [msh_el; el_ids];
end
%% write mesh file
writeMeshfile('CVHM_mesh_modif.npsat', msh_nd, msh_el);
%}
%% ================ ELEVATION DATA based on CVHM grid =====================
% CVHM is divided into 10 layers therefore there are 11 2D grids. However
% the top elevation that descibes the surface is not used because the top
% elevation is set equal to the average hydraulic head.
%% expand the nodes by a mile
buff = shaperead('gis_data/CVHM_Mesh_outline_buffer');
%% or
load('BufferPnts.mat')
% make unique list of buffer outline points
buff_pnt = [buff(1,1).X(1) buff(1,1).Y(1)];
for ii = 1:length(buff(1,1).X) - 1
    dst = sqrt((buff_pnt(:,1) - buff(1,1).X(ii)).^2 + (buff_pnt(:,2) - buff(1,1).Y(ii)).^2);
    if min(dst) > 0.1
        buff_pnt = [buff_pnt; buff(1,1).X(ii) buff(1,1).Y(ii)];
    end
end
%% Read the Elevation from the DIS shapefile
dis = shaperead('/home/giorgk/Documents/UCDAVIS/cvhm_data/CVHM_database/DIS.shp');
clear BOT_ELEV
BOT_ELEV{10,1} = [];
for ii = 1:size(dis,1)
    ii
    xc = mean(dis(ii,1).X(1:4));
    yc = mean(dis(ii,1).Y(1:4));
    
    for jj = 1:10
        if jj == 10
            el = dis(ii,1).cvr2lay10b*0.3048;
        else
            el = dis(ii,1).(['cvr2lay' num2str(jj+1) 't'])*0.3048;
        end
        if el ~= 0
            BOT_ELEV{jj,1} = [BOT_ELEV{jj,1}; xc yc el];
        end
    end
end
%% Add the buffer points
for ii = 1:length(BOT_ELEV)
    Felev = scatteredInterpolant(BOT_ELEV{ii,1}(:,1), BOT_ELEV{ii,1}(:,2), BOT_ELEV{ii,1}(:,3));
    Felev.Method = 'nearest';
    Felev.ExtrapolationMethod = 'nearest';
    buff_pnt_Z = Felev(buff_pnt(:,1), buff_pnt(:,2));
    BOT_ELEV{ii,1} = [BOT_ELEV{ii,1}; buff_pnt(:,1), buff_pnt(:,2) buff_pnt_Z];
end
%% OLD dont use
ids = [2:10 10]';
for ii = 1:9 tp{ii,1} = 'Top';end; tp{10,1} = 'Bot';
for ii = 1:length(ids)
    [ BOT, info ] = readArcGisASCIIfile( ['gis_data/PP1766_DIS_Layer' num2str(ids(ii)) '_' tp{ii,1} '_ASCII.txt'] );
    BOT(BOT == BOT(1,1)) = nan;
    % find cells with elevation information
    [II, JJ] = find(BOT ~= 0 & ~isnan(BOT));
    if ii == 1
        % create the raster grid the first time an info file is read
        Xgrid = info.xllcorner+info.cellsize/2:...
                info.cellsize:...
                info.xllcorner+info.cellsize/2+info.cellsize*(info.ncols-1);
    
        Ygrid = info.yllcorner+info.cellsize/2:...
                info.cellsize:...
                info.yllcorner+info.cellsize/2+info.cellsize*(info.nrows-1);
    end
    
    % Interpolate from the BOT raster the buffer points
    Fbot = scatteredInterpolant(Xgrid(JJ)', Ygrid(II)', ...
        BOT(sub2ind(size(BOT),II,JJ)));
    Fbot.Method = 'nearest';
    Fbot.ExtrapolationMethod = 'nearest';
    buff_pnt_Z = Fbot(buff_pnt(:,1), buff_pnt(:,2));
    BOT_ELEV{ii,1} = [ Xgrid(JJ)', Ygrid(II)', ...
                       BOT(sub2ind(size(BOT),II,JJ));...
                       buff_pnt buff_pnt_Z ];
    
end
%% plot elevations
clf
hold on
for ii = 1:10
    plot3(BOT_ELEV{ii,1}(:,1), BOT_ELEV{ii,1}(:,2), BOT_ELEV{ii,1}(:,3),'.')
end

%% ============ BOTTOM Elevation file based on CVHM grid ==================
% write Bottom elevation file
BOT_mod = BOT_ELEV{10,1};
BOT_mod(BOT_mod(:,3)>0,3) = 0;
writeScatteredData('CVHM_Bot_elev.npsat', struct('PDIM',2,'TYPE','HOR','MODE','SIMPLE'), BOT_mod)
%% ============ Hydraulic Conductivity data based on CVHM grid ============
% read the data from C. Faunt
for ii = 1:10
    fid = fopen(['hyd_cond/HK_lyr' num2str(ii) '.txt'],'r');
    A = fscanf(fid,'%f', 441*98);
    HK{ii,1} = reshape(A,98,441)';
    fclose(fid);
    
    fid = fopen(['hyd_cond/VK_lyr' num2str(ii) '.txt'],'r');
    A = fscanf(fid,'%f', 441*98);
    VK{ii,1} = reshape(A,98,441)';
    fclose(fid);
end
%% plot HK
ii = 3;
if ii < 1; ii = 1;end
if ii > 10; ii = 10; end
M = HK{ii,1}(368:378,55:65);
surf(M,'edgecolor','none');
axis equal
axis ij
view(0,90);
%% tweak parts of the HK field if needed
for ii = 1:10
   %
   M = HK{ii,1}(320:332,17:22);
   M(M==0.053) = M(8,5);
   HK{ii,1}(320:332,17:22) = M;
   %
   M = HK{ii,1}(310:322,76:90);
   temp_val = min(M(M>1));
   if ~isempty(temp_val)
       M(M<1 & M > 0)=temp_val;
       HK{ii,1}(310:322,76:90) = M;
   end
end
%% 
% loop through the cells that contain HK and VK information and identify
% their coordinates in the rotated system.
HKALL = nan(50000,12);
VKALL = nan(50000,12);
cnt = 1;
for ii = 1:size(HK{1,1},1)
    for jj = 1:size(HK{1,1},2)
        ii
        if HK{1,1}(ii,jj) == 0
            continue;
        end
        % find the row in the BAS shapefile
        id = find([bas.ROW]' == ii & [bas.COLUMN_]' == jj);
        % calculate the center coordinate of the cell
        xc = mean(bas(id,1).X(1:end-1));
        yc = mean(bas(id,1).Y(1:end-1));
        hk = nan(1,10);
        vk = nan(1,10);
        for kk = 1:10
            hk(kk) = HK{kk,1}(ii,jj);
            vk(kk) = VK{kk,1}(ii,jj);
        end
      HKALL(cnt,:) = [xc yc hk];
      VKALL(cnt,:) = [xc yc vk];
      cnt = cnt + 1;
    end
end
HKALL(cnt:end,:)=[];
VKALL(cnt:end,:)=[];
%% 
% so far we have made a list of nodes with conductivity information.
% next we will extrapolate the values for the buffer nodes
for ii = 1:10
    ii
    % first remove those that are 0
    id_zero = find(HKALL(:,2+ii) == 0);
    if ~isempty(id_zero)
        id_nonzero = find(HKALL(:,2+ii) ~= 0);
        FHK = scatteredInterpolant(HKALL(id_nonzero,1), HKALL(id_nonzero,2), HKALL(id_nonzero,2+ii));
        FHK.Method = 'linear';
        FHK.ExtrapolationMethod = 'nearest';
        HKALL(id_zero,2+ii) = FHK(HKALL(id_zero,1), HKALL(id_zero,2));
    end
    
    FHK = scatteredInterpolant(HKALL(:,1), HKALL(:,2), HKALL(:,2+ii));
    FHK.Method = 'linear';
    FHK.ExtrapolationMethod = 'nearest';
    HK_buff(:,ii) = FHK(buff_pnt(:,1), buff_pnt(:,2));
    
    
    id_zero = find(VKALL(:,2+ii) == 0);
    if ~isempty(id_zero)
        id_nonzero = find(VKALL(:,2+ii) ~= 0);
        FVK = scatteredInterpolant(VKALL(id_nonzero,1), VKALL(id_nonzero,2), VKALL(id_nonzero,2+ii));
        FVK.Method = 'linear';
        FVK.ExtrapolationMethod = 'nearest';
        VKALL(id_zero,2+ii) = FVK(VKALL(id_zero,1), VKALL(id_zero,2));
    end
    
    FVK = scatteredInterpolant(VKALL(:,1), VKALL(:,2), VKALL(:,2+ii));
    FVK.Method = 'nearest';
    FVK.ExtrapolationMethod = 'nearest';
    VK_buff(:,ii) = FVK(buff_pnt(:,1), buff_pnt(:,2));
end
%%
% Next we will calculate the bottom elevation on the nodes where the
% hydraulic conductivity is defined. For the hydraulic conductivity
% functions we dont need the elevation of the last layer
XY_ALL = [HKALL(:,1) HKALL(:,2);buff_pnt(:,1), buff_pnt(:,2)];
HKbuff = [HKALL(:,3:end);HK_buff];
VKbuff = [VKALL(:,3:end);VK_buff];
for ii = 1:9
    FELEV = scatteredInterpolant(BOT_ELEV{ii,1}(:,1), BOT_ELEV{ii,1}(:,2),BOT_ELEV{ii,1}(:,3));
    FELEV.Method = 'linear';
    FELEV.ExtrapolationMethod = 'nearest';
    HVK_ELEV(:,ii) = FELEV(XY_ALL(:,1),XY_ALL(:,2));
end
%%
% Last, print to file
HKDATA = [XY_ALL HKbuff(:,1)];
VKDATA = [XY_ALL VKbuff(:,1)];
for ii = 1:9
    HKDATA = [HKDATA HVK_ELEV(:,ii) HKbuff(:,ii+1)];
    VKDATA = [VKDATA HVK_ELEV(:,ii) VKbuff(:,ii+1)];
end
%%
writeScatteredData('CVHM_HK.npsat', struct('PDIM',2,'TYPE','FULL','MODE','STRATIFIED'), HKDATA);
writeScatteredData('CVHM_VK.npsat', struct('PDIM',2,'TYPE','FULL','MODE','STRATIFIED'), VKDATA);

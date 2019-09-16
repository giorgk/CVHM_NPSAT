%% The stream data exists in the CVHM database.
% The original shapefile is the cvhm_sfr_network.
% In this repository we uploaded the same geometry split into the grid
% defined by the bas shapefile and deleting any fields not releated to
% streams.
%
cvhm_stream = shaperead('gis_data/CVHM_streams');
%% process each river independently.
% group the segments for each river
% initialize the custom structure with the first river segment
CVHMSTRM(1,1).Name = cvhm_stream(1,1).NAME;
CVHMSTRM(1,1).segments(1,1).X = cvhm_stream(1,1).X;
CVHMSTRM(1,1).segments(1,1).Y = cvhm_stream(1,1).Y;
CVHMSTRM(1,1).segments(1,1).CELLNUM = cvhm_stream(1,1).CELLNUM;
CVHMSTRM(1,1).segments(1,1).row = cvhm_stream(1,1).ROW;
CVHMSTRM(1,1).segments(1,1).col = cvhm_stream(1,1).COLUMN_;
for ii = 2:size(cvhm_stream,1)
    new_river = true;
    for jj = 1:size(CVHMSTRM,1)
        if strcmp(CVHMSTRM(jj,1).Name, cvhm_stream(ii,1).NAME)
            % add a segment to an existing river
            Nseg = size(CVHMSTRM(jj,1).segments,1);
            CVHMSTRM(jj,1).segments(Nseg+1,1).X = cvhm_stream(ii,1).X;
            CVHMSTRM(jj,1).segments(Nseg+1,1).Y = cvhm_stream(ii,1).Y;
            CVHMSTRM(jj,1).segments(Nseg+1,1).CELLNUM = cvhm_stream(ii,1).CELLNUM;
            CVHMSTRM(jj,1).segments(Nseg+1,1).row = cvhm_stream(ii,1).ROW;
            CVHMSTRM(jj,1).segments(Nseg+1,1).col = cvhm_stream(ii,1).COLUMN_;
            new_river = false;
            break;
        end
    end
    if new_river
        Nriv = size(CVHMSTRM,1);
        CVHMSTRM(Nriv+1,1).Name = cvhm_stream(ii,1).NAME;
        CVHMSTRM(Nriv+1,1).segments(1,1).X = cvhm_stream(ii,1).X;
        CVHMSTRM(Nriv+1,1).segments(1,1).Y = cvhm_stream(ii,1).Y;
        CVHMSTRM(Nriv+1,1).segments(1,1).CELLNUM = cvhm_stream(ii,1).CELLNUM;
        CVHMSTRM(Nriv+1,1).segments(1,1).row = cvhm_stream(ii,1).ROW;
        CVHMSTRM(Nriv+1,1).segments(1,1).col = cvhm_stream(ii,1).COLUMN_;
    end
end
%
%% for each river make a unique list of nodes
for ii = 1:size(CVHMSTRM,1)
    CVHMSTRM(ii,1).ND = [];
    for jj = 1:size(CVHMSTRM(ii,1).segments,1)
        [Xs, Ys] = polysplit(CVHMSTRM(ii,1).segments(jj,1).X, CVHMSTRM(ii,1).segments(jj,1).Y);
        for kk = 1:size(Xs,1)
            for mm = 1:length(Xs{kk,1})
                if isempty(CVHMSTRM(ii,1).ND)
                    CVHMSTRM(ii,1).ND = [Xs{kk,1}(mm) Ys{kk,1}(mm)];
                else
                    dst = sqrt((CVHMSTRM(ii,1).ND(:,1) - Xs{kk,1}(mm)).^2 + ...
                               (CVHMSTRM(ii,1).ND(:,2) - Ys{kk,1}(mm)).^2);
                    if min(dst) > 0.1
                        CVHMSTRM(ii,1).ND = [CVHMSTRM(ii,1).ND;...
                                             Xs{kk,1}(mm) Ys{kk,1}(mm)];
                    end
                end
            end
        end
    end
end
%% Create a Graph for each river
for ii = 1:size(CVHMSTRM,1)
    CVHMSTRM(ii,1).G = graph;
    for jj = 1:size(CVHMSTRM(ii,1).segments,1)
        [Xs, Ys] = polysplit(CVHMSTRM(ii,1).segments(jj,1).X, CVHMSTRM(ii,1).segments(jj,1).Y);
        for kk = 1:size(Xs,1)
            for mm = 1:length(Xs{kk,1})-1
                dst = sqrt((CVHMSTRM(ii,1).ND(:,1) - Xs{kk,1}(mm)).^2 + ...
                               (CVHMSTRM(ii,1).ND(:,2) - Ys{kk,1}(mm)).^2);
                id1 = find(dst < 0.1);
                dst = sqrt((CVHMSTRM(ii,1).ND(:,1) - Xs{kk,1}(mm+1)).^2 + ...
                               (CVHMSTRM(ii,1).ND(:,2) - Ys{kk,1}(mm+1)).^2);
                id2 = find(dst < 0.1);
                if id1 ~= id2
                    Newedge = table([id1 id2], jj, 'VariableNames',{'EndNodes','seg_id'});
                    CVHMSTRM(ii,1).G = addedge(CVHMSTRM(ii,1).G, Newedge);
                end
            end
        end
    end
end
%% Calculate the maximum numbender of connections for each river node
% for ii = 1:46
%     mx_conn(ii,1) = max(CVHMSTRM(ii,1).G.degree);
% end
%% plot a specific river
 ii = 8;
 p = CVHMSTRM(ii,1).G.plot;
 p.XData = CVHMSTRM(ii,1).ND(:,1);
 p.YData = CVHMSTRM(ii,1).ND(:,2);
 p.NodeLabel = [1:size(CVHMSTRM(ii,1).ND,1)];


%% Calculate Node normals for each river segment
% This section, besides node normals identifies the main rivers and the
% tributaries and makes separate paths for each.
for ii = 1:size(CVHMSTRM,1)
    clear PTHS NRM
    % find the terminal nodes
    id_term = find(CVHMSTRM(ii,1).G.degree == 1);
    if length(id_term) == 2 % there is only one path
        PTHS{1,1} = shortestpath(CVHMSTRM(ii,1).G, id_term(1), id_term(2));
    else
        nds = [];
        for m = 1:length(id_term)
            for n = m+1:length(id_term)
                if m ~= n
                    nds = [nds; id_term(m) id_term(n) length(shortestpath(CVHMSTRM(ii,1).G, id_term(m), id_term(n)))];
                end
            end
        end
        [cc, dd] = sort(nds(:,3),'descend');
        PTHS{1,1} = shortestpath(CVHMSTRM(ii,1).G, nds(dd(1),1), nds(dd(1),2));
        nd_started_from = [nds(dd(1),1); nds(dd(1),2)];
        for m = 2:length(dd)
            tmp1 = nds(dd(m),1);
            tmp2 = nds(dd(m),2);
            if isempty(find(nd_started_from == tmp1, 1))
                PTHS{m,1} = shortestpath(CVHMSTRM(ii,1).G, tmp1, tmp2);
                nd_started_from = [nd_started_from; tmp1];
            else 
                if isempty(find(nd_started_from == tmp2, 1))
                    PTHS{m,1} = shortestpath(CVHMSTRM(ii,1).G, tmp2, tmp1);
                    nd_started_from = [nd_started_from; tmp2];
                end
            end
        end
    end
    
    for jj = 1:size(PTHS,1)
        for k = 1:length(PTHS{jj,1})
           if k == 1
               a = [CVHMSTRM(ii,1).ND(PTHS{jj,1}(k),:) 0];
               b = [CVHMSTRM(ii,1).ND(PTHS{jj,1}(k+1),:) 0];
               dm = b; dm(3) = pdist([a;b]);
               aver = false;
           elseif k == length(PTHS{jj,1})
               a = [CVHMSTRM(ii,1).ND(PTHS{jj,1}(k-1),:) 0];
               b = [CVHMSTRM(ii,1).ND(PTHS{jj,1}(k),:) 0];
               dm = b; dm(3) = pdist([a;b]);
               aver = false;
           else
               a = [CVHMSTRM(ii,1).ND(PTHS{jj,1}(k-1),:) 0];
               b = [CVHMSTRM(ii,1).ND(PTHS{jj,1}(k),:) 0];
               c = [CVHMSTRM(ii,1).ND(PTHS{jj,1}(k+1),:) 0];
               aver = true;
           end
           if jj > 1 && CVHMSTRM(ii,1).G.degree(PTHS{jj,1}(k)) == 3
               aver = false;
           end
           if aver
               dm = b; dm(3) = pdist([a;b]);
               crs = cross(b-a,dm-a);
               N1 = crs/sqrt(sum(crs.^2));
               dm = c; dm(3) = pdist([b;c]);
               crs = cross(c-b,dm-b);
               N2 = crs/sqrt(sum(crs.^2));
               N = (N1+N2)/2;
           else
               crs = cross(b-a,dm-a);
               N = crs/sqrt(sum(crs.^2));
           end
           NRM{jj,1}(k,:) = N;
           
           if jj > 1 && CVHMSTRM(ii,1).G.degree(PTHS{jj,1}(k)) == 3
               PTHS{jj,1}(k+1:end) = [];
               break;
           end
        end
    end
    CVHMSTRM(ii,1).PATHS = PTHS;
    CVHMSTRM(ii,1).NORMS = NRM;
end
%
%% Make a unique list of stream cells and calculate the total river length 
% per cell
stream_cell_unique = [];
cell_riv_len = [];
for ii = 1:size(CVHMSTRM,1)
    for jj = 1:size(CVHMSTRM(ii,1).segments,1)
        r = CVHMSTRM(ii,1).segments(jj,1).row;
        c = CVHMSTRM(ii,1).segments(jj,1).col;
        if ~isempty(stream_cell_unique)
            id = find(stream_cell_unique(:,1) ==  r & stream_cell_unique(:,2) == c);
            if isempty(id)
                stream_cell_unique = [stream_cell_unique;r c];
                cell_riv_len = [cell_riv_len; 0];
                id = length(cell_riv_len);
            end
        else
            stream_cell_unique = [r c];
            cell_riv_len = 0;
            id = 1;
        end
        [Xs, Ys] = polysplit(CVHMSTRM(ii,1).segments(jj,1).X, CVHMSTRM(ii,1).segments(jj,1).Y);
        len = 0;
        for k = 1:size(Xs,1)
            for kk = 1:length(Xs{k,1})-1
                len = len + sqrt((Xs{k,1}(kk+1) - Xs{k,1}(kk))^2 + (Ys{k,1}(kk+1) - Ys{k,1}(kk))^2);
            end
        end
        CVHMSTRM(ii,1).segments(jj,1).riv_len = len;
        cell_riv_len(id) = cell_riv_len(id) + len;
    end
end
%% For each unique cell assign a stream flow volume
load('AvStresses.mat', 'STRMS')
CellFlow = zeros(size(stream_cell_unique,1),1);
for ii = 1:length(stream_cell_unique)
    id = find(STRMS(:,2) == stream_cell_unique(ii,1) & ...
              STRMS(:,3) == stream_cell_unique(ii,2));
    for jj = 1:length(id)
        CellFlow(ii) = CellFlow(ii) + STRMS(id(jj),4);
    end
   
end
%% Assign the rate and width to the cell according to the length
for ii = 1:size(CVHMSTRM,1)
    for jj = 1:size(CVHMSTRM(ii,1).segments,1)
        r = CVHMSTRM(ii,1).segments(jj,1).row;
        c = CVHMSTRM(ii,1).segments(jj,1).col;
        id = find(stream_cell_unique(:,1) ==  r & stream_cell_unique(:,2) == c);
        CVHMSTRM(ii,1).segments(jj,1).Qcell = CellFlow(id)*(CVHMSTRM(ii,1).segments(jj,1).riv_len/cell_riv_len(id));
        
        q_rate = 1000;
        width = 50;
        while abs(q_rate) > 0.4
            width = width + 50;
            q_rate = CVHMSTRM(ii,1).segments(jj,1).Qcell/(CVHMSTRM(ii,1).segments(jj,1).riv_len*width);
        end
        CVHMSTRM(ii,1).segments(jj,1).Width = width;
    end
end
%% loop through the segments and set the width and Q rate with the order 
% the graph nodes of the rivers are.
for ii = 1:size(CVHMSTRM,1)
    for mm = 1:size(CVHMSTRM(ii,1).PATHS,1)
        for k = 1:length(CVHMSTRM(ii,1).PATHS{mm,1})
            idm = 0;
            if k == 1
                id1 = CVHMSTRM(ii,1).PATHS{mm,1}(k);
                id2 = CVHMSTRM(ii,1).PATHS{mm,1}(k+1);
            elseif k == length(CVHMSTRM(ii,1).PATHS{mm,1})
                id1 = CVHMSTRM(ii,1).PATHS{mm,1}(k-1);
                id2 = CVHMSTRM(ii,1).PATHS{mm,1}(k);
            else
                id1 = CVHMSTRM(ii,1).PATHS{mm,1}(k-1);
                idm = CVHMSTRM(ii,1).PATHS{mm,1}(k);
                id2 = CVHMSTRM(ii,1).PATHS{mm,1}(k+1);
            end
            if idm == 0
                id_edge = CVHMSTRM(ii,1).G.findedge(id1,id2);
                jj = CVHMSTRM(ii,1).G.Edges.seg_id(id_edge);
                w = CVHMSTRM(ii,1).segments(jj,1).Width;
            else
                id_edge1 = CVHMSTRM(ii,1).G.findedge(id1,idm);
                id_edge2 = CVHMSTRM(ii,1).G.findedge(idm,id2);
                jj1 = CVHMSTRM(ii,1).G.Edges.seg_id(id_edge1);
                jj = CVHMSTRM(ii,1).G.Edges.seg_id(id_edge2);
                w = max(CVHMSTRM(ii,1).segments(jj1,1).Width, CVHMSTRM(ii,1).segments(jj,1).Width);
            end
            CVHMSTRM(ii,1).Path_nd_width{mm,1}(k) = w;
            if k ~= length(CVHMSTRM(ii,1).PATHS{mm,1})
                len = sqrt((CVHMSTRM(ii,1).ND(CVHMSTRM(ii,1).PATHS{mm,1}(k+1),1) - CVHMSTRM(ii,1).ND(CVHMSTRM(ii,1).PATHS{mm,1}(k),1) ).^2 + ...
                           (CVHMSTRM(ii,1).ND(CVHMSTRM(ii,1).PATHS{mm,1}(k+1),2) - CVHMSTRM(ii,1).ND(CVHMSTRM(ii,1).PATHS{mm,1}(k),2) ).^2);
                CVHMSTRM(ii,1).Path_Flow{mm,1}(k) = CVHMSTRM(ii,1).segments(jj,1).Qcell*(len / CVHMSTRM(ii,1).segments(jj,1).riv_len);
            end
            
        end
    end
end
%% check the total flow
TOT_check = 0;
for ii = 1:size(CVHMSTRM,1)
    for mm = 1:size(CVHMSTRM(ii,1).PATHS,1)
        TOT_check = TOT_check + sum(CVHMSTRM(ii,1).Path_Flow{mm,1});
    end
end
%% print the information for reading in Houdini
ii = 4;
cnt = 1;
for ii = 1:size(CVHMSTRM,1)
    for mm = 1:size(CVHMSTRM(ii,1).PATHS,1)
        fid = fopen(['stream_temp_files/temp_stream_per2_in_' num2str(cnt) '.txt'], 'w');
        fprintf(fid, '%d\n', length(CVHMSTRM(ii,1).PATHS{mm, 1}));
        fprintf(fid, '%f %f %f %f %f %f %f %f\n',...
            [CVHMSTRM(ii,1).ND(CVHMSTRM(ii,1).PATHS{mm, 1}',:)/100000 zeros(length(CVHMSTRM(ii,1).PATHS{mm, 1}),1)... % P {x,y,z}
             CVHMSTRM(ii,1).NORMS{mm, 1}(:,1:2) zeros(length(CVHMSTRM(ii,1).PATHS{mm, 1}),1) ... % N{x,y,z}
             CVHMSTRM(ii,1).Path_nd_width{mm,1}' [CVHMSTRM(ii,1).Path_Flow{mm,1} 0]'  ...
             ]');
        fclose(fid);
        cnt = cnt + 1;
    end
end
%
%% Read stream polygons from houdini.
% Run this part after the animation in houdini file "CVHM_HOU.hipnc" has
% been played. The animation creates one file per river path.
% Here we read the files and create the stream input file
cnt = 1;
clear STRM_POLY
STRM_POLY(1,1).Q = [];
STRM_POLY(1,1).Geo = [];
for ii = 1:53 % check the number of generated files
    fid = fopen(['stream_temp_files/temp_stream_per2_out' num2str(ii) '.txt'], 'r');
    Npoly = textscan(fid,'%d',1);
    for jj = 1:Npoly{1,1}
        C = textscan(fid,'%d %f',1);
        STRM_POLY(cnt,1).Q = C{1,2};
        poly = zeros(C{1,1},2);
        for k = 1:C{1,1}
            C = textscan(fid,'%f %f',1);
            poly(k,:) = [C{1,1} C{1,2}] ;
        end
        poly = poly*100000;
        D = pdist(poly);
        
        if min(D) < 0.5
            while 1
                D = squareform(D);
                break_in = false;
                for r = 1:size(D,1)
                    for c = 1:size(D,2)
                        if r == c
                            continue;
                        else
                            if D(r,c) == 0
                                poly(c,:) = [];
                                break_in = true;
                            end
                            if break_in
                                break;
                            end
                        end
                        if break_in
                            break;
                        end
                    end
                    if break_in
                        break;
                    end
                end
                D = pdist(poly);
                if min(D) > 0.5
                    break;
                end
            end
        end
        STRM_POLY(cnt,1).Geo = poly;
        cnt = cnt + 1;
    end
    fclose(fid);
end
STRM_POLY(find([STRM_POLY.Q]' == 0),:) = [];
%% remove the segments without area
dlt_id = [];
for ii = 1:size(STRM_POLY,1)
    Np = size(STRM_POLY(ii,1).Geo,1);
    if Np <=2
        dlt_id = [dlt_id; ii];
    end
end
STRM_POLY(dlt_id,:) = [];
%% Write the stream file (finally!)
fid = fopen('CVHM_streams_per2.npsat','w');
fprintf(fid, '%d\n', size(STRM_POLY,1));
for ii = 1:size(STRM_POLY,1)
    Np = size(STRM_POLY(ii,1).Geo,1);
    if Np <=2
        continue;
    end
    
    if Np == 3
        A = polyarea(STRM_POLY(ii,1).Geo(:, 1), STRM_POLY(ii,1).Geo(:, 2));
        fprintf(fid, '%d %f\n', [size(STRM_POLY(ii,1).Geo,1) STRM_POLY(ii,1).Q/A]);
        fprintf(fid, '%f %f\n', STRM_POLY(ii,1).Geo');
    elseif Np == 4
        A = polyarea(STRM_POLY(ii,1).Geo([1 2 4 3],1), STRM_POLY(ii,1).Geo([1 2 4 3],2));
        fprintf(fid, '%d %f\n', [size(STRM_POLY(ii,1).Geo,1) STRM_POLY(ii,1).Q/A]);
        fprintf(fid, '%f %f\n', STRM_POLY(ii,1).Geo([1 2 4 3],:)');
    else
        chech_that = true;
    end
end
fclose(fid);


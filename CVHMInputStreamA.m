function CVHMInputStreamA(STRMS, startTime, endTime, opt)

load('CVHMSTRM_Graph.mat', 'CVHMSTRM');

timestring = [num2str(startTime(1)) 'm' num2str(startTime(2)) '_'  num2str(endTime(1)) 'm' num2str(endTime(2))];

% Make a unique list of stream cells and calculate the total river length 
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

% For each unique cell assign a stream flow volume
CellFlow = zeros(size(stream_cell_unique,1),1);
for ii = 1:length(stream_cell_unique)
    CellFlow(ii) = CellFlow(ii) + STRMS(stream_cell_unique(ii,1), stream_cell_unique(ii,2));
end

% Assign the rate and width to the cell according to the length
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

% loop through the segments and set the width and Q rate with the order 
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

if ~exist([opt.simFolder filesep 'stream_temp_files'], 'dir')
    mkdir([opt.simFolder filesep 'stream_temp_files']);
end

% print the information for reading in Houdini
cnt = 1;
for ii = 1:size(CVHMSTRM,1)
    for mm = 1:size(CVHMSTRM(ii,1).PATHS,1)
        fid = fopen([opt.simFolder filesep 'stream_temp_files/temp_' timestring '_in_' num2str(cnt) '.txt'], 'w');
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

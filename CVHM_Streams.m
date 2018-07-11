%% The stream data exists in the CVHM database.
% The original shapefile is the cvhm_sfr_network.
% In this repository we uploaded the same geometry split into the grid
% defined by the bas shapefile and deleting any fiels not releated to
% streams.
%{
cvhm_stream = shaperead('gis_data/CVHM_streams');
%% process each river independently.
% group the segments for each river
% initialize the custom structure with the first river segments
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
                    CVHMSTRM(ii,1).G = addedge(CVHMSTRM(ii,1).G, id1, id2);
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
% ii = 9;
% p = CVHMSTRM(ii,1).G.plot;
% p.XData = CVHMSTRM(ii,1).ND(:,1);
% p.YData = CVHMSTRM(ii,1).ND(:,2);
% p.NodeLabel = [1:size(CVHMSTRM(ii,1).ND,1)];


%% Calculate Node normals for each river segment
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
               break;
           end
        end
    end
    CVHMSTRM(ii,1).PATHS = PTHS;
    CVHMSTRM(ii,1).NORMS = NRM;
end
%}
%% Make a unique list of stream cells and calculate the total river length 
% per cell
stream_cell_unique = [];
for ii = 1:size(CVHMSTRM,1)
    for jj = 1:size(CVHMSTRM(ii,1).segments,1)
        
    end
end
%% Add rate and width information
for ii = 1:size(CVHMSTRM,1)
    for jj = 1:size(CVHMSTRM(ii,1).segments,1)
        r = CVHMSTRM(ii,1).segments(jj,1).row;
        c = CVHMSTRM(ii,1).segments(jj,1).col;
        CVHMSTRM(ii,1).segments(jj,1).Qm3day = 0;
        id = find(STRMS(:,2) == r & STRMS(:,3) == c);
        for k = 1:length(id)
            CVHMSTRM(ii,1).segments(jj,1).Qm3day = CVHMSTRM(ii,1).segments(jj,1).Qm3day + STRMS(id(k),4);
        end
    end
end
%% print the information for reading in Houdini
fid = fopen('temp.txt', 'w');
fprintf(fid, '%f %f %f %f %f %f\n',...
    [CVHMSTRM(1,1).ND(CVHMSTRM(1).PATHS{1, 1}',:)/100000 zeros(147,1)...
     CVHMSTRM(1).NORMS{1, 1}(:,1:2) zeros(147,1)]');
fclose(fid);


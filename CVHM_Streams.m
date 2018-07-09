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
    ii
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
%}
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
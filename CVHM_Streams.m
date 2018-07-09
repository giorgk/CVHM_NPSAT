%% The stream data exists in the CVHM database.
% The original shapefile is the cvhm_sfr_network.
% In this repository we uploaded the same geometry split into the grid
% defined by the bas shapefile and deleting any fiels not releated to
% streams.
cvhm_stream = shaperead('gis_data/CVHM_streams');
%% process each river independently.
% first make a unique list of rivers based on their name
unique_names = [];
% initialize the custom structure with the first river segments
CVHMSTRM(1,1).Name = cvhm_stream(1,1).NAME;
CVHMSTRM(1,1).segments(1,1).X = cvhm_stream(1,1).X;
CVHMSTRM(1,1).segments(1,1).Y = cvhm_stream(1,1).Y;
CVHMSTRM(1,1).segments(1,1).CELLNUM = cvhm_stream(1,1).CELLNUM;
CVHMSTRM(1,1).segments(1,1).row = cvhm_stream(1,1).ROW;
CVHMSTRM(1,1).segments(1,1).col = cvhm_stream(1,1).COLUMN_;
for ii = 2:size(cvhm_stream,1)
    for jj = 1:size(CVHMSTRM,1)
        if strcmp(CVHMSTRM(jj,1).Name, cvhm_stream(ii,1).NAME)
            % add a segment to an existing river
            Nseg = size(CVHMSTRM(jj,1).segments,1);
            CVHMSTRM(jj,1).segments(Nseg+1,1).X = cvhm_stream(ii,1).X;
            CVHMSTRM(jj,1).segments(Nseg+1,1).Y = cvhm_stream(ii,1).Y;
            CVHMSTRM(jj,1).segments(Nseg+1,1).CELLNUM = cvhm_stream(ii,1).CELLNUM;
            CVHMSTRM(jj,1).segments(Nseg+1,1).row = cvhm_stream(ii,1).ROW;
            CVHMSTRM(jj,1).segments(Nseg+1,1).col = cvhm_stream(ii,1).COLUMN_;
        end
    end
end
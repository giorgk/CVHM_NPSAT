function Wells_Gen = writeWells(Wells_Gen, botFile,topFile, opt)
% This function calculates and adjust the absolute top and bottom well
% screen elevations of the wells and writes them into file
%
% Wells_Gen are the generated wells
% botFile and topFile are the files that define the bottom and top elevation.
% It is assumed that they are written as scattered interpolation file

if ~isfield(opt, 'above_thres')
    opt.above_thres = 10;
end

if ~isfield(opt, 'below_thres')
    opt.below_thres = 30;
end

if ~isfield(opt, 'minScreen')
    opt.minScreen = 10;
end

if ~isfield(opt, 'wellVersion')
    opt.wellVersion = 1;
end


topelev = read_Scattered(topFile, 2);
Ftop = scatteredInterpolant(topelev.p(:,1), topelev.p(:,2), topelev.v, 'natural');

% abd the bottom of the aquifer
botelev = read_Scattered(botFile, 2);
Fbot = scatteredInterpolant(botelev.p(:,1), botelev.p(:,2), botelev.v, 'natural');

WW = [[Wells_Gen.X]' [Wells_Gen.Y]' [Wells_Gen.D]' [Wells_Gen.SL]' [Wells_Gen.Q]'];

cnvrs = 0.3048;
for ii = 1:size(WW,1)
    Wt = Ftop(WW(ii,1), WW(ii,2)); % Water table
    Bs = Fbot(WW(ii,1), WW(ii,2)); %Aquifer Base
    aquif_height = Wt - Bs;
    depth = WW(ii,3)*cnvrs;
    screen = WW(ii,4)*cnvrs;
    if aquif_height < 0
        warning(['The aquifer height for the ' num2str(ii) 'well is negative'])
    elseif aquif_height < opt.above_thres + opt.below_thres + opt.minScreen
        % This is an extreme case where the aquifer is too shallow
        % split the aquifer height in 3/3
        topw = Wt - aquif_height/3;
        botw = topw - aquif_height/3;
    elseif aquif_height < max(depth,opt.above_thres) + screen + opt.below_thres
        % if the combination of screen depth doesnt fit in the aquifer we
        % will sacrifice the larger
        if depth > screen
            botw = Bs + opt.below_thres;
            topw = botw + screen;
            if topw > Wt - opt.above_thres
                topw = Wt - opt.above_thres;
            end
        else
            topw = Wt - max(depth,opt.above_thres);
            botw = Bs + opt.below_thres;
            if topw - botw < opt.minScreen
                topw = Wt - aquif_height/3;
                botw = topw - aquif_height/3;
            end
        end
    else
        topw = Wt - max(depth,opt.above_thres);
        botw = topw - screen;
    end
    
    if any([Wt-topw<0 topw-botw<0 botw-Bs<0])
        topw = Wt - aquif_height/3;
        botw = topw - aquif_height/3;
    end
    
    WW(ii,8) = topw;
    WW(ii,9) = botw;
    WW(ii,10) = Wt;
    WW(ii,11) = Bs;
    Wells_Gen(ii,1).Top = topw;
    Wells_Gen(ii,1).Bot = botw;
    Wells_Gen(ii,1).WT = Wt;
    Wells_Gen(ii,1).Base = Bs;
end

% Print wells to file
fid = fopen([opt.simFolder filesep opt.prefix '_' opt.timestring  '_Wells_' num2str(opt.wellVersion) '.npsat'],'w');
fprintf(fid, '%d\n', size(WW, 1));
fprintf(fid, '%f %f %f %f -%f\n', [WW(:,1) WW(:,2) WW(:,8) WW(:,9) WW(:,5)]');
fclose(fid);
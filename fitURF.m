function fittedParam = fitURF( prefix, nproc_sim)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


 
% fit options


m(nproc_sim,1).data = []; 
s(nproc_sim,1).data = [];
wgh(nproc_sim,1).data = [];
parfor ii = 0:nproc_sim-1
    ft = fittype( '(1/(x*b*sqrt(2*pi)))*exp((-(log(x)-a)^2)/(2*b^2))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.499112701571756 0.336174051321482];
    X = 1:200;
    
    m_temp = nan(20000,100);
    s_temp = nan(20000,100);
    wgh_temp = nan(20000,100);
    w = load([prefix num2str(ii,'%04d') '.mat']);
    
    for k = 1:length(w.wellurfs)
        disp([ii k]);
        for j = 1:length(w.wellurfs(k,1).DATA)
            tmp = w.wellurfs(k,1).DATA(j,1);
            [xData, yData] = prepareCurveData( X, tmp.URF );
            [fitresult, gof] = fit( xData, yData, ft, opts );
            ab = coeffvalues(fitresult);
            m_temp(tmp.Eid+1,tmp.Sid+1) = ab(1);
            s_temp(tmp.Eid+1,tmp.Sid+1) = ab(2);
            wgh_temp(tmp.Eid+1,tmp.Sid+1) = tmp.v_cds;
        end
    end
    m(ii+1).data = m_temp;
    s(ii+1).data = s_temp;
    wgh(ii+1).data = wgh_temp;
end

m_all = nan(20000,100);
s_all = nan(20000,100);
wght_all = nan(20000,100);

for ii = 1:length(m)
    ind = find(~isnan(m(ii,1).data));
    m_all(ind) = m(ii,1).data(ind);
    s_all(ind) = s(ii,1).data(ind);
    wght_all(ind) = wgh(ii,1).data(ind);
end


fittedParam.mu = m_all;
fittedParam.std = s_all;
fittedParam.wgh = wght_all;
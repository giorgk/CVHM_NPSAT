function PDF = Calc_2D_PDF(X, Y)
%% generated 2D propability density functions

in = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(in);
Y = Y(in);

xlm = [min(X) max(X)];
ylm = [min(Y) max(Y)];
dx = linspace(xlm(1), xlm(2), 100);
dy = linspace(ylm(1), ylm(2), 100);
dx_n = linspace(0, 100, 100);
dy_n = linspace(0,100, 100);
% rescale
X_n = interp1(dx, dx_n, X);
Y_n = interp1(dy, dy_n, Y);
%
%%
w = 20;
V = zeros(length(dx_n),length(dy_n));
for ii = 1:length(dx_n)
    for jj = 1:length(dy_n)
        % find the number of points inside the kernel
        dst = sqrt((dx_n(ii) - X_n).^2 + (dy_n(jj) - Y_n).^2);
        V(jj,ii) = sum(dst < w);
    end
end
V = V./max(max(V));
PDF.X = dx;
PDF.Y = dy;
PDF.V = V;
[nx, ny] = ndgrid(PDF.X, PDF.Y);
PDF.F = griddedInterpolant(nx, ny, PDF.V','linear','linear');



% clf
% surf(dx,dy,PDF,'edgecolor','none')
% hold on 
% plot(X,Y,'.')
% alpha(0.5)
% view(0,90)
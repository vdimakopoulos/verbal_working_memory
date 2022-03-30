function createfigure(X1, YMatrix1, X2, Y1, X3, Y2)
%CREATEFIGURE(X1, YMatrix1, X2, Y1, X3, Y2)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  X2:  vector of x data
%  Y1:  vector of y data
%  X3:  vector of x data
%  Y2:  vector of y data

%  Auto-generated by MATLAB on 19-Mar-2020 18:32:02

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to semilogx
semilogx1 = semilogx(X1,YMatrix1,'LineWidth',2,'Parent',axes1);
set(semilogx1(1),'Color',[0 0 0]);
set(semilogx1(2),'Color',[1 0 0]);
set(semilogx1(3),'Color',[0 1 0]);

% Create semilogx
semilogx(X2,Y1,'LineWidth',3,'Color',[0 0 0]);

% Create semilogx
semilogx(X3,Y2,'LineWidth',3,'Color',[1 0 0]);

% Create ylabel
ylabel({''});

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[3 100]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 1]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',30,'XMinorTick','on','XScale','log','XTickLabel',...
    {'10','100'});
% Create textbox
annotation(figure1,'textbox',...
    [0.323214285714286 0.802380952380954 0.628571428571428 0.109523809523811],...
    'String',{'ECoG C3 - AHL2-3'},...
    'FontSize',26,...
    'FitBoxToText','off',...
    'EdgeColor','none');

function createfigure(ymatrix1)
%CREATEFIGURE(ymatrix1)
%  YMATRIX1:  bar matrix data

%  Auto-generated by MATLAB on 12-Dec-2019 13:17:22

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to bar
bar1 = bar(ymatrix1,'Horizontal','on');
set(bar1(1),'FaceColor',[1 0 0]);

% Create xlabel
xlabel('Percentage Changes on output due to -% red and +% blue change in input');

% Create title
title('Sensitivities');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'YTick',[1 2 3 4 5],'YTickLabel',...
    {'  Amount of Electricity, at grid in Electricity, at grid','  Amount of Electricity, at grid in Transport, pipeline, unspecified petroleum products','  Amount of Transport, pipeline, unspecified petroleum products in Transport, school bus, diesel powered','  Amount of Transport, pipeline, unspecified petroleum products in Transport, pipeline, unspecified petroleum products','  Amount of Electricity, bituminous coal, at power plant in Electricity, biomass, at power plant'});
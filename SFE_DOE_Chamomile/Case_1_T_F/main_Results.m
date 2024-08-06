startup;
delete(gcp('nocreate'));
%p = Pushbullet(pushbullet_api);
%initParPool

%addpatch('casadi_folder')
%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
%addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
%import casadi.*

%% Load data
Parameters_table        = readtable('Parameters.csv') ;                     % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});                  % Parameters within the model + (m_max), m_ratio, sigma
r                       = Parameters{3};                                    % Radius of the extractor  [m]
epsi                    = Parameters{4};                                    % Fullness [-]
dp                      = Parameters{5};                                    % Paritcle diameter
L                       = Parameters{6};                                    % Total length of the extractor [m]

V                       = L  * pi * r^2;                                    % Total volume of the extractor [m3]
A                       = pi *      r^2;                                    % Extractor cross-section

%--------------------------------------------------------------------

N_exp                   = 50;
COLORS                  = ['b','r','k','m','g'];
COST_I   = []; COST_F   = []; GROUP    = []; Yield      = [];
PP                      = [100, 125, 150, 175, 200];
%hold on
for ii = 1:numel(PP)
    PRES = PP(ii);
    AA       = readlines(['Cost_',num2str(PRES),'.txt']);
    InitCost = str2num(AA(1));
    FinaCost = str2num(AA(2));
    COST_I   = [ COST_I , InitCost];
    COST_F   = [ COST_F , FinaCost];
    GROUP    = [GROUP; repmat(['P = ', num2str(PRES),' bar'], numel(FinaCost),1) ];
    
    AA       = readlines(['Control_',num2str(PRES),'.txt']);
    TempCont = str2num(AA(1));
    TempCont = reshape(TempCont,[],N_exp);
    FlowCont = str2num(AA(2));
    FlowCont = reshape(FlowCont,[],N_exp);
    
    Time     = linspace(0,300,size(TempCont,1)+1);

    % Plots of controls
    ind      = find( FinaCost == min(FinaCost));
    
    %{\
    subplot(3,1,1)
    hold on
    stairs(Time, [TempCont(:, ind); TempCont(end,ind)]-273, 'LineWidth', 2 ,'Color',COLORS(ii));
    hold off
    ylabel('$T^{in}~^\circ$C')
    xlabel('Time [min]')
    %legend('P = 100 bar','P = 125 bar','P = 150 bar','P = 175 bar','P = 200 bar', 'Location','best','NumColumns',5)
    %legend box off
    %axis square
    set(gca,'FontSize',16)
    
    subplot(3,1,2)
    hold on
    stairs(Time, [FlowCont(:, ind); FlowCont(end,ind)], 'LineWidth', 2, 'Color',COLORS(ii));
    hold off
    ylabel('$F\cdot 10^{-5}$ $kg/s$')
    xlabel('Time [min]')
    %axis square
    set(gca,'FontSize',16)
    %}

    Yield = [Yield; Yield_Plot(PRES,TempCont(:, ind),FlowCont(:, ind))];

    Time     = linspace(0,300,size(Yield,2));

    subplot(3,1,3)
    hold on
    plot(Time, Yield(ii,:), 'LineWidth',3, 'Color', COLORS(ii));
    hold off
    xlabel('Time min')
    ylabel('y gram')
    set(gca,'FontSize',16);
end
%exportgraphics(figure(1), ['Profiles_all.png'], "Resolution",300); close all

%%
%{
figure()
s = scatterhist(COST_F, COST_I, 'Group', GROUP, 'Kernel', 'on', 'LineWidth',3, 'MarkerSize',6, 'Color',COLORS, 'LineStyle',{'-','-','-','-','-','-'} );
s(1).Children(5).MarkerFaceColor = 'b';
s(1).Children(4).MarkerFaceColor = 'r';
s(1).Children(3).MarkerFaceColor = 'k';
s(1).Children(2).MarkerFaceColor = 'm';
s(1).Children(1).MarkerFaceColor = 'g';

legend box off

ylabel('Inital cost of the objective function [-]');
xlabel('Final cost of the objective function [-]');
set(gca,'FontSize',10);

%exportgraphics(figure(1), ['scatter.png'], "Resolution",300);
%close all

%}
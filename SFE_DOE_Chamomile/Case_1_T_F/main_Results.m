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

N_exp    = 50;
COST_I   = []; COST_F   = []; GROUP    = [];
%hold on
for ii = [100, 125, 150, 175, 200]
    AA       = readlines(['Cost_',num2str(ii),'.txt']);
    InitCost = str2num(AA(1));
    FinaCost = str2num(AA(2));
    COST_I   = [ COST_I , InitCost];
    COST_F   = [ COST_F , FinaCost];
    GROUP    = [GROUP; repmat(['P = ', num2str(ii),' bar'], numel(FinaCost),1) ];
    
    AA       = readlines(['Control_',num2str(ii),'.txt']);
    TempCont = str2num(AA(1));
    TempCont = reshape(TempCont,[],N_exp);
    FlowCont = str2num(AA(2));
    FlowCont = reshape(FlowCont,[],N_exp);
    
    Time     = linspace(0,300,size(TempCont,1)+1);

    % Plots of controls
    ind      = find( FinaCost == min(FinaCost));
    
    figure(1)
    hold on
    stairs(Time, [TempCont(:, ind); TempCont(end,ind)]-273, 'LineWidth', 2);
    hold off
    
    figure(2)
    hold on
    stairs(Time, [FlowCont(:, ind); FlowCont(end,ind)], 'LineWidth', 2);
    hold off

    %exportgraphics(figure(1), ['Profile_',num2str(ii),'.png'], "Resolution",300);
    %close all   
end

figure(1)
ylabel('Temperature [$^\circ$C]')
xlabel('Time [min]')
set(gca,'FontSize',12)
legend('P = 100 bar','P = 125 bar','P = 150 bar','P = 175 bar','P = 200 bar', 'Location','best')
legend box off
exportgraphics(figure(1), ['Profile_T.png'], "Resolution",300);
close figure 1

figure(2)
ylabel('Mass flow rate $\times 10^{-5}$ [$kg/s$]')
xlabel('Time [min]')
set(gca,'FontSize',12)
legend('P = 100 bar','P = 125 bar','P = 150 bar','P = 175 bar','P = 200 bar', 'Location','bestoutside')
legend box off
%exportgraphics(figure(2), ['Profile_F.png'], "Resolution",300);
%close figure 2

%%
figure(1)
s = scatterhist(COST_F, COST_I, 'Group', GROUP, 'Kernel', 'on', 'LineWidth',3, 'MarkerSize',6, 'Color',['b','r','k','m','g'], 'LineStyle',{'-','-','-','-','-','-'} );
s(1).Children(5).MarkerFaceColor = 'b';
s(1).Children(4).MarkerFaceColor = 'r';
s(1).Children(3).MarkerFaceColor = 'k';
s(1).Children(2).MarkerFaceColor = 'm';
s(1).Children(1).MarkerFaceColor = 'g';

legend box off

ylabel('Inital cost of the objective function [-]');
xlabel('Final cost of the objective function [-]');
set(gca,'FontSize',12);

exportgraphics(figure(1), ['scatter.png'], "Resolution",300);
close all

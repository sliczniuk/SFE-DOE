startup;
delete(gcp('nocreate'));
% %p = Pushbullet(pushbullet_api);

%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

%%
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

N_exp                   = 32;
COLORS                  = ['b','r','k','m','g'];
COST_I   = []; COST_F   = []; GROUP    = []; Yield      = [];
PP                      = [100, 125, 150, 175, 200];
%PP                      = [150, 175];

figure(3)
tiledlayout(numel(PP),2)

%hold on
for ii = 1:numel(PP)
    PRES = PP(ii);
    AA       = readlines(['COST_',num2str(PRES),'.txt']);
    
    InitCost = str2num(AA(1,:));
    FinaCost = str2num(AA(2,:));
    COST_I   = [ COST_I , InitCost];
    COST_F   = [ COST_F , FinaCost];
    GROUP    = [GROUP; repmat(['P = ', num2str(PRES),' bar'], numel(FinaCost),1) ];
    
    AA       = readlines(['CONTROL_',num2str(PRES),'.txt']);

    TempCont = str2num(AA(1, :));
    TempCont = reshape(TempCont,[],N_exp);
    FlowCont = str2num(AA(2, :)) * 1e-5;
    FlowCont = reshape(FlowCont,[],N_exp);

    Y_FP     = readlines(['FP_',num2str(PRES),'.txt']);
    YY_FP_0  = reshape(str2num(Y_FP(1,:)),[],N_exp);
    YY_FP    = reshape(str2num(Y_FP(2,:)),[],N_exp); YY_FP(end,:)

    Time     = linspace(0,600,size(TempCont,1)+1);

    % Plots of controls
    selected_exp  = sort(FinaCost);
    ind      = find( FinaCost == selected_exp(1) );
    
    %{\
    %subplot(3,1,1)
    figure(1)
    hold on
    stairs(Time, [TempCont(:, ind); TempCont(end,ind)]-273, 'LineWidth', 2 ,'Color',COLORS(ii));
    hold off
    ylabel('$T^{in}~^\circ$C')
    xlabel('Time min')
    %legend('P = 100 bar','P = 125 bar','P = 150 bar','P = 175 bar','P = 200 bar', 'Location','best','NumColumns',5)
    %legend box off
    %axis square
    set(gca,'FontSize',16)
    
    figure(2)
    %subplot(3,1,2)
    hold on
    stairs(Time, [FlowCont(:, ind); FlowCont(end,ind)], 'LineWidth', 2, 'Color',COLORS(ii));
    hold off
    ylabel('$F~kg/s$')
    xlabel('Time min')
    %axis square
    set(gca,'FontSize',16)
    %}
    %{'
    %Yield = [Yield; Yield_Plot(PRES,TempCont(:, ind),FlowCont(:, ind))];

    Time     = linspace(0,600,size(YY_FP,1));

    figure(3)
    %subplot(5,1,ii)
    %nexttile
    hold on
    %plot(Time, YY_RBF(:,ii), 'LineWidth',3, 'Color', 'r', 'LineStyle',':');
    plot(Time, YY_FP( :,ind), 'LineWidth',3, 'Color', COLORS(ii), 'LineStyle','-');
    %plot(Time, YY_RBF_0(:,ii), 'LineWidth',3, 'Color', 'b', 'LineStyle',':');
    %plot(Time, YY_FP_0( :,ii), 'LineWidth',3, 'Color', 'b', 'LineStyle','-');
    hold off
    xlabel('Time min')
    ylabel('y gram')
    %ylim([0 2.8])
    %title(['P =',num2str(round(PRES)),' bar'])
    set(gca,'FontSize',16);

    figure(4)
    %subplot(5,1,ii)
    %nexttile
    hold on
    plot(Time(1:end-1), diff(YY_FP(:,ind)), 'LineWidth',3, 'Color', COLORS(ii),'LineStyle','-');
    %plot(Time(1:end-1), diff(YY_RBF(:,ii)), 'LineWidth',3, 'Color', 'r','LineStyle',':');
    %plot(Time(1:end-1), diff(YY_FP(:,ii)), 'LineWidth',3, 'Color', 'r','LineStyle','-');
    %plot(Time(1:end-1), diff(YY_RBF_0(:,ii)), 'LineWidth',3, 'Color', 'b','LineStyle',':');
    %plot(Time(1:end-1), diff(YY_FP_0(:,ii)), 'LineWidth',3, 'Color', 'b','LineStyle','-');
    hold off
    xlabel('Time min')
    ylabel('$\frac{dy}{dt}$ gram/s')
    %ylim([0 0.13])
    %title(['P =',num2str(round(PRES)),' bar'])
    set(gca,'FontSize',16);
%}

end
%{\
figure(1); legend({'100 bar','125 bar','150 bar','175 bar','200 bar'}, 'Location', 'southwest' ); legend('boxoff')
exportgraphics(figure(1), ['Profile_T.png'], "Resolution",300); 

figure(2); legend({'100 bar','125 bar','150 bar','175 bar','200 bar'}, 'Location', 'southeast' ); legend('boxoff')
exportgraphics(figure(2), ['Profile_F.png'], "Resolution",300); 

figure(3); legend({'100 bar','125 bar','150 bar','175 bar','200 bar'}, 'Location', 'northwest' ); legend('boxoff')
exportgraphics(figure(3), ['yield.png'], "Resolution",300); 

figure(4); legend({'100 bar','125 bar','150 bar','175 bar','200 bar'}, 'Location', 'northeast' ); legend('boxoff')
exportgraphics(figure(4), ['diff_yield.png'], "Resolution",300); 
%}
%%
%{\
figure(5)
s = scatterhist(COST_F, COST_I, 'Group', GROUP, 'Kernel', 'on', 'LineWidth',3, 'MarkerSize',6, 'Color',COLORS, 'LineStyle',{'-','-','-','-','-','-'} );
%s = scatter(COST_F, COST_I, 'Group', GROUP, 'LineWidth',3, 'MarkerSize',6, 'Color',COLORS, 'LineStyle',{'-','-','-','-','-','-'} );
s(1).Children(5).MarkerFaceColor = 'b';
s(1).Children(4).MarkerFaceColor = 'r';
s(1).Children(3).MarkerFaceColor = 'k';
s(1).Children(2).MarkerFaceColor = 'm';
s(1).Children(1).MarkerFaceColor = 'g';

legend box off

ylabel('Inital value of $-\ln j_D$ [-]');
xlabel('Final value of $-\ln j_D$ [-]');
set(gca,'FontSize',16);

exportgraphics(figure(5), ['scatter.png'], "Resolution",300);
close all

%}
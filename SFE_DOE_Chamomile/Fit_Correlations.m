startup;
delete(gcp('nocreate'));

%%
AA = xlsread('Regression.xlsx');
DI = [0.71383013	1.076801997	2.179470155	2.475532632	1.390707877	1.336111172	1.882954204	2.457886055	0.564935512	1.542106938	0.835725102	0.87349666];
GG = [4.229739602	3.091520556	2.359538225	1.132795818	2.204975712	2.739220425	1.868538631	1.69935869	3.452308202	1.995905641	3.012676539	2.596460037];
RE = [0.4632, 0.3783, 0.3029, 0.2619, 0.3579, 0.3140, 0.2635, 0.2323, 0.1787, 0.1160, 0.1889, 0.1512];
TT = [313, 313, 313, 313, 303, 303, 303, 303, 303, 303, 313, 313];
Tr = TT ./ 304;
RHO= [630, 691, 793, 840, 772, 802, 856, 891, 772, 891, 691, 793];

delta = 2.79 .* sqrt(7.38) .* Tr.^(1/4) .* (RHO./470);

FF = [6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 3.33, 3.33, 3.33, 3.33];

%%
x1 = RE(1:8) ; y1 = DI(1:8);
x2 = RE(9:12); y2 = DI(9:12);

[mdl1,stuff1] = fit(x1',y1','poly1');
[mdl2,stuff2] = fit(x2',y2','poly1');

p1 = coeffvalues(mdl1);
p2 = coeffvalues(mdl2);

xtest  = linspace(0.1,0.5);
ytest1 = polyval(p1,xtest);
ytest2 = polyval(p2,xtest);

hold on
plt1 = plot(x1,y1,'ko', 'LineWidth',2);
plt2 = plot(x2,y2,'kd', 'LineWidth',2);
plot(xtest,ytest1, 'LineWidth',2);
plot(xtest,ytest2, 'LineWidth',2);
yline(0, 'k', 'LineWidth',2);
hold off

%text(0.3, 4.0, sprintf('~$~\\Upsilon = %.3f \\cdot + Re %.3f$ \n $R^2=%.3f$',[p1, stuff1.rsquare]))
%text(0.15, 6.0, sprintf('~$~ \\Upsilon = %.3f \\cdot Re + %.3f$ \n $R^2=%.3f$',[p2, stuff2.rsquare]))

text(0.35, 2.0, sprintf('$~ D_i^R = %.3f \\cdot Re + %.3f$ \n $R^2=%.3f$',[p1, stuff1.rsquare]))
text(0.15, -1.0, sprintf('$~ D_i^R = %.3f \\cdot Re + %.3f$ \n $R^2=%.3f$',[p2, stuff2.rsquare]))

legend([plt1(1),plt2(1)],'$6.67\times 10^{-5}[kg/s]$','$3.33\times 10^{-5}[kg/s]$','Location','best');
legend box off

ylabel('$D_i^R \times 10^{-13} [m^2/s]$')
%ylabel('$\Upsilon [-]$')
xlabel('Re [-]')
set(gca,'FontSize',12)

%exportgraphics(figure(1), ['Correlation_Di_Re.png'], "Resolution",300);
%close all

%% gamma function plot
%{\
C0 = 1;
XX = C0:-0.01:0;

figure()
for ii=1:12

    %RE_temp = sort(RE, 'descend');
    %jj = find(RE_temp(ii)==RE);

    f = @(x) (DI(ii) * 1e-13) * exp( -GG(ii) * (1-(x/(C0))) );    
    hold on
    %plot(XX, f(XX) ,'LineWidth', 2, 'DisplayName',['Re = ', num2str(round(RE(ii),2)), ', $\rho$ =', num2str(RHO(ii)), '$[kg/m^3]$'])
    plot(XX, f(XX) ,'LineWidth', 5, 'DisplayName', ['$\delta_H$ = ', num2str(round(delta(ii),2))])
    hold off
    text(XX(1)+0.01, f(XX(1)), ['$\delta_H$ = ', num2str(round(delta(ii),2)),'$[MPa^{1/2}]$, Re = ', num2str(round(RE(ii),2)), '[-]'], 'FontSize',28 );
end

%legend('Location','northwest')
%legend box off

xlim([0 1.6])
ylabel('$D_i^R \gamma(\Upsilon, C_s) [m^2/s] $')
xlabel('Normalized $C_s [kg/m^3]$')
set(gca,'FontSize',42)
%exportgraphics(figure(1), ['Gamma_function.png'], "Resolution",300);
close all
%}

%% Fit a surface in RE and F space
%{\
figure()
tiledlayout(1,2)

pbaspect([2 1 1])
[ml3, sf3] = fit([RE', FF'],DI','poly11');
h3 = plot(ml3, [RE', FF'],DI');
set(h3,'linestyle','none');
alpha(.5)

p3 = coeffvalues(ml3);
title(sprintf('$D_i^R = %.3f %.3f \\cdot Re + %.3f \\cdot F$ \n $R^2 = %.3f $',[p3, sf3.rsquare]))
xlabel('Re[-]', 'rotation', 15)
ylabel('$F \cdot 10^{-5}$ kg/s', 'rotation',-25, 'Position', [-0.1 3.0])
zlabel('$D_i^R \cdot 10^{-13}~ m^2/s$')
axis tight
grid off
view([-37.5 30]);
set(gca,'FontSize',20)
exportgraphics(figure(1), ['Di_Re_F.png'], "Resolution",300); close all

figure()
tiledlayout(1,2)

pbaspect([2 1 1])
[ml4, sf4] = fit([RE', FF'],GG','poly11');
h4 = plot(ml4, [RE', FF'],GG');
set(h4,'linestyle','none');
alpha(.5)

p4 = coeffvalues(ml4);
title(sprintf('$\\Upsilon ~~= %.3f + %.3f \\cdot Re %.3f \\cdot F$ \n $R^2 = %.3f $',[p4, sf4.rsquare]))
xlabel('Re[-]', 'rotation', 15)
ylabel('$F \cdot 10^{-5}$ kg/s', 'rotation',-25, 'Position', [0 4.5])
zlabel('$\Upsilon[-]$')
axis tight
grid off
view([-37.5 30]);
set(gca,'FontSize',20)
exportgraphics(figure(1), ['Gamma_Re_F.png'], "Resolution",300); close all





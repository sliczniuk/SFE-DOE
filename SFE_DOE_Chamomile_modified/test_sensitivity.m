startup;
delete(gcp('nocreate'));
% %p = Pushbullet(pushbullet_api);

%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

excel_file = 'Chamomile_Di_Gamma_2.xls';
%rng(69)

%%
Parameters_table        = readtable('Parameters.csv') ;                     % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});                  % Parameters within the model + (m_max), m_ratio, sigma

LabResults              = xlsread('dataset_2.xlsx');
%Q_RBF = readmatrix('Q_RBF.txt');
%Q_FP  = readmatrix('Q_FP.txt');

%% Create the solver
Iteration_max               = 0;                                         % Maximum number of iterations for optimzer
%Time_max                    = 17;                                           % Maximum time of optimization in [h]

nlp_opts                    = struct;
%nlp_opts.ipopt.max_iter     = Iteration_max;
nlp_opts.ipopt.max_cpu_time = Time_max*3600;
nlp_opts.ipopt.hessian_approximation ='limited-memory';

%% Load paramters
m_total                 = 3.0;

% Bed geometry
before                  = 0.04;                                             % Precentage of length before which is empty
bed                     = 0.92;                                              % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 600;
timeStep                = 5;                                                % Minutes
OP_change_Time          = 10; 
Sample_Time             = 5;

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                    % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
Time                    = [0 Time_in_sec/60];                               % Minutes

N_Time                  = length(Time_in_sec);

SAMPLE                  = Sample_Time:Sample_Time:ExtractionTime;

OP_change               = 0:OP_change_Time:ExtractionTime-OP_change_Time;

% Check if the number of data points is the same for both the dataset and the simulation
N_Sample                = [];
for i = 1:numel(SAMPLE)
    N_Sample            = [N_Sample ; find(round(Time,3) == round(SAMPLE(i))) ];
end
if numel(N_Sample) ~= numel(SAMPLE)
    keyboard
end

%% Specify parameters to estimate
nstages                 = Parameters{1};

nstagesbefore           = 1:floor(before*nstages);
nstagesbed              = nstagesbefore(end)+1 : nstagesbefore(end) + floor(bed*nstages);
nstagesafter            = nstagesbed(end)+1:nstages;

bed_mask                = nan(nstages,1);
bed_mask(nstagesbefore) = 0;
bed_mask(nstagesbed)    = 1;
bed_mask(nstagesafter)  = 0;

%% Number of variables
Nx_FP                   = (1+3)*(3 * nstages+3);                                % 3*Nstages(C_f, C_s, H) + P(t) + yield
Nu                      = 3 + numel( Parameters );                          % T_in, P, F + numel(Parameters)

%% Extractor geometry
r                       = Parameters{3};                                    % Radius of the extractor  [m]
epsi                    = Parameters{4};                                    % Fullness [-]
dp                      = Parameters{5};                                    % Paritcle diameter
L                       = Parameters{6};                                    % Total length of the extractor [m]

L_nstages               = linspace(0,L,nstages);
V                       = L  * pi * r^2;                                    % Total volume of the extractor [m3]
A                       = pi *      r^2;                                    % Extractor cross-section

%--------------------------------------------------------------------
V_slice                 = (L/nstages) * pi * r^2;

V_before                = V_slice * numel(nstagesbefore);
V_after                 = V_slice * numel(nstagesafter);
V_bed                   = V_slice * numel(nstagesbed);                      % Volume of the fixed bed [m3]

V_before_solid          = repmat(V_before * 0          / numel(nstagesbefore), numel(nstagesbefore),1);
V_bed_solid             = repmat(V_bed    * epsi       / numel(nstagesbed)   , numel(nstagesbed)   ,1);
V_after_solid           = repmat(V_after  * 0          / numel(nstagesbed)   , numel(nstagesafter) ,1);

V_solid                 = [V_before_solid; V_bed_solid; V_after_solid];

V_before_fluid          = repmat(V_before * 1          / numel(nstagesbefore), numel(nstagesbefore),1);
V_bed_fluid             = repmat(V_bed    * (1 - epsi) / numel(nstagesbed)   , numel(nstagesbed)   ,1);
V_after_fluid           = repmat(V_after  * 1          / numel(nstagesafter) , numel(nstagesafter) ,1);

V_fluid                 = [V_before_fluid; V_bed_fluid; V_after_fluid];

L_bed_after_nstages     = L_nstages(nstagesbed(1):end);
L_bed_after_nstages     = L_bed_after_nstages - L_bed_after_nstages(1);
L_end                   = L_bed_after_nstages(end);

%% symbolic variables
x_FP                    = MX.sym('x', Nx_FP);
u                       = MX.sym('u', Nu);

%% Set inital state and inital conditions
msol_max                = m_total;                                          % g of product in solid and fluid phase
mSol_ratio              = 1;

mSOL_s                  = msol_max*mSol_ratio;                              % g of product in biomass
mSOL_f                  = msol_max*(1-mSol_ratio);                          % g of biomass in fluid

C0solid                 = mSOL_s * 1e-3 / ( V_bed * epsi)  ;                % Solid phase kg / m^3
Parameters{2}           = C0solid;

G                       =@(x) -(2*mSOL_f / L_end^2) * (x-L_end) ;

m_fluid                 = G(L_bed_after_nstages)*( L_bed_after_nstages(2) ); % Lienarly distirubuted mass of solute in fluid phase, which goes is zero at the outlet. mass*dz
m_fluid                 = [zeros(1,numel(nstagesbefore)) m_fluid];
C0fluid                 = m_fluid * 1e-3 ./ V_fluid';

%% Set Integrator - RBF
f_corr                   = @(x, u) modelSFE_Corr_sensitivity(x, u, bed_mask, timeStep_in_sec);

% Integrator
[F_corr]                 = buildIntegrator(f_corr, [Nx_FP,Nu] , timeStep_in_sec);

%%

% Set operating conditions
%{\
OPT_solver                  = casadi.Opti();
ocp_opts                    = {'nlp_opts', nlp_opts};
OPT_solver.solver(             'ipopt'   , nlp_opts)

T0homog                 = OPT_solver.variable(numel(OP_change))';
                          OPT_solver.subject_to( 30+273 <= T0homog <= 40+273 );
    
Flow                    = OPT_solver.variable(numel(OP_change))';
                          OPT_solver.subject_to( 3.33 <= Flow <= 6.67 );

feedPress               = 100;               % MPa -> bar

%Pressure                = OPT_solver.variable(numel(OP_change_P))';
%                          OPT_solver.subject_to( 100 <= Pressure <= 200 );

feedTemp                = repmat(T0homog,OP_change_Time/timeStep,1);
feedTemp                = feedTemp(:)';
feedTemp                = [ feedTemp, T0homog(end)*ones(1,N_Time - numel(feedTemp)) ];   
T_0                     = feedTemp(1)       ;

feedPress               = feedPress * ones(1,numel(feedTemp)) + 0 ;     % Bars
%feedPress                = repmat(Pressure,OP_change_Time_P/timeStep,1);
%feedPress                = feedPress(:)';
%feedPress                = [ feedPress, Pressure(end)*ones(1,N_Time - numel(feedPress)) ];    

feedFlow                = repmat(Flow,OP_change_Time/timeStep,1);
feedFlow                = feedFlow(:)';
feedFlow                = [ feedFlow, Flow(end)*ones(1,N_Time - numel(feedFlow)) ];    
F_0                     = feedFlow(1) ;
%}

%{
T0homog                 = 35+273;                    % K
feedPress               = 150;                       % bar
Flow                    = 5 ;                        % kg/s
feedTemp                = T0homog   * ones(1,length(Time_in_sec)) + 0 ;     % Kelvin
feedTemp                = linspace(30+273,40+273,N_Time);
%T0homog                 = feedTemp(1)  ;

feedPress               = feedPress * ones(1,N_Time) + 0 ;     % Bars
    
feedFlow                = Flow * ones(1,length(Time_in_sec));               % kg/s
%feedFlow                = linspace(3.3,6.6, N_Time);
%}

Z                       = Compressibility( T0homog, feedPress,         Parameters );
rho                     = rhoPB_Comp(      T0homog, feedPress, Z,      Parameters );
    
enthalpy_rho            = rho.*SpecificEnthalpy(T0homog, feedPress, Z, rho, Parameters );

uu                      = [feedTemp', feedPress', feedFlow'];

Z                       = Compressibility(      feedTemp(1), feedPress(1),         Parameters );
rho                     = rhoPB_Comp(           feedTemp(1), feedPress(1), Z,      Parameters );
enthalpy_rho            = rho.*SpecificEnthalpy(feedTemp(1), feedPress(1), Z, rho, Parameters ) ;

x0_FP                   = [ C0fluid'                         ;
                            C0solid         * bed_mask       ;
                            enthalpy_rho    * ones(nstages,1);
                            feedPress(1)                     ;
                            0                                ;
                            0                                ;
                            zeros(Nx_FP  - (3*nstages+3), 1) ;
                            ];

Parameters_init_time   = [uu repmat(cell2mat(Parameters),1,N_Time)'];

[xx_F]                 = simulateSystem(F_corr , [], x0_FP, Parameters_init_time );

%%
sigma = 0.01;

Q_FP_1 = readmatrix('Q_FP_1.txt');
Q_FP_1 = Q_FP_1(1:3,1:3);
Q_FP_2 = readmatrix('Q_FP_2.txt');
Q_FP_2 = Q_FP_2(1:3,1:3);

XX_FP    = xx_F(3*nstages+2:3*nstages+3:Nx_FP,:);
yy_FP    = XX_FP(1,N_Sample);
S_FP     = XX_FP(2:end,N_Sample) ;

SS = S_FP * S_FP' ;
JJ = -log(myDet(SS./sigma + pinv( Q_FP_1 ) ));

%%
T0 = ( (30-40).*rand(1,numel(Flow)) + 40 ) + 273;
F0 = ( (3.33-6.67).*rand(1,numel(Flow)) + 3.33 );

%{\
%% Solve the optimization problem
OPT_solver.minimize(JJ);
OPT_solver.set_initial([T0homog, Flow], [T0, F0] );

%%
tic
try
    sol  = OPT_solver.solve();
    KOUT = full(sol.value([T0homog, Flow])) ;
catch
    KOUT = OPT_solver.debug.value([T0homog, Flow]);
end
toc

%%
AA     = Function('AA', {[T0homog,Flow]}, {xx_F} );
res_0  = full(AA([T0, F0]));
res    = full(AA(KOUT));

%%
%{'
subplot(4,1,1); 
hold on; 
plot(Time, res_0(3*nstages+2,:)); 
plot(Time, res(3*nstages+2,:)); 
hold off; 

legend('Inital yield curve', 'Optimized yield curve', 'box', 'off')
subplot(4,1,2)
hold on;
stairs(OP_change, T0-273)
stairs(OP_change, KOUT(1:numel(OP_change))-273)
hold off

subplot(4,1,3)
hold on;
stairs(OP_change, F0)
stairs(OP_change, KOUT(numel(OP_change)+1:end))
hold off
legend('Inital', 'Final', 'box', 'off')

subplot(4,1,4)
plot([0:OPT_solver.stats.iter_count], OPT_solver.stats.iterations.obj)

%}
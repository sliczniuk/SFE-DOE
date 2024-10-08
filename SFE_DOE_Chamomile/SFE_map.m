%startup;
%delete(gcp('nocreate'));
%p = Pushbullet(pushbullet_api);
initParPool

addpatch('casadi_folder')
%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
%addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

%% Create the solver
tic
Iteration_max               = 100;                                           % Maximum number of iterations for optimzer
Time_max                    = 12;                                             % Maximum time of optimization in [h]

nlp_opts                    = struct;
nlp_opts.ipopt.max_iter     = Iteration_max;
nlp_opts.ipopt.max_cpu_time = Time_max*3600;
%nlp_opts.ipopt.hessian_approximation ='limited-memory';
%nlp_opts.ipopt.acceptable_tol  = 1e-4;
%nlp_opts.ipopt.acceptable_iter = 5;

OPT_solver                  = casadi.Opti();
ocp_opts                    = {'nlp_opts', nlp_opts};
OPT_solver.solver(             'ipopt'   , nlp_opts)

%%
Parameters_table        = readtable('Parameters.csv') ;                     % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});                  % Parameters within the model + (m_max), m_ratio, sigma

LabResults              = xlsread('wpd_datasets.xlsx');
which_dataset           = 4;

SAMPLE                  = LabResults(21:34,1);
data_org                = LabResults(21:34,which_dataset+1)';

%% Load paramters
N_exp                   = 12;
N_core                  = 14;

m_total                 = 3.0;

% Bed geometry
before                  = 0.04;                                             % Precentage of length before which is empty
bed                     = 0.92;                                              % Percentage of length occupied by fixed bed

% Set time of the simulation
PreparationTime         = 0;
ExtractionTime          = 300;
timeStep                = 10;                                                % Minutes
SamplingTime            = 10;                                                % Minutes
OP_change_Time          = 10;                                                % Minutes

simulationTime          = PreparationTime + ExtractionTime;

timeStep_in_sec         = timeStep * 60;                                    % Seconds
Time_in_sec             = (timeStep:timeStep:simulationTime)*60;            % Seconds
Time                    = [0 Time_in_sec/60];                               % Minutes

N_Time                  = length(Time_in_sec);
SAMPLE                  = SamplingTime:SamplingTime:ExtractionTime;
OP_change               = OP_change_Time:OP_change_Time:ExtractionTime;

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
Nx                      = 3 * nstages+2;                                    % 3*Nstages(C_f, C_s, H) + P(t) + yield
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
x                       = MX.sym('x', Nx);
u                       = MX.sym('u', Nu);

%% Set Integrator
f                       = @(x, u) modelSFE(x, u, bed_mask, timeStep_in_sec);
xdot                    = modelSFE(x, u, bed_mask, timeStep_in_sec);

% Integrator
F                       = buildIntegrator(f, [Nx,Nu] , timeStep_in_sec);

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

%% Set operating conditions
%T0homog                 = LabResults(1,which_dataset+1) * ones(1,numel(OP_change));                    % K

T0homog                 = OPT_solver.variable(numel(OP_change))';
                          OPT_solver.subject_to( 30+273 <= T0homog <= 40+273 );
%T0homog                 = linspace(30,40,numel(OP_change))+273;

feedPress               = 150;                                               % bar

%Flow                    = LabResults(3,which_dataset+1) * 1e-5 ;            % kg/s
Flow                    = OPT_solver.variable(numel(OP_change))';
                          OPT_solver.subject_to( 3.33 <= Flow <= 6.67 );

%feedTemp                = T0homog   * ones(1,length(Time_in_sec)) + 0 ;     % Kelvin

feedTemp                = repmat(T0homog,OP_change_Time/timeStep,1);
feedTemp                = feedTemp(:)';
feedTemp                = [ feedTemp, T0homog(end)*ones(1,N_Time - numel(feedTemp)) ];    
T_0                     = feedTemp(1);   

feedPress               = feedPress * ones(1,length(Time_in_sec)) + 0 ;     % Bars

Z                       = Compressibility( T_0, feedPress(1),         Parameters );
rho                     = rhoPB_Comp(      T_0, feedPress(1), Z,      Parameters );
enthalpy_rho            = rho.*SpecificEnthalpy(T_0, feedPress(1), Z, rho, Parameters ) ;

%feedFlow                = Flow * ones(1,length(Time_in_sec));               % kg/s
feedFlow                = repmat(Flow,OP_change_Time/timeStep,1);
feedFlow                = feedFlow(:)';
feedFlow                = [ feedFlow, Flow(end)*ones(1,N_Time - numel(feedFlow)) ];    

uu                      = [feedTemp', feedPress', feedFlow'];

% Initial conditions
x0                      = [ C0fluid'                         ;
                            C0solid         * bed_mask       ;
                            enthalpy_rho    * ones(nstages,1);
                            feedPress(1)                     ;
                            0                                ;
                            ];

%%
%{\
sigma                   = 0.12;
which_theta             = [44:46];

% Store symbolic results of the simulation
Parameters_sym          = MX(cell2mat(Parameters));
Parameters_sym_t        = OPT_solver.parameter(numel(which_theta));
Parameters_sym(which_theta)   = Parameters_sym_t;

X                       = MX(Nx,N_Time+1);
X(:,1)                  = x0;
% Symbolic integration
for jj=1:N_Time
    X(:,jj+1)=F(X(:,jj), [uu(jj,:)'; Parameters_sym] );
end

Yield_estimate          = X(Nx,[1; N_Sample]);
data_obs                = diff(Yield_estimate);

S                       = jacobian(data_obs, Parameters_sym_t);

S_norm                  = S' .* abs(cell2mat(Parameters(which_theta))) ./ repmat(data_obs,numel(which_theta) );
FI                      = (S_norm * S_norm');
FI                      = FI ./ (sigma^2);  

if numel(FI) ~= numel(which_theta)^2
    keyboard
end

D_opt = -log(myDet(FI));

GG = Function('GG',{T0homog,Flow,Parameters_sym_t}, {D_opt} )

%T0 = linspace(30,40,numel(T0homog))+273;
T0 = ( (40-30).*rand(1,numel(T0homog)*N_exp) + 30 ) + 273;

%F0 = linspace(3.33,6.67,numel(Flow));
F0 = ( (6.67-3.33).*rand(1,numel(Flow)*N_exp) + 3.33 );
%GG(T0,F0,cell2mat(Parameters(which_theta))')

OPT_solver.set_value(Parameters_sym_t, cell2mat(Parameters(which_theta)));
OPT_solver.minimize(D_opt);
%OPT_solver.set_initial([T0homog, Flow], [T0, F0] );
toc 
%%
%{
tic 
try
    sol  = OPT_solver.solve();
    KOUT = full(sol.value([T0homog, Flow])) 
catch
    KOUT = OPT_solver.debug.value([T0homog, Flow])
end
toc 
%}
%%
tic
H = OPT_solver.to_function('H',{T0homog,Flow},{T0homog,Flow});
toc
%%
tic
%[HH_1, HH_2] = H(T0,F0)
toc
%%
%parpool()
tic 
H_map = H.map(N_exp,'thread',N_core);
toc

%%
tic
[HH_map_1, HH_map_2] = full( H_map(T0, F0) );
toc

%%
tic
COST_0 = GG([T0, T0] ,  [F0, F0]  , cell2mat(Parameters(which_theta)) );
toc
tic 
COST   = GG(HH_map_1, HH_map_2, cell2mat(Parameters(which_theta)) );
toc 

save DOE_T.mat

writematrix([COST_0; COST], 'COST.txt')
writematrix([HH_map_1; HH_map_2], 'CONTROL.txt')






























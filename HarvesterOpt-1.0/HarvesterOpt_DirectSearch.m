% -----------------------------------------------------------------
%  HarvesterOpt_DirectSearch.m
% ----------------------------------------------------------------- 
%  This is the main file for a program which solve an optimization
%  problem that seaks to maximize the output power of a 
%  piezo-magneto-elastic beam, which evolves according to the 
%  follwing system of ordinary differential equations
%
%    d2x/dt2 + 2*ksi*dx/dt - 0.5*x*(1-x^2) - chi*v = f*cos(Omega*t)
%
%    dv/dt + lambda*v + kappa*dx/dt = 0
%
%        +
%
%    initial conditions,
%  
%  where
%  
%   x(t)   - dimensionless displacement of the beam tip
%   v(t)   - dimensionless voltage across the load resistance
%   t      - dimensionless time
%   ksi    - mechanical damping ratio
%   chi    - dimensionless piezoeletric coupling term (mechanical)
%   f      - dimensionless excitation amplitude
%   Omega  - dimensionless excitation frequency
%   lambda - dimensionless time constant reciprocal
%   kappa  - dimensionless piezoeletric coupling term (eletrical)
%  
%  Reference:
%  A. Cunha Jr
%  Enhancing the performance of a bistable energy harvesting 
%  device via the cross-entropy method (2020)
%  
% -----------------------------------------------------------------
%  programmer: Americo Cunha
%              americo@ime.uerj.br
%
%  last update: March 31, 2020
% -----------------------------------------------------------------

clc
clear
close all


% program header
% -----------------------------------------------------------
disp(' ---------------------------------------------------')
disp(' Piezo-Magneto-Elastic Harvester Optimization       ')
disp(' (direct search)                                    ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Cunha (UERJ)                               ')
disp(' americo@ime.uerj.br                                ')
disp(' ---------------------------------------------------')
% -----------------------------------------------------------



% simulation information
% -----------------------------------------------------------
case_name = 'HarvesterOpt_DS_case1';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------



% define physical parameters
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- defining model parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

ksi    = 0.01;  % mechanical damping ratio
chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)
lambda = 0.05;  % dimensionless time constant reciprocal
kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)
f      = 0.0;   % dimensionless excitation amplitude
Omega  = 0.0;   % dimensionless excitation frequency

x0    = 1.0;  % dimensionless initial displacement
xdot0 = 0.0;  % dimensionless initial velocity
v0    = 0.0;  % dimensionless initial voltage

toc
%---------------------------------------------------------------


% 0-1 test for chaos parameters
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- defining 0-1 test for chaos parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

% tolerance for 0-1 test
tol01 = 0.1;

% parameter c for 0-1 test
cmin = 0.0;
cmax = 2*pi;

% number of test repeats
Nc = 100;

% oversampling flag
OSflag = 0;

toc
% -----------------------------------------------------------


% define ODE solver parameters
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- defining ODE solver parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

% initial time of analysis
t0 = 0.0;

% final time of analysis
t1 = 2500.0;

% time interval of analysis
tspan = [t0 t1];

% inicial condition vector
IC = [x0; xdot0; v0];

% ODE solver optional parameters
%opt = odeset('RelTol',1.0e-9,'AbsTol',1.0e-6);

toc
% -----------------------------------------------------------



% solve optimization problem via direct search
%---------------------------------------------------------------
tic
disp(' '); 
disp(' --- optimization via direct search --- ');
disp(' ');
disp('    ... ');
disp(' ');

% number of points in each direction
Np1 = 10;
Np2 = 10;

% penalty parameter vector (K01 p1max p1min p2max p2min)
Hpenalty = 10.0;

% control parameter 1 (excitation amplitude)
p1_min = 0.08;
p1_max = 0.10;

% control parameter 2 (excitation frequency)
p2_min = 0.75;
p2_max = 0.85;

% vector of p1 samples 
p1 = linspace(p1_min,p1_max,Np1);

% vector of p2 samples
p2 = linspace(p2_min,p2_max,Np2);

% preallocate memory for PerfFunc
% (x = columns and y = lines)
S = zeros(Np2,Np1);

% preallocate memory for power function
% (x = columns and y = lines)
power = zeros(Np2,Np1);

% preallocate memory for 0-1 classifier
% (x = columns and y = lines)
K01 = zeros(Np2,Np1);

% initialize PerfFunc maximum/mininum value
S_max_overall = -Inf;
S_min_overall = +Inf;

% initialize power maximum/minimum value
power_max_overall = -Inf;
power_min_overall = +Inf;

% initialize K01 maximum/minimum value
K01_max_overall = -Inf;
K01_min_overall = +Inf;

% initialize power maximum/minimum points
p1_max_overall    = 0.0;
p2_max_overall    = 0.0;
p1_min_overall    = 0.0;
p2_min_overall    = 0.0;

% define data file name
file_name = [case_name,'.dat'];

% open data file
fileID = fopen(file_name,'w');

% define data format
formatSpec = '%.4f %.4f %+0.4E %0.4E %1.4f \n';

for np1 = 1:Np1
    
    % print loop indicador
    disp(['Np1 = ',num2str(np1)]);
    
    for np2 = 1:Np2
        
        % update parameters
        f     = p1(np1);
        Omega = p2(np2);
        
        % define physical paramters vector
        phys_param = [ksi chi f Omega lambda kappa];
                          
        % penalized performance function S(x)
        [S(np2,np1),power(np2,np1),K01(np2,np1)] = ...
          PiezoMagBeam_PerfFunc(phys_param,tspan,IC,...
                                cmin,cmax,Nc,tol01,OSflag,Hpenalty,...
                                p1_min,p1_max,p2_min,p2_max);
        
        % update global maximum
        if S(np2,np1) > S_max_overall
            
              S_max_overall = S(np2,np1);
          power_max_overall = power(np2,np1);
            K01_max_overall = K01(np2,np1);
             p1_max_overall = p1(np1);
             p2_max_overall = p2(np2);
        end
        
        % update global minimum
        if S(np2,np1) < S_min_overall
            
             S_min_overall = S(np2,np1);
         power_min_overall = power(np2,np1);
           K01_min_overall = K01(np2,np1);
            p1_min_overall = p1(np1);
            p2_min_overall = p2(np2);
        end
        
        % print local values in data file
        fprintf(fileID,formatSpec,...
                 p1(np1),p2(np2),S(np2,np1),power(np2,np1),K01(np2,np1));
        
    end
end

% print global maximum
fprintf(fileID,formatSpec,...
        p1_max_overall,p2_max_overall,...
        S_max_overall,power_max_overall,K01_max_overall);
    
% print global minimum
fprintf(fileID,formatSpec,...
        p1_min_overall,p2_min_overall,...
        S_min_overall,power_min_overall,K01_min_overall);

% close data file
fclose(fileID);

% display global extreme values on screen
fprintf('\n\nGlobal maximum (f Omega PerfFunc power K01):\n\n');

fprintf(formatSpec,...
        p1_max_overall,p2_max_overall,...
        S_max_overall,power_max_overall,K01_max_overall);
    
fprintf('\n\nGlobal minimum (f Omega PerfFunc power K01):\n\n');

fprintf(formatSpec,...
        p1_min_overall,p2_min_overall,...
        S_min_overall,power_min_overall,K01_min_overall);

fprintf('\n');

disp('   ');
time_elapsed = toc
%---------------------------------------------------------------



% save simulation results
% -----------------------------------------------------------
tic
disp(' ')
disp(' --- saving simulation results --- ');
disp(' ');
disp('    ... ');
disp(' ');

save([num2str(case_name),'.mat']);

toc
% -----------------------------------------------------------


% post processing
% -----------------------------------------------------------
tic
disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');


% plot performance function countourmap
% ......................................................
gtitle = ' performance function contour map';
xlab   = ' excitation amplitude';
ylab   = ' excitation frequency';
xmin   = p1_min;
xmax   = p1_max;
ymin   = p2_min;
ymax   = p2_max;
gname  = [case_name,'__perf_func'];
flag   = 'eps';
fig1   = graph_contour_pnt(p1,p2,S',...
                            p1_max_overall,p2_max_overall,...
                            gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1);
% ......................................................


% plot mean power countourmap
% ......................................................
gtitle = ' mean power contour map';
xlab   = ' excitation amplitude';
ylab   = ' excitation frequency';
xmin   = p1_min;
xmax   = p1_max;
ymin   = p2_min;
ymax   = p2_max;
gname  = [case_name,'__mean_power'];
flag   = 'eps';
fig2   = graph_contour_pnt(p1,p2,power',...
                            p1_max_overall,p2_max_overall,...
                            gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);

%close(fig2);
% ......................................................


% plot 0-1 identification parameter
% ......................................................
xlab   = 'excitation amplitude';
ylab   = 'excitation frequency';
label0 = 'regular';
label1 = 'chaos';
gtitle = '0-1 test for chaos';
xmin   = p1_min;
xmax   = p1_max;
ymin   = p2_min;
ymax   = p2_max;
zmin   = 0;
zmax   = 1;
gname  = [case_name,'__01test'];
flag   = 'eps';
fig3   = graph_binarymap(p1,p2,K01,...                          
                          p1_max_overall,p2_max_overall,...
                          gtitle,xlab,ylab,label0,label1,...
                          xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig3);
% ......................................................

toc
% -----------------------------------------------------------

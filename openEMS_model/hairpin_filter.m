%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hairpin filter openEMS model
%
% Simulation time: >30 min
%
% Article at:
% http://charleslabs.fr/en/project-Hairpin+filter+design
%
% Created with
%  - Octave 5.2.0
%  - openEMS v0.0.35
%
% (c) 2020 Charles Grassin
% This code is licensed under MIT license.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc
physical_constants;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% START OF CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hairpin setup
hairpin.length = 38;    % Elements length
hairpin.feed = 10;      % Feeder length
hairpin.margin = 10;    % Left/Right/Top margin
hairpin.width = 2.54;   % Width of the traces
hairpin.D2_D3 = 2;      % D2 and D3 distances
hairpin.D1 = 0.254;     % D1 distance

% Substrate setup
substrate.thickness = 1.6;  % Thickness of substrate
substrate.epsilon   = 4.34; % Diaelectric constant
substrate.mu   = 1;         % Magnetic permeability
substrate.cells = 4;        % Number of cells for meshing substrate
substrate.kappa  = 1e-3 * 2*pi*2.45e9 * EPS0 * substrate.epsilon; % Electric conductivity
substrate.x_size  = hairpin.width*8+hairpin.D2_D3*5+hairpin.D1*2+hairpin.margin*2;
substrate.y_size = hairpin.length + hairpin.feed + hairpin.margin;

% Copper setup
copper.no_losses = true;    % Boolean : use perfect conductor?
copper.conductivity = 56e6; % Metal conductivity
copper.thickness = 35e-6;   % Metal thickness in m

% Feeding setup
feed.R = 50;         % Feed resistance
feed.f_center = 1e9; % Feed center frequency
feed.f_cut = 250e6;  % Feed -20 dB cut frequency

% Misc setup
show_AppCSXCAD = 1; %open AppCSXCAD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% END OF CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup the simulation
unit = 1e-3; % all length in mm

%% Simulation box sizing
SimBox = [substrate.x_size substrate.y_size substrate.thickness*5];

%% Setup FDTD parameter & excitation function
f_max = feed.f_center+feed.f_cut;
f_min = feed.f_center-feed.f_cut;
FDTD = InitFDTD('NrTS',  400000 );
FDTD = SetGaussExcite( FDTD, feed.f_center, feed.f_cut );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'PEC' 'MUR'}; % boundary conditions (PEC on zmin is GND plane)
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
CSX = InitCSX();

% Initialize the mesh
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [0 SimBox(3)];

% -------------------------------

%% Create substrate
CSX = AddMaterial( CSX, 'FR4');
CSX = SetMaterialProperty( CSX, 'FR4', 'Epsilon',substrate.epsilon, ...
                           'Mue', substrate.mu, 'Kappa', substrate.kappa);
start = [-substrate.x_size/2  -substrate.y_size/2                    0];
stop  = [ substrate.x_size/2   substrate.y_size/2  substrate.thickness];
CSX = AddBox( CSX, 'FR4', 1, start, stop );
% add extra cells to discretize the substrate thickness
mesh.z = [linspace(0,substrate.thickness,substrate.cells+1) mesh.z];

%% Create metal
if(copper.no_losses)
  CSX = AddMetal( CSX, 'Copper' ); % create a perfect electric conductor (PEC)
else
  CSX = AddConductingSheet( CSX, 'Copper', copper.conductivity , copper.thickness ) % Lossy conductor
end

% -------------------------------

%% --- CREATE HAIRPIN STRUCTURE ----
pt_feed_port = [0 hairpin.width hairpin.width 0;0 0 hairpin.length+hairpin.feed hairpin.length+hairpin.feed];
pt_U = [0 0 hairpin.width hairpin.width+hairpin.D2_D3 hairpin.width*2+hairpin.D2_D3 hairpin.width*2+hairpin.D2_D3 hairpin.width+hairpin.D2_D3 hairpin.width+hairpin.D2_D3 hairpin.width hairpin.width;0 hairpin.length hairpin.length+hairpin.width hairpin.length+hairpin.width hairpin.length 0 0 hairpin.length hairpin.length 0];
pt_U_inverted = [0 hairpin.width hairpin.width+hairpin.D2_D3 hairpin.width*2+hairpin.D2_D3 hairpin.width*2+hairpin.D2_D3 hairpin.width+hairpin.D2_D3 hairpin.width+hairpin.D2_D3 hairpin.width hairpin.width 0;0 -hairpin.width -hairpin.width 0 hairpin.length hairpin.length 0 0 hairpin.length hairpin.length];

% translate vector for relative positionning
tl = [-substrate.x_size/2+hairpin.margin;-substrate.y_size/2];  

% Hairpin Feed
CSX = AddPolygon( CSX, 'Copper', 2, 2, substrate.thickness, pt_feed_port + tl);
%% Excitation port (1)
start = [0 0 0] + [tl(1) tl(2) substrate.thickness];
stop  = start + [hairpin.width 1 0];
[CSX port{1}] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 1 0], true);
% translate
tl = tl + [hairpin.width+hairpin.D1;hairpin.feed]; 
% Hairpin first U
CSX = AddPolygon( CSX, 'Copper', 2, 2, substrate.thickness, pt_U + tl);
% translate
tl = tl + [hairpin.width*2+hairpin.D2_D3*2;0];
% Hairpin second U
CSX = AddPolygon( CSX, 'Copper', 2, 2, substrate.thickness, pt_U_inverted + tl);
% translate
tl = tl + [hairpin.width*2+hairpin.D2_D3*2;0];
% Hairpin third U
CSX = AddPolygon( CSX, 'Copper', 2, 2, substrate.thickness, pt_U + tl);
% translate
tl = tl + [hairpin.width*2+hairpin.D2_D3+hairpin.D1;-hairpin.feed];
% Hairpin output
CSX = AddPolygon( CSX, 'Copper', 2, 2, substrate.thickness, pt_feed_port + tl);
% Output port (2)
start = [0 0 0] + [tl(1) tl(2) substrate.thickness];
stop  = start + [hairpin.width 1 0];
[CSX port{2}] = AddLumpedPort(CSX, 5 ,2 ,feed.R, start, stop, [0 -1 0], false);

% -------------------------------

% Measure Box (for Paraview)
CSX = AddDump(CSX, 'Et');
CSX = AddBox(CSX, 'Et', 0, [min(mesh.x) min(mesh.y) min(mesh.z)], [max(mesh.x) max(mesh.y) max(mesh.z)]);

% -------------------------------

%% finalize the mesh
% generate a smooth mesh with max. cell size: lambda_min / 20
mesh.x = SmoothMeshLines(mesh.x, 2);
mesh.y = SmoothMeshLines(mesh.y, 5);

mesh = DetectEdges(CSX, mesh);
mesh = SmoothMesh(mesh, c0 / f_max / unit / 20);
CSX = DefineRectGrid(CSX, unit, mesh);

% -------------------------------

%% prepare simulation folder
Sim_Path = 'tmp_Hairpin';
Sim_CSX = 'Hairpin.xml';

try confirm_recursive_rmdir(false,'local'); end
 
[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

% -------------------------------

%% show the structure
if (show_AppCSXCAD == 1)
  CSXGeomPlot( [Sim_Path '/' Sim_CSX] );
end

%% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);

% -------------------------------

%% POST PROCESSING
% Plot s11 and s21
figure
f = linspace( f_min, f_max, 1000 );
p = calcPort( port, Sim_Path, f, 'RefImpedance', feed.R);
s11 = p{1}.uf.ref./ p{1}.uf.inc;
s21 = p{2}.uf.ref./ p{1}.uf.inc;

plot(f/1e9,20*log10(abs(s11)),'k-','LineWidth',2);
hold on;
grid on;
plot(f/1e9,20*log10(abs(s21)),'r--','LineWidth',2);
title( 'S_{11} and S_{21} parameters' );
legend('S_{11}','S_{21}');
ylabel('S-Parameter (dB)','FontSize',12);
xlabel('frequency (GHz) \rightarrow','FontSize',12);
ylim([-40 2]);
drawnow
clc; clear; close all;
syms s K1 K2 K3 K4;

Jm = 1.03 * 10^-4 ;% kgm^2
Jp = 3.10 * 10^-4 ;% kgm^2 
Bm = 0.268;% Nms/rad
Brp = Bm * 80 / 100 ;

Jeq = Jm + Jp; % kgm^2 
Beq = Brp + Bm; % Nms/rad % Friction considered as %80 of Bm and added to the Bm and found Beq
Kb = 0.242; % Vs/rad
Km = 0.242; % Nm/A
L = 0.000632; % Hennry
R = 0.218; % Ohm

A = [-Beq/Jeq, Km/Jeq; -Kb/L, -R/L];
B = [0; 1/L];
C = [1, 0];
D = 0;

motor = ss(A,B,C,D);
gear_ratio = tf(1,[5,0]);
plant = motor * gear_ratio;
%step(plant)
A = plant.A;
B = plant.B;
C = [0,0,1];
D = 0;

plant = ss(A,B,C,D,'Inputname','v', 'Outputname', 'y', 'Statename',{'w','i','x'});
%characteristicEquation 

characteristicEquation = det(s * eye(size(A)) - A);
poles = eig(A);
%desired poles and new characteristic equation
F = [0 1 0;0 0 1;0 -6.2727e+05 -1.5130e+03];
G = [0;0;1];

DCE = det(s*eye(3)- F + G * [K1 K2 K3]);% = s^3 + 150s^2 + 7.78 s + 1167
%DCE = K1 + (627270 + K2)s + (K3 + 1513)s^2 + s^3
%K1 = 1167, K2 = -627262.22, K3 = -1363

ctrb_matrix = ctrb(plant);
iscontrollable = det(ctrb_matrix); % Different than zero which means controllable

pc = [-4+2.79i, -4-2.79i,-150];
K = place(A,B,pc);
Ntilda = rscale(plant, K);

Acl = A - B*K;
Bcl = B*Ntilda;
syscl = ss(Acl, Bcl, C, D);
%step(syscl)

%Integral Action
Abar = [plant.A, zeros(3,1);-plant.C, 0];
Bbar = [plant.B;0];
newcharacteristicEquation = det(s * eye(size(Abar)) - Abar - Bbar * [K1 K2 K3 K4]);
eig(Abar);
pc_integral = [-4+2.79i, -4-2.79i,-150 ,-200];
Kint = acker(Abar, Bbar, pc_integral);
F = Kint(1:3); % Closed loop coefficient 
H = Kint(4);% Integral gain


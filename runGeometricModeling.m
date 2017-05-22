% This script generates a struct TACR (Tendon Actuated Continuum Robot)
% which contains specific parameters of the robot, it executes the geometric
% model and includes a function which can be used to visualize the robot.
%
% Copyright: 2016 Leibniz Universit�t Hannover, All rights reserved
%
% Email: continuumrobotics@lkr.uni-hannover.de
%
% Version: 1
% Date: 11/16/2016
%

% 2 segment tendon actuated continuum robot
TACR.ndisks = [10;10];              % number of disks per segment
TACR.diskRadius = [8;8];            % disk radius
TACR.diskHeight = 3;                % heigth of the disks
TACR.diskPitchRadius = [6.5;5];     % pitch circle radius of disks
TACR.segmentLength = [92;102];      % segment length = length of first backbone per segment

q_1 = [-3 1.5 1.5; -4 2 2];        % actuation parameters (delta l per tendon);
q_2 = [-3 1.5 1.5; 2 -1 -1];       % Remember: only 2 tendons can be
q_3 = [0 0 0; -4 2 2];             % retracted at once, the 3rd tendon has
q_4 = [0 0 0; -2.5 -2.5 5];        % retracted at once, the 3rd tendon has
q_5 = [3 -2 -1; -3 1 2];           % to extend
q_6 = [3.6 -1.2 -2.4; 3 -1 -2];
q_7 = [3.6 -1.2 -2.4; 0 0 0];
q_8 = [-3 1.5 1.5; 0 0 0];
q_9 = [-3 1.5 1.5; -1 -2 3];
q_10 = [1.4 -2.2 0.8; 1.6 0.8 -2.4];
q_0 = [0 0 0;0 0 0];


dP_1= [0.00	0.00	3.00;
0.44	0.05	13.71;
1.38	0.39	23.60;
2.66	0.19	32.79;
4.78	0.56	43.50;
6.57	2.62	53.13;
10.14	2.08	62.43;
12.93	1.54	71.59;
17.00	1.73	81.75;
20.98	1.73	90.51;
25.43	1.50	99.50;
30.13	2.67	109.01;
35.43	1.25	117.73;
40.60	1.99	126.67;
46.57	1.18	134.93;
51.96	0.72	143.06;
57.09	1.82	151.17;
62.88	0.98	158.98;
70.25	3.85	166.24;
78.50	5.33	176.59];

dP_2= [0.00	0.00	3.00;
0.28	-0.12	13.62;
1.23	0.19	23.66;
2.42	0.58	32.81;
3.94	0.55	43.68;
6.04	0.60	53.29;
8.29	0.38	62.84;
10.81	0.75	71.89;
14.61	0.08	82.66;
18.40	0.19	91.38;
21.28	0.06	100.61;
24.48	0.43	111.42;
26.36	-0.28	121.21;
28.76	0.02	131.20;
29.53	0.88	141.59;
30.60	0.57	150.93;
31.09	0.47	160.98;
30.65	1.49	170.78;
30.07	0.40	180.89;
28.80	2.56	193.70;
];

dP_3=[0.00	0.00	3.00;
0.21	-0.10	13.68;
0.39	0.07	23.81;
1.07	0.71	32.96;
1.51	0.76	43.97;
2.23	1.65	53.76;
3.00	1.62	63.48;
4.35	1.88	73.13;
5.04	1.63	83.99;
6.30	1.87	94.02;
7.89	0.78	103.29;
10.51	1.11	113.97;
12.92	0.95	123.70;
16.20	-0.28	133.45;
19.76	-0.24	142.99;
23.88	-0.14	151.48;
27.18	-0.03	161.21;
32.71	-0.69	168.82;
38.06	-2.04	177.06;
46.69	-0.99	187.38;
];
dP_4=[0.00	0.00	3.00;
0.21	0.36	13.71;
-0.31	0.65	23.79;
0.30	1.17	33.03;
0.28	1.69	44.03;
0.63	2.10	53.82;
0.27	2.91	63.53;
0.59	3.41	73.15;
1.03	4.17	85.02;
1.29	4.36	93.81;
1.42	5.54	103.52;
2.62	6.52	114.55;
4.17	9.22	124.29;
5.74	12.18	134.22;
7.59	15.49	143.57;
9.78	19.30	151.98;
12.57	23.21	160.61;
14.87	27.68	168.70;
18.19	33.46	176.78;
22.71	40.27	185.53;
];
dP_5=[0.00	0.00	3.00;
0.09	0.42	13.71;
-0.51	0.86	23.79;
-0.90	1.61	33.06;
-2.26	2.26	44.01;
-3.59	3.07	53.73;
-4.98	3.69	63.42;
-7.03	4.52	72.93;
-9.45	5.21	83.84;
-11.44	7.04	93.51;
-13.58	7.38	102.61;
-15.01	7.87	113.55;
-15.52	9.13	123.49;
-15.19	9.19	133.83;
-14.31	10.59	143.92;
-13.26	11.12	153.21;
-11.72	11.96	162.91;
-9.03	12.63	172.28;
-5.42	13.66	181.86;
-1.40	16.18	192.94;
];
dP_6=[0.00	0.00	3.00;
-0.17	-0.03	13.71;
-1.45	0.39	23.71;
-3.52	0.40	32.74;
-5.37	-0.18	43.60;
-8.32	-0.32	52.96;
-11.14	-0.97	62.40;
-14.34	-1.24	71.46;
-18.36	-1.86	81.74;
-22.89	-2.53	90.39;
-27.96	-3.44	98.78;
-33.82	-4.59	108.17;
-38.76	-4.99	116.67;
-43.72	-6.13	125.73;
-48.98	-6.77	134.62;
-53.68	-7.49	142.62;
-59.41	-7.86	150.97;
-63.86	-9.26	159.45;
-69.99	-8.76	167.90;
-78.08	-8.47	178.08;
];
dP_7=[0.00	0.00	3.00;
-0.12	-0.01	13.58;
-1.28	-0.05	23.71;
-2.07	0.15	32.88;
-4.33	-0.38	43.76;
-5.97	-0.20	53.64;
-9.24	-1.21	62.76;
-11.85	-1.62	71.98;
-15.09	-1.61	82.55;
-18.85	-2.63	92.01;
-22.79	-3.65	100.59;
-26.68	-4.87	110.82;
-29.46	-5.45	120.53;
-32.67	-6.75	130.21;
-35.10	-8.94	140.12;
-37.64	-7.35	149.24;
-39.96	-8.01	158.99;
-41.35	-9.00	168.54;
-42.57	-8.33	178.62;
-44.22	-7.32	191.90;
];
dP_8=[0.00	0.00	3.00;
0.29	1.69	13.67;
1.06	1.52	23.75;
2.19	1.45	32.81;
4.06	1.50	43.59;
5.46	2.13	53.22;
7.80	1.74	62.81;
10.83	2.02	71.91;
13.96	2.08	82.30;
17.31	2.32	91.49;
20.81	1.25	100.71;
23.38	1.38	111.41;
26.80	1.92	120.99;
29.13	1.31	131.05;
31.82	1.39	141.04;
34.69	2.96	149.92;
34.10	2.94	160.55;
35.78	2.43	170.01;
36.10	3.84	180.03;
38.66	5.68	192.96;
];
dP_9=[0.00	0.00	3.00;
0.62	0.44	13.68;
1.62	1.41	23.65;
2.84	1.52	32.78;
4.65	1.99	43.64;
6.13	2.40	53.37;
8.52	4.07	62.77;
11.67	4.30	71.94;
15.66	4.03	82.41;
18.95	5.66	91.36;
22.05	6.55	100.60;
26.16	8.82	110.94;
28.62	10.04	120.52;
32.28	13.06	129.76;
35.08	17.15	139.10;
37.38	20.38	147.62;
39.70	24.12	156.76;
41.72	27.72	165.37;
43.64	31.84	174.51;
46.04	38.95	184.57;
];
dP_10=[0.00	0.00	3.00;
-0.37	-0.14	13.65;
-0.56	1.05	23.82;
-1.69	1.29	32.95;
-2.71	2.14	43.96;
-3.91	3.48	53.66;
-5.19	4.10	63.40;
-7.23	5.86	72.68;
-9.62	7.07	83.60;
-11.82	8.80	92.72;
-14.55	10.67	102.00;
-17.31	11.83	112.60;
-19.85	11.94	122.24;
-22.62	12.17	132.10;
-25.58	11.85	141.93;
-28.40	11.12	150.75;
-31.42	9.36	160.33;
-34.40	8.87	169.46;
-36.88	6.82	178.55;
-43.23	3.69	190.54;
];
q = q_10;
robotShape_mess.diskPoints = dP_10;

% compute robot's space curve                                    
[robotShape] = GeometricModel(TACR,q);

% visualize the robot
drawRobotOctave(robotShape,TACR);

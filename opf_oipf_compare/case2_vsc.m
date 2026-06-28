function mpc = case2_vsc
mpc.version = '2';
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm      Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1       3       0	0	0   0   1       1	0	345     1       1.1     0.9;
	2       1       0	0	0   0   1       1       0	345     1       1.1     0.9;
];

%% generator data
%	bus	Pg      Qg	Qmax	Qmin	Vg	mBase       status	Pmax	Pmin
mpc.gen = [
	1	0       0	500      -500    1	100       1       250     10;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status
mpc.branch = [
    1   2   0.015    0.12    0.00    100   100   100     0       0       1;
];
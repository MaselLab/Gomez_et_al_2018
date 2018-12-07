cd ~
cd MATLAB/

%% Simulations for new figures (2,4,5,6,7)

N = 1e9;
s = 1e-2;
u = 1e-5;
steps = 2e4;
times = 1:steps;

start_time = 5e3;
end_time = steps;
collect_data = true;
outputfile = '2dwave_data_time_series_stats_ml-01';

tic
[v,v1,v2,varx,vary,cov] = stochastic_simulation2(N,s,u,steps,collect_data,start_time,end_time,outputfile);
toc

%% Simulations for new figures (4)
N = 1e9;
s = 1e-2;
U = 1e-5;
steps = 1.5e6;

start_time = 1;
end_time = 1;
collect_data = false;
outputfile = '2dwave_data_time_series_stats_ml-01';

data_pts = 50;
x = ones(data_pts+1,1);

% sampling points for varying N, s and u (even spacing on log scale)
Narry = (1e7)*(1e11/1e7).^((0:1:data_pts)./data_pts);
sarry = (1e-3)*(1e-1/1e-3).^((0:1:data_pts)./data_pts);
Uarry = (1e-6)*(1e-4/1e-6).^((0:1:data_pts)./data_pts);

NsU = [Narry' s*x U*x; N*x sarry' U*x; N*x s*x Uarry'];
v = zeros(size(NsU,1),1);
vDF = zeros(size(NsU,1),1);
v1 = zeros(size(NsU,1),1);
v2 = zeros(size(NsU,1),1);
varx = zeros(size(NsU,1),1);
vary = zeros(size(NsU,1),1);
cov = zeros(size(NsU,1),1);

tic
for i=1:size(NsU,1)
    [v(i),v1(i),v2(i),varx(i),vary(i),cov(i)] = stochastic_simulation2(NsU(i,1),NsU(i,2),NsU(i,3),collect_data,steps,start_time,end_time,outputfile);
end
toc

csvwrite('~/Documents/kgrel2d/data/2dwave_data_time_avg_stats_ml-01-0.dat',NsU);
csvwrite('~/Documents/kgrel2d/data/2dwave_data_time_avg_stats_ml-01-1.dat',[v v1 v2 varx vary cov]);

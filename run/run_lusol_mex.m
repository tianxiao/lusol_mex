
%% generate random sparse matrix
rand('twister',0);
m = 5;
n = 5;
%density = 1;
%A = sprand(m,n,density);

A = rand(m,n);
A = sparse(A);

%% build input data
[i j s] = find(A);

m = int32(m);
n = int32(n);
nelem = int32(length(i));
lena = int32(300);

indc = zeros(lena,1,'int32');
indr = zeros(lena,1,'int32');
a = zeros(lena,1,'double');

indc(1:nelem) = int32(i);
indr(1:nelem) = int32(j);
a(1:nelem) = s;

luparm = zeros(30,1,'int32');
parmlu = zeros(30,1,'double');
ip = zeros(m,1,'int32');
iq = zeros(n,1,'int32');
lenc = zeros(n,1,'int32');
lenr = zeros(m,1,'int32');
locc = zeros(n,1,'int32');
locr = zeros(m,1,'int32');
iploc = zeros(n,1,'int32');
iqloc = zeros(m,1,'int32');
ipinv = zeros(m,1,'int32');
iqinv = zeros(n,1,'int32');
w = zeros(n,1,'double');
inform = zeros(1,1,'int32');

%% set the options
luparm(1) = 6;
luparm(2) = -1;
luparm(3) = 5;
luparm(6) = 0;
luparm(8) = 1;

parmlu(1) = 10.0;
parmlu(2) = 10.0;
parmlu(3) = eps^0.8;
parmlu(4) = eps^0.67;
parmlu(5) = eps^0.67;
parmlu(6) = 3.0;
parmlu(7) = 0.3;
parmlu(8) = 0.5;

%% make the call
lusol_mex('lu1fac',m,n,nelem,lena,luparm,parmlu,a,indc,indr,ip,iq,lenc,lenr,locc,locr,iploc,iqloc,ipinv,iqinv,w,inform);

%% extract U
lenU = luparm(24);
ui = zeros(lenU,1);
uj = zeros(lenU,1);
for i = 1:m
  ui(locr(i):(locr(i)+lenr(i)-1)) = i;
end
uj = double(indr(1:lenU));
ua = a(1:lenU);
U = sparse(ui,uj,ua)

%% extract L
% number of columns
numL0 = luparm(20)
% number of entries
lenL0 = luparm(21)

la = a(end-lenL0+1:end)
li = indc(end-lenL0+1:end)
lj = zeros(lenL0,1);
pos = 1;
for i = numL0:-1:1
    lj(pos:pos+lenc(i)-1) = i;
    pos = pos + lenc(i);
end

L = sparse(double(li),lj,la,double(m),double(n))







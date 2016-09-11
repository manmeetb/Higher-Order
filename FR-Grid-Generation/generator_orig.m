function generator(n,N)

%Number of elements in the x, y, and z directions

Lx = 2*pi();
Ly = 2*pi();
Lz = 2*pi();

fname = horzcat(horzcat(horzcat(num2str(n),'_'),num2str(N)),'.msh');

ni  = n;
nj  = n;
nk  = n;

npi = ni+1;
npj = nj+1;
npk = nk+1;

x   = zeros(npi*npj*npk,1);
y   = zeros(npi*npj*npk,1);
z   = zeros(npi*npj*npk,1);

BSSONIC = [];
BSYMM   = [];
BWALL   = [];

%Build the list of nodes
for i=1:npi
for j=1:npj
for k=1:npk
  x((i-1)*npj*npk+(j-1)*npk+k) = Lx/(npi-1)*((i-1))-Lx/2;
  y((i-1)*npj*npk+(j-1)*npk+k) = Ly/(npj-1)*((j-1))-Ly/2;
  z((i-1)*npj*npk+(j-1)*npk+k) = Lz/(npk-1)*((k-1))-Lz/2;
  
  if(i>1)
  if(j>1)
  if(k>1)
  if(i<npi)
  if(i<npj)
  if(i<npk)
      x((i-1)*npj*npk+(j-1)*npk+k) = x((i-1)*npj*npk+(j-1)*npk+k)+Lx/(npi-1)*(0.0*(rand()-0.5));
      y((i-1)*npj*npk+(j-1)*npk+k) = y((i-1)*npj*npk+(j-1)*npk+k)+Ly/(npj-1)*(0.0*(rand()-0.5));
      z((i-1)*npj*npk+(j-1)*npk+k) = z((i-1)*npj*npk+(j-1)*npk+k)+Lz/(npk-1)*(0.0*(rand()-0.5));
  end
  end
  end
  end
  end
  end
  
%   %Boundary conditions
%   if i==1
%       BSYMM = [BSYMM (i-1)*npj*npk+(j-1)*npk+k];
%   end
%   if i==npi
%       BSYMM = [BSYMM (i-1)*npj*npk+(j-1)*npk+k];
%   end
%   if j==1
%       BSYMM = [BSYMM (i-1)*npj*npk+(j-1)*npk+k];
%   end
%   if j==npj
%       BSYMM = [BSYMM (i-1)*npj*npk+(j-1)*npk+k];
%   end
%   if k==1
%       BSYMM = [BSYMM (i-1)*npj*npk+(j-1)*npk+k];
%   end
%   if k==npk
%       BSYMM = [BSYMM (i-1)*npj*npk+(j-1)*npk+k];
%   end
end
end
end

BSSONIC = unique(BSSONIC);
BSYMM   = unique(BSYMM);
BWALL   = unique(BWALL);

%Build the connectivity
E = zeros(ni*nj*nk,8);

for i=1:ni
for j=1:nj
for k=1:nk
  E((i-1)*nj*nk+(j-1)*nk+k,1) = (i-1)*npj*npk+(j-1)*npk+k;
  E((i-1)*nj*nk+(j-1)*nk+k,2) = (i-1)*npj*npk+(j-1)*npk+k+1;
  E((i-1)*nj*nk+(j-1)*nk+k,3) = (i-1)*npj*npk+(j-1+1)*npk+k;
  E((i-1)*nj*nk+(j-1)*nk+k,4) = (i-1)*npj*npk+(j-1+1)*npk+k+1;
  E((i-1)*nj*nk+(j-1)*nk+k,5) = (i-1+1)*npj*npk+(j-1)*npk+k;
  E((i-1)*nj*nk+(j-1)*nk+k,6) = (i-1+1)*npj*npk+(j-1)*npk+k+1;
  E((i-1)*nj*nk+(j-1)*nk+k,7) = (i-1+1)*npj*npk+(j-1+1)*npk+k;
  E((i-1)*nj*nk+(j-1)*nk+k,8) = (i-1+1)*npj*npk+(j-1+1)*npk+k+1;
end
end
end

[x y z];
E;

%BUILD PERIODIC BOUNDARIES
PER = [];

%k=1 and k=nk
for i=1:ni
for j=1:nj
  k  = 1;
  E1 = (i-1)*nj*nk+(j-1)*nk+k-1;
  k  = nk;
  E2 = (i-1)*nj*nk+(j-1)*nk+k-1;
  
  F1 = 4;
  F2 = 5;
  
  PER = [PER;E1 E2 F1 F2];
end
end

%j=1 and j=nj
for i=1:ni
for k=1:nk
  j  = 1;
  E1 = (i-1)*nj*nk+(j-1)*nk+k-1;
  j  = nj;
  E2 = (i-1)*nj*nk+(j-1)*nk+k-1;
  
  F1 = 2;
  F2 = 3;
  
  PER = [PER;E1 E2 F1 F2];
end
end

%i=1 and i=ni
for j=1:nj
for k=1:nk
  i  = 1;
  E1 = (i-1)*nj*nk+(j-1)*nk+k-1;
  i  = ni;
  E2 = (i-1)*nj*nk+(j-1)*nk+k-1;
  
  F1 = 0;
  F2 = 1;
  
  PER = [PER;E1 E2 F1 F2];
end
end
PER

plot3(x,y,z,'o')

%Flag == 0 METIS, else do simple cartesian partitioning
%Write the Connectivity info for METIS
fid = fopen('connectivity.dat', 'w');
fprintf(fid, '%d \n',  ni*nj*nk);
for i=1:ni*nj*nk
    fprintf(fid, '%d %d %d %d %d %d %d %d \n', E(i,1), E(i,2), E(i,3), E(i,4), E(i,5), E(i,6), E(i,7), E(i,8));
end
fclose(fid);

%Call METIS
if N>1
    ['mpmetis -contig connectivity.dat' num2str(N)]
    system(['mpmetis -contig connectivity.dat ' num2str(N)]); 

    horzcat('connectivity.dat.epart.',num2str(N))
    fid = fopen(horzcat('connectivity.dat.epart.',num2str(N)));
    P = fscanf(fid,'%i\n');
    fclose(fid)
else
    P = zeros(max(size(E)),1);
end

%write to file
fid = fopen(fname, 'w');
fprintf(fid, '%s \n', 'Number of grid points: ');
fprintf(fid, '%d \n',  npi*npj*npk);
fprintf(fid, '%s \n', 'Number of HEXAHEDRAL:  ');
fprintf(fid, '%d \n',  ni*nj*nk);
fprintf(fid, '%s \n', 'Nodes coordinates: ');
for i=1:npi*npj*npk
    fprintf(fid, '%12.12e %12.12e %12.12e \n', x(i), y(i), z(i));
end
fprintf(fid, '%s \n', 'Connectivity HEXAHEDRAL:  ');
for i=1:ni*nj*nk
    fprintf(fid, '%d %d %d %d %d %d %d %d %d %d\n', E(i,1), E(i,2), E(i,3), E(i,4), E(i,5), E(i,6), E(i,7), E(i,8), P(i), 0);
end
% fprintf(fid, '%s %d \n', 'SYMMETRY ',max(size(BSYMM)));
% for i=1:max(size(BSYMM))
%     fprintf(fid, '%d\n', BSYMM(i));
% end
fprintf(fid, '%s %d \n', 'PERIODIC ',max(size(PER)));
for i=1:max(size(PER))
    fprintf(fid, '%d %d %d %d\n',PER(i,1),PER(i,2),PER(i,3),PER(i,4));
end
fclose(fid);

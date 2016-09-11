function generator(ni,nj,N)

%Number of elements in the x, y, and z directions

Lx = 16;
Ly = 16;

fname = horzcat(horzcat(horzcat(horzcat(horzcat(num2str(ni),'x'),num2str(nj)),'_'),num2str(N)),'.msh');

npi = ni+1;
npj = nj+1;

x   = zeros(npi*npj,1);
y   = zeros(npi*npj,1);

BSSONIC = [];

%Build the list of nodes
for i=1:npi
for j=1:npj
  x((i-1)*npj+j) = Lx/(npi-1)*((i-1)+0.0*(rand()-0.5))-Lx/2;
  y((i-1)*npj+j) = Ly/(npj-1)*((j-1)+0.0*(rand()-0.5))-Ly/2;
  
  %Boundary conditions
  if i==1
      BSSONIC = [BSSONIC (i-1)*npj+j];
  end
  if i==npi
      BSSONIC = [BSSONIC (i-1)*npj+j];
  end
  if j==1
      BSSONIC = [BSSONIC (i-1)*npj+j];
  end
  if j==npj
      BSSONIC = [BSSONIC (i-1)*npj+j];
  end
end
end

BSSONIC = unique(BSSONIC);

%Build the connectivity
E = zeros(ni*nj,4);

for i=1:ni
for j=1:nj
  E((i-1)*nj+j,1) = (i-1+1)*npj+j;
  E((i-1)*nj+j,2) = (i-1+1)*npj+j+1;
  E((i-1)*nj+j,3) = (i-1  )*npj+j+1;
  E((i-1)*nj+j,4) = (i-1  )*npj+j;
end
end

%Write the Connectivity info for METIS
fid = fopen('connectivity.dat', 'w');
fprintf(fid, '%d \n',  ni*nj);
for i=1:ni*nj
    fprintf(fid, '%d %d %d %d \n', E(i,1), E(i,2), E(i,3), E(i,4));
end
fclose(fid);

%Call METIS
if N>1
    ['mpmetis connectivity.dat' num2str(N)]
    system(['mpmetis connectivity.dat ' num2str(N)]); 

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
fprintf(fid, '%d \n',  npi*npj);
fprintf(fid, '%s \n', 'Number of QUADS:  ');
fprintf(fid, '%d \n',  ni*nj);
fprintf(fid, '%s \n', 'Nodes coordinates: ');
for i=1:npi*npj
    fprintf(fid, '%12.12e %12.12e \n', x(i), y(i));
end
fprintf(fid, '%s \n', 'Connectivity QUAD:  ');
for i=1:ni*nj
    fprintf(fid, '%d %d %d %d %d %d \n', E(i,1), E(i,2), E(i,3), E(i,4), P(i), 0);
end
fprintf(fid, '%s %d \n', 'RIEMANN ',max(size(BSSONIC)));
for i=1:max(size(BSSONIC))
    fprintf(fid, '%d\n', BSSONIC(i));
end
fclose(fid);
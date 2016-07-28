function generator(ni,nj,N)

%Number of elements in the x, y, and z directions

Lx = 16;
Ly = 16;
Lz = 1;

nk = 1;

fname = horzcat(horzcat(horzcat(horzcat(horzcat(num2str(ni),'x'),num2str(nj),'x',num2str(nk)),'_'),num2str(N)),'.msh');

npi = ni+1;
npj = nj+1;
npk = nk+1;

x   = zeros(npi*npj*npk,1);
y   = zeros(npi*npj*npk,1);
z   = zeros(npi*npj*npk,1);

BSSONIC = [];
BSYMM   = [];

%Build the list of nodes
for i=1:npi
for j=1:npj
for k=1:npk
  x((i-1)*npj*npk+(j-1)*npk+k) = Lx/(npi-1)*((i-1)+0.0*(rand()-0.5))-Lx/2;
  y((i-1)*npj*npk+(j-1)*npk+k) = Ly/(npj-1)*((j-1)+0.0*(rand()-0.5))-Ly/2;
  z((i-1)*npj*npk+(j-1)*npk+k) = Lz/(npk-1)*((k-1)+0.0*(rand()-0.5))-Lz/2;
  
  %Boundary conditions
  if i==1
      BSSONIC = [BSSONIC (i-1)*npj*npk+(j-1)*npk+k];
  end
  if i==npi
      BSSONIC = [BSSONIC (i-1)*npj*npk+(j-1)*npk+k];
  end
  if j==1
      BSSONIC = [BSSONIC (i-1)*npj*npk+(j-1)*npk+k];
  end
  if j==npj
      BSSONIC = [BSSONIC (i-1)*npj*npk+(j-1)*npk+k];
  end
  if k==1
      BSYMM = [BSYMM (i-1)*npj*npk+(j-1)*npk+k];
  end
  if k==npk
      BSYMM = [BSYMM (i-1)*npj*npk+(j-1)*npk+k];
  end
end
end
end

BSSONIC = unique(BSSONIC);
BSYMM   = unique(BSYMM);

%Build the connectivity
E = zeros(ni*nj*nk,8);

% for i=1:ni
% for j=1:nj
% for k=1:nk
%   E((i-1)*nj*nk+(j-1)*nk+k,1) = (i-1)*npj*npk+(j-1)*npk+k;
%   E((i-1)*nj*nk+(j-1)*nk+k,2) = (i-1)*npj*npk+(j-1)*npk+k+1;
%   E((i-1)*nj*nk+(j-1)*nk+k,3) = (i-1)*npj*npk+(j-1+1)*npk+k;
%   E((i-1)*nj*nk+(j-1)*nk+k,4) = (i-1)*npj*npk+(j-1+1)*npk+k+1;
%   E((i-1)*nj*nk+(j-1)*nk+k,5) = (i-1+1)*npj*npk+(j-1)*npk+k;
%   E((i-1)*nj*nk+(j-1)*nk+k,6) = (i-1+1)*npj*npk+(j-1)*npk+k+1;
%   E((i-1)*nj*nk+(j-1)*nk+k,7) = (i-1+1)*npj*npk+(j-1+1)*npk+k;
%   E((i-1)*nj*nk+(j-1)*nk+k,8) = (i-1+1)*npj*npk+(j-1+1)*npk+k+1;
% end
% end
% end

for i=1:ni
for j=1:nj
for k=1:nk
  E((i-1)*nj*nk+(j-1)*nk+k,1) = (i-1  )*npj*npk+(j-1  )*npk+k+1;
  E((i-1)*nj*nk+(j-1)*nk+k,2) = (i-1+1)*npj*npk+(j-1  )*npk+k+1;
  E((i-1)*nj*nk+(j-1)*nk+k,3) = (i-1  )*npj*npk+(j-1+1)*npk+k+1;
  E((i-1)*nj*nk+(j-1)*nk+k,4) = (i-1+1)*npj*npk+(j-1+1)*npk+k+1;
  E((i-1)*nj*nk+(j-1)*nk+k,5) = (i-1  )*npj*npk+(j-1  )*npk+k;
  E((i-1)*nj*nk+(j-1)*nk+k,6) = (i-1+1)*npj*npk+(j-1  )*npk+k;
  E((i-1)*nj*nk+(j-1)*nk+k,7) = (i-1  )*npj*npk+(j-1+1)*npk+k;
  E((i-1)*nj*nk+(j-1)*nk+k,8) = (i-1+1)*npj*npk+(j-1+1)*npk+k;
end
end
end

% %BUILD PERIODIC BOUNDARIES
PER = [];

%Vertical Faces
for j=0:nj-1
    E1 = j;
    E2 = j+(nj-1)*ni;
    F1 = 4;
    F2 = 5;
    
    PER = [PER;E1 E2 F1 F2];
end

%Horizontal Faces
for j=0:nj-1
    E1 = (0   )+j*(ni);
    E2 = (ni-1)+j*(ni);
    F1 = 2;
    F2 = 3;
    
    PER = [PER;E1 E2 F1 F2];
end
PER

%Write the Connectivity info for METIS
fid = fopen('connectivity.dat', 'w');
fprintf(fid, '%d \n',  ni*nj*nk);
for i=1:ni*nj*nk
    fprintf(fid, '%d %d %d %d %d %d %d %d \n', E(i,1), E(i,2), E(i,3), E(i,4), E(i,5), E(i,6), E(i,7), E(i,8));
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

%RANDOM IMEX DISTRIBUTION
IMEX = zeros(max(size(E)),1);
% for i=1:max(size(E))
%     val = rand();
%     
%     if(val<0.1)
%         IMEX(i) = 1;
%     else
%         IMEX(i) = 0;
%     end
%     xc = (x(E(i,1))+x(E(i,2))+x(E(i,3))+x(E(i,4))+x(E(i,5))+x(E(i,6))+x(E(i,7))+x(E(i,8)))/8;
%     yc = (y(E(i,1))+y(E(i,2))+y(E(i,3))+y(E(i,4))+y(E(i,5))+y(E(i,6))+y(E(i,7))+y(E(i,8)))/8;
%     
%     r = sqrt(xc^2+yc^2)
%     IMEX(i) = 0;
%     if(r<1)
%         if(r>0.5)
%             IMEX(i) = 1;
%         end
%     end
% end

if(max(size(sort(unique(P))))~=N)
    disp('PARTITIONING ERROR!!!');
end

%L/R split
% for i=1:max(size(E))
%     if((x(E(i,1))+x(E(i,2))+x(E(i,3))+x(E(i,4))+x(E(i,5))+x(E(i,6))+x(E(i,7))+x(E(i,8)))/8<0)
%         P(i) = 1;
%     else
%         P(i) = 0;
%     end
% end
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
    fprintf(fid, '%d %d %d %d %d %d %d %d %d %d \n', E(i,1), E(i,2), E(i,3), E(i,4), E(i,5), E(i,6), E(i,7), E(i,8), P(i), IMEX(i));
end
fprintf(fid, '%s %d \n', 'SYMMETRY ',max(size(BSYMM)));
for i=1:max(size(BSYMM))
    fprintf(fid, '%d\n', BSYMM(i));
end
fprintf(fid, '%s %d \n', 'RIEMANN ',max(size(BSSONIC)));
for i=1:max(size(BSSONIC))
    fprintf(fid, '%d\n', BSSONIC(i));
end
% fprintf(fid, '%s %d \n', 'PERIODIC ',max(size(PER)));
% for i=1:max(size(PER))
%     fprintf(fid, '%d %d %d %d\n',PER(i,1),PER(i,2),PER(i,3),PER(i,4));
% end
fclose(fid);

plot3(x,y,z,'k.')
hold on
for i=1:max(size(E))
    if(IMEX(i)) == 1
        xc = (x(E(i,1))+x(E(i,2))+x(E(i,3))+x(E(i,4))+x(E(i,5))+x(E(i,6))+x(E(i,7))+x(E(i,8)))/8;
        yc = (y(E(i,1))+y(E(i,2))+y(E(i,3))+y(E(i,4))+y(E(i,5))+y(E(i,6))+y(E(i,7))+y(E(i,8)))/8;
        plot(xc,yc,'ro')
    end
    if(IMEX(i)) == 0
        xc = (x(E(i,1))+x(E(i,2))+x(E(i,3))+x(E(i,4))+x(E(i,5))+x(E(i,6))+x(E(i,7))+x(E(i,8)))/8;
        yc = (y(E(i,1))+y(E(i,2))+y(E(i,3))+y(E(i,4))+y(E(i,5))+y(E(i,6))+y(E(i,7))+y(E(i,8)))/8;
        plot(xc,yc,'bo')
    end
end

PWX = [];
PWY = [];
PWZ = [];
for i=1:max(size(BSSONIC))
    PWX = [PWX x(BSSONIC(i))];
    PWY = [PWY y(BSSONIC(i))];
    PWZ = [PWZ z(BSSONIC(i))];
end
plot3(PWX,PWY,PWZ,'or')
hold off
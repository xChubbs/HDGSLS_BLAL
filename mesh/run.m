close all
% clear
% clc

malla=load_gmsh('cubeD.msh');

Elements=malla.TETS(:,1:end-1); %numero de test via nodo
nel=malla.nbTets;
nver=malla.nbNod;    %numero de nodos en la malla
Coordinates=malla.POS(:,1:end);   %Coordenadas de cada nodo

Dirichlet=find(malla.TRIANGLES(:,4)==50);
NodosDirichlet=unique(malla.TRIANGLES(Dirichlet,1:3));

Neumann=find(malla.TRIANGLES(:,4)==60);
NodosNeumann=malla.TRIANGLES(Neumann,1:3);

% Organicemos los resultados que tenemos en el momento en una estructura
% inicial que nos genere un orden similar para el procesamiento posterior.

% -------------- Procesamiento de Estructura aumentada ---------------- %

% Organizaremos como Sayas, la estructura en un valor T
T.coordinates = Coordinates;
T.elements    = Elements;

% Processing to get the dirichlet and Neuman surfaces
% - Dirichlet : identifier (50)
dirichlet_mask = find(malla.TRIANGLES(:, 4) == 50);
dirichlet_surfaces = malla.TRIANGLES(dirichlet_mask, 1:3);

% - Neumann : identifier (60)
neumann_mask = malla.TRIANGLES(:, 4) == 60;
neumann_surfaces = malla.TRIANGLES(neumann_mask, 1:3);

% Addition to structure
T.dirichlet = dirichlet_surfaces;
T.neumann   = neumann_surfaces;

% Construction of a list of all faces
shape=[1 3 2;1 2 4;1 4 3;2 3 4];
nelts=size(T.elements,1);
faces=zeros(4*nelts,3);
 for k=1:nelts
     nodes=T.elements(k,:); %1x4
     faces(4*(k-1)+(1:4),:)=nodes(shape); % 1x3
 end

% Processing of faces
copyoffaces=faces;  % 4 nelts x 3 (in local positive orientation)
faces=sort(faces,2);

% Storage of unique values for faces (No repetition)
[allfaces,~,j]=unique(faces,'rows');

% Difference and Intersection of Boundary Faces
bdfaces=sort([T.dirichlet;T.neumann],2);
[intfaces,i]=setdiff(allfaces,bdfaces,'rows');
[bdfaces,ii,jj]=intersect(allfaces,bdfaces,'rows');

nintfaces=size(intfaces,1);
ndirfaces=size(T.dirichlet,1);
nneufaces=size(T.neumann,1);
nbdfaces =size(bdfaces,1);
nfaces   =nintfaces+nbdfaces;

% Interior Faces ∩ Neumann Faces ∩ Dirichlet Faces = Total Faces
T.faces=[intfaces zeros(nintfaces,1);...
         T.dirichlet ones(ndirfaces,1);...
         T.neumann 2*ones(nneufaces,1)];
T.dirfaces=(nintfaces+1):(nintfaces+ndirfaces);
T.neufaces=(nintfaces+ndirfaces+1):(nintfaces+ndirfaces+nneufaces);

% Backward referencing to construct T.facebyele
u=zeros(nfaces,1);
v=nintfaces+1:nintfaces+nbdfaces;
u(i)=1:nintfaces;
u(ii)=v(jj);
j=u(j);                % pointer from T.faces to the copyoffaces
faces=T.faces(j,1:3);  % 4 nelts x 3, with global numbering of nodes
j=reshape(j,[4 nelts]);
T.facebyele=j';

% Matrix with orientations
A=T.facebyele';
faces=T.faces(A(:),1:3);
t=sum(faces==copyoffaces,2)==ones(4*nelts,1);
t=1-2*t;t=reshape(t,[4,nelts]);
T.orientation=t';

% Matrix with permutation order
eq = @(u,v) sum(u==v,2)==3; % checks what rows are equal
rot=[1 2 3;...              % permutations 
     1 3 2;...
     3 1 2;...
     3 2 1;...
     2 3 1;...
     2 1 3];
pattern=[1 2 3;...    % (s,t,0)
         1 2 4;...    % (s,0,t)
         1 3 4;...    % (0,s,t)
         4 2 3];      % (s,t,1-s-t)
orient=zeros(nelts,4);

for f=1:4    % counter over faces
    faceGlobal=T.faces(T.facebyele(:,f),1:3);
    faceLocal =T.elements(:,pattern(f,:));
    for j=1:6
        orient(:,f)=orient(:,f)+j*eq(faceGlobal,faceLocal(:,rot(j,:)));
    end
end
T.perm=orient;

% Volumes 
T.volume=(1/6)*...
 ((T.coordinates(T.elements(:,2),1)-T.coordinates(T.elements(:,1),1)).*...
     (T.coordinates(T.elements(:,3),2)-T.coordinates(T.elements(:,1),2)).*...
     (T.coordinates(T.elements(:,4),3)-T.coordinates(T.elements(:,1),3))-...
  (T.coordinates(T.elements(:,2),1)-T.coordinates(T.elements(:,1),1)).*...
     (T.coordinates(T.elements(:,3),3)-T.coordinates(T.elements(:,1),3)).*...
     (T.coordinates(T.elements(:,4),2)-T.coordinates(T.elements(:,1),2))-...
  (T.coordinates(T.elements(:,2),2)-T.coordinates(T.elements(:,1),2)).*...
     (T.coordinates(T.elements(:,3),1)-T.coordinates(T.elements(:,1),1)).*...
     (T.coordinates(T.elements(:,4),3)-T.coordinates(T.elements(:,1),3))+...
  (T.coordinates(T.elements(:,2),2)-T.coordinates(T.elements(:,1),2)).*...
     (T.coordinates(T.elements(:,3),3)-T.coordinates(T.elements(:,1),3)).*...
     (T.coordinates(T.elements(:,4),1)-T.coordinates(T.elements(:,1),1))+...
  (T.coordinates(T.elements(:,2),3)-T.coordinates(T.elements(:,1),3)).*...
     (T.coordinates(T.elements(:,3),1)-T.coordinates(T.elements(:,1),1)).*...
     (T.coordinates(T.elements(:,4),2)-T.coordinates(T.elements(:,1),2))-...
  (T.coordinates(T.elements(:,2),3)-T.coordinates(T.elements(:,1),3)).*...
     (T.coordinates(T.elements(:,3),2)-T.coordinates(T.elements(:,1),2)).*...
     (T.coordinates(T.elements(:,4),1)-T.coordinates(T.elements(:,1),1)));

% Areas definition
x_coordinate = ((T.coordinates(T.faces(:,2),2)-T.coordinates(T.faces(:,1),2)).*...
  (T.coordinates(T.faces(:,3),3)-T.coordinates(T.faces(:,1),3))-...
  (T.coordinates(T.faces(:,2),3)-T.coordinates(T.faces(:,1),3)).*...
  (T.coordinates(T.faces(:,3),2)-T.coordinates(T.faces(:,1),2))).^2;

y_coordinate = ((T.coordinates(T.faces(:,2),1)-T.coordinates(T.faces(:,1),1)).*...
  (T.coordinates(T.faces(:,3),3)-T.coordinates(T.faces(:,1),3))-...
  (T.coordinates(T.faces(:,2),3)-T.coordinates(T.faces(:,1),3)).*...
  (T.coordinates(T.faces(:,3),1)-T.coordinates(T.faces(:,1),1))).^2;

z_coordinate = ((T.coordinates(T.faces(:,2),1)-T.coordinates(T.faces(:,1),1)).*...
  (T.coordinates(T.faces(:,3),2)-T.coordinates(T.faces(:,1),2))-...
  (T.coordinates(T.faces(:,2),2)-T.coordinates(T.faces(:,1),2)).*...
  (T.coordinates(T.faces(:,3),1)-T.coordinates(T.faces(:,1),1))).^2;

T.area=(1/2)*sqrt( x_coordinate + y_coordinate + z_coordinate );
 

% Normals to the faces 
for f=1:4
    oneface=T.faces(T.facebyele(:,f),1:3);
    x=T.coordinates(oneface(:),1);
    y=T.coordinates(oneface(:),2);
    z=T.coordinates(oneface(:),3);
%
    x12=x(nelts+1:2*nelts)-x(1:nelts); %x_2-x_1
    x13=x(2*nelts+1:end)-x(1:nelts);   %x_3-x_1
    y12=y(nelts+1:2*nelts)-y(1:nelts); %y_2-y_1
    y13=y(2*nelts+1:end)-y(1:nelts);   %y_3-y_1
    z12=z(nelts+1:2*nelts)-z(1:nelts); %z_2-z_1
    z13=z(2*nelts+1:end)-z(1:nelts);   %z_3-z_1
 %   
    normals=(1/2)*[y12.*z13-y13.*z12,...
                   z12.*x13-z13.*x12,...
                   x12.*y13-x13.*y12];
    normals=bsxfun(@times,normals,T.orientation(:,f)); % Give the normals correctly orientated
%
    T.normals(:,(f-1)*3+(1:3))=normals;
end

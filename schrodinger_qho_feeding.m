 % ------------------------------------------------------------------------
%   @ autor : alujan
% ------------------------------------------------------------------------
%   @ title : HDG Scheme to eigenvalues of poisson equation
% ------------------------------------------------------------------------
%   @ about : Resource of bachelor's thesis
% ------------------------------------------------------------------------
%   @ credits : Universidad Nacional de Colombia, sede medellín
%               blopezr - Thesis advisor
% ------------------------------------------------------------------------

% [0] Preparation of files
% Clear of previous definitions
clear; format short; syms x y z real

disp('>> Lectura de directorios')

% Let's add to the path the src and drivers: Original Article
addpath('./src/')

% Let's add the path of the gmsh structure: Load Mesh structure
addpath('./mesh/')

disp('>> Lectura de mallas')

% Importing of a singular structure for first case scenario: Unitary cube
listT = {};     % Storage unit
for i=1:6
    % Addition of element
    listT{i} = load_mesh('L10-cubeD' + string(i) + '.msh');
end

disp('>> Última malla considerada:')
disp(listT{end})

disp('>> Definición de solución teórica')

% Let's define a simple value for U: harmonic vector
% - Characteristics of the oscillator particle
m = 1; k = 2; omega = sqrt(k / m); hbar = 1; x0 = 20/2;

% - Definition per coordinate of the oscillator
n = 2; Ar = 1 / sqrt(2^n * prod(1:n)) * (m*omega/(pi*hbar))^0.25;
Ur(x) = Ar * exp(-m*omega*(x-x0)^2/(2*hbar)) * ...
                          hermiteH(n, sqrt(m*omega/hbar)*(x-10));

% - Caracteristics of the function in space
U(x, y, z) = Ur(x)*Ur(y)*Ur(z); vpa(U, 2)

% Partial derivarives of the component
Ux = diff(U, x); Uxx = diff(U, x, 2);
Uy = diff(U, y); Uyy = diff(U, y, 2);
Uz = diff(U, z); Uzz = diff(U, z, 2);

Q = -[Ux; Uy; Uz]; vpa(Q, 2)   % Gradient definition
% Definition of potential taking action
F(x, y, z) = - (Uxx + Uyy +Uzz); vpa(F, 2)

% Processing of functions
u = matlabFunction(U, Vars={x, y, z});

qx = matlabFunction(-Ux, Vars={x, y, z});
qy = matlabFunction(-Uy, Vars={x, y, z});
qz = matlabFunction(-Uz, Vars={x, y, z});

% It's possible to change the values below to alter the polynomial and
% cuadrature precision for more general viewing
TablesQuadForm3d        % Importing of cuadrature for thetrahedron
TablesQuadForm          % Importing of cuadrature for triangles

% Norms of unknowns for relative error
Tmax=listT{end};                 % Listing of T elements

% [1] \forall k \in [1, .., 4]
for k=1:4

    disp('>> Recorrido k = ' + string(k))

    % \forall T \in listT
    for i=1:size(listT, 2)

        disp('>> Recorrido de malla T = ' + string(i))

        % Definition of current T
        T = listT{i};                    % Current mesh
    
        % Error elements definition for display handle
        ErrorU = []; ErrorUhat = []; ErrorLambda = [];

        disp('>> Definición de cuadratura y métricas de error')
        
        % Switch case for the selection of formulas: cuadratures
        switch k
            case 0 
                formulas={tetra3,tetra1,matrix0,matrix4};
            case 1
                formulas={tetra5,tetra3,matrix4,matrix9};
            case 2
                formulas={tetra7,tetra5,matrix9,matrix11};
            case 3
                formulas={tetra9,tetra7,matrix11,matrix14};
            case 4
                formulas={tetra9,tetra9,matrix14,matrix16};
        end
        
        % Reduction of weights: Sum of weights it's defined to be equal to 2, for
        % normalization it's reduced to 1.
        formulas{3}(:,4)=formulas{3}(:,4)/2;
        formulas{4}(:,4)=formulas{4}(:,4)/2;

        % [3] Definitions associated with DG
        disp('>> Definiciones de incognitas y enumeraciones de DG')
        
        % -------------------------------- HDG 3D main function --------------------------------------- %
        
        % Definitions asociated with the tetrahedrization
        d2 = nchoosek(k+2,2);    
        d3 = nchoosek(k+3,3); 
        
        block3 = @(x) (1+(x-1)*d3):(x*d3);
        
        Nelts  = size(T.elements,1);
        Nfaces = size(T.faces,1);
        Ndir   = size(T.dirichlet,1);
        Nneu   = size(T.neumann,1);
        
        % 4th intervention
        % - Faces
        face = T.facebyele';                  % 4 x Nelts
        face = (face(:)-1)*d2;                % First degree of freedom of each face by element                                               
        face = bsxfun(@plus,face,1:d2);       % 4*Nelts x d2 (d.o.f. for each face)                                              
        face = reshape(face',4*d2,Nelts);     % d.o.f. for the 4 faces of each element
        
        [J,I] = meshgrid(1:4*d2);
        
        R_d2 = face(I(:),:); R_d2=reshape(R_d2,4*d2,4*d2,Nelts);   % Mesh grid definition inside T_i element "Row-Axis"
        C_d2 = face(J(:),:); C_d2=reshape(C_d2,4*d2,4*d2,Nelts);   % Mesh grid definition inside T_i element "Column-Axis"
        
        % - Volumes
        element = (1:Nelts);                          % Nelts
        element = (element(:)-1)*4*d3;                % First degree of freedom of each face by element                                               
        element = bsxfun(@plus,element,1:4*d3);       % 4*Nelts x d2 (d.o.f. for each face)                                              
        element = reshape(element',4*d3,Nelts);       % d.o.f. for the 4 faces of each element
        
        [J,I] = meshgrid(1:4*d3); 
        
        R_d3 = element(I(:),:); R_d3 = reshape(R_d3,4*d3,4*d3,Nelts);   % Mesh grid definition inside T_i element "Row-Axis"
        C_d3 = element(J(:),:); C_d3 = reshape(C_d3,4*d3,4*d3,Nelts);   % Mesh grid definition inside T_i element "Column-Axis"
        
        % R_ij^K d.o.f. for local (i,j) d.o.f. in element K ; R_ij^K=C_ji^K
        RowsRHS_d2 = reshape(   face,4*d2*Nelts,1);            % Definition of solution shape   
        RowsRHS_d3 = reshape(element,4*d3*Nelts,1);              % Definition of solution shape   
        
        % 4th intervention
        % Dirichlet faces: defintion
        dirfaces = (T.dirfaces(:)-1)*d2;          % First degree of freedom of each face by element
        dirfaces = bsxfun(@plus,dirfaces,1:d2);   % Bitwise sum of dirfaces and 1:d2 > dirfaces + 1:d2
        dirfaces = reshape(dirfaces',d2*Ndir,1);  % Reshape of dirfaces
        
        % Reduced version of free
        free = (1: d2*Nfaces);
        
        % Cleaning of dirfaces: Positions not needed
        free(dirfaces) = [];      % Empty of dirfaces positions on free
        
        % Neumann faces: defintion
        neufaces = (T.neufaces(:)-1)*d2;          % First degree of freedom of each face by element
        neufaces = bsxfun(@plus,neufaces,1:d2);   % Bitwise sum of neufaces and 1:d2 > neufaces + 1:d2
        neufaces = reshape(neufaces',d2*Nneu,1);
        
        % Reordering of values find indexes
        indexes = arrayfun(@(i) find(T.facebyele == T.dirfaces(i)), 1:size(T.dirfaces, 2));
        
        % - Map of values indexes
        map = T.orientation(indexes) == -1; T.faces(T.dirfaces(map), [2 3]) = T.faces(T.dirfaces(map), [3 2]);
        
        % Correction of orientation
        T.orientation(indexes) = 1;

        % ----------------------------- Local solvers definition function ------------------------------ %
        % Definition of the action of tau: Identity
        tau = 1*ones(4,Nelts);

        % Computation of the differentials form Matrix: Return in the already desired form
        disp('>> Computo de matrices asociadas al sistema')
        disp('   - Matriz de convección')
        % - Computation of Convection Integrals
        [CMx, CMy, CMz] = ConvMatrix(T,k,formulas{1});
        
        disp('   - Matriz de masa: Mk y Mv')
        % - Computation of the Mass Integrals
        Mass = MassMatrix(T,{@(x, y, z) 0*x.*y.*z + 1, ...
                     @(x, y, z) 0.5*m*omega^2*((x - x0).^2 + (y - x0).^2 + (z - x0).^2)}, k,formulas{1});

        % - ALR: Reduce this computation into one single line
        Mi=Mass{1}; Mv=Mass{2};
        
        disp('   - Matriz productos en caras')
        % - Computation of Faces Integrals
        [tauPP,tauDP,nxDP,nyDP,nzDP,tauDD]=matricesFace(T,tau,k,formulas{3});
        
        % Condensation of matrices known
        nDP = [nxDP, nyDP, nzDP]; divPP = [CMx, CMy, CMz];
        
        % - Matrices defintions: A solvers local
        A1=zeros(4*d3,4*d3,Nelts);    % Solver storage reservation
        A2=zeros(4*d3,4*d2,Nelts);    % Solver storage reservation
        A3=zeros(4*d3,4*d3,Nelts);    % Solver storage reservation
        A4=zeros(4*d2,4*d3,Nelts);    % Solver storage reservation
        
        % Matrices related to solvers
        % - Definition of the Mass matrix asociated with a_1: 6*d3 x 6*d3
        O = zeros(d3,d3,Nelts);   % Big O cero matrix: Shapes juntions
        
        disp('   - Definición de matrices aumentadas')
        % - Block definition of matrix
        Mk = [Mi ,O  ,O   ; ...
              O  ,Mi ,O   ; ...
              O  ,O  ,Mi ];
        
        A1 = [       Mk, -permute(divPP, [2 1 3]); ...
              0.5*divPP,  tauPP + Mv            ];
        
        % - Definition of  the Mass matrix asociated with a_2: 6*d3 x 2*(4*d2)
        %   Pendiente la revisión de errores propuesta
        A2 = [-permute(nDP  , [2 1 3])  ; ...
               permute(tauDP, [2 1 3]) ];
        
        O = zeros(3*d3,3*d3,Nelts);             % Big O cero matrix: Shapes juntions
        Ov = zeros(3*d3,d3,Nelts);
        
        A3_w = [                  Mk,  Ov;....
                permute(Ov, [2 1 3]), Mi];
        
        A3 = [                   O,  Ov;....
              permute(Ov, [2 1 3]), Mi];
        
        O = zeros(d3,d3,Nelts); Ov = [O, O, O]; % Big O cero matrix: Shapes juntions
        
        A4 = [Mk , permute(Ov, [2, 1, 3]) ; ...
              Ov , tauPP + Mv            ];
        
        O = zeros(4*d2,d3,Nelts);
        
        A5 = [O, O, O, tauDP];
        
        % Matrices defintions: C solvers Flux
        CM = zeros(4*d2,4*d2,Nelts);         % Solver storage reservation
        CK = zeros(4*d2,4*d2,Nelts);         % Solver storage reservation
        
        % Operator definition
        Q_w_U_w = zeros(4*d3, 4*d3, Nelts);  % Solver operator storage reservation
        
        Q_w     = zeros(4*d3, 4*d3, Nelts);  % Solver Q_w storage reservation
        U_w     = zeros(4*d3, 4*d3, Nelts);  % Solver U_w storage reservation
        
        % Parallel creation of flux operators
        % Revisión del sistema final a generar
        % - Element separation block
        block_3d3 = [block3(1) block3(2) block3(3)];

        disp('   - Generación de sistema')
        parfor i=1:Nelts            
            CM(:,:,i) = A2(:, :, i)' * (A1(:, :, i)' \ A4(:, :, i)) * (A1(:, :, i) \ A2(:, :, i)) - ...
                        A2(:, :, i)' * (A1(:, :, i)' \ A5(:, :, i)') - ...
                        A5(:, :, i)  * (A1(:, :, i)  \ A2(:, :, i)) + tauDD(:, :, i);
            CK(:,:,i) = A2(:, :, i)' * (A1(:, :, i)' \ A3(:, :, i)) * (A1(:, :, i) \ A2(:, :, i));
        
            Q_w_U_w(:, :, i) = A1(:, :, i) \ A3_w(:, :, i);
        end
        
        % Separation of operators 
        U_w(block3(4),block3(4), :) = Q_w_U_w(block3(4), block3(4), :);
        Q_w(block_3d3,block_3d3, :) = Q_w_U_w(block_3d3, block_3d3, :);

        disp('   - Sistema completo de caras')
        % Recovery of all elements
        M = sparse(R_d2(:),C_d2(:), CM(:)); K = sparse(R_d2(:), C_d2(:),CK(:));
        
        % Boundary conditions evaluation
        % The resulting terms are:
        %  - uhartD : Reference term evaluation on the Dirichlet faces
        %  - qhatN  : Reference term evaluation on the Neumann faces
        [etahatD, qhatN] = BC3d(@(x, y, z) 0*x,@(x, y, z) 0*x,@(x, y, z) 0*x,@(x, y, z) 0*x,T,k,formulas{3});
        
        % Dirichlet BC: Storage
        etahatD=reshape(etahatD,d2*Ndir,1); % uhatD stored as a vector: (2*d2*Ndir) x 1 : todo > Adjust space
        Etahatv=zeros(d2*Nfaces, 1);        % Storage reservation
        Etahatv(dirfaces)=etahatD;          % Uhat stored as a vector: d2*Nfaces
        Etahatv(free)    = 1;
        
        % Export of the system
        system={M, K, free, dirfaces};
        solvers={A1, A2, A3, A4, A5};
        Uh=[]; Qxh=[]; Qyh=[]; Qzh=[]; Uhat=[];
        
        disp('   - Solución de problema de autovalores')
        % RHS currently has the information of the frontier conditions 
        % and the information of inner skeleton: solution of uhat
        [Etav_free, lambda] = eigs(M(free, free), K(free, free), 1, (hbar*omega^2*(0 + 3/2)));

        % ---------------------------------- Local solutions definitions ---------------------------------- %

        % Insertion of BC
        Etav = zeros(d2*Nfaces, 1);
        Etav(dirfaces, :) = etahatD;
        Etav(free, :)     = Etav_free - M(free,dirfaces)*Etav(dirfaces, :);
        
        % Reorganization of values
        Eta=zeros(d2, Nfaces);
        
        Uh=zeros(d3,Nelts);                                 % Reservation of storage
        
        Qxh=zeros(d3,Nelts);                                % Reservation of storage
        Qyh=zeros(d3,Nelts);                                % Reservation of storage
        Qzh=zeros(d3,Nelts);                                % Reservation of storage
        
        sol = zeros(4*d3, 1);                                   % Temporal reservation of space
        
        % Parallel solution of system: Newton algorithm
        parfor N = 1:1
        
            % Assignation of variables
            Etav_n = Etav(:, N);
        
            % Energies found
            lambda_0 = lambda(N); lambda_n = lambda_0; 
        
            % Tolerance till convergence
            tolerance = 1E-5; delta_refinement = 1E-3;
            
            scale = 1; delta_lambda = 1; 
        
            % Update of Eta
            % - Reordering of faces
            Eta_n = reshape(Etav_n, [d2, Nfaces]);
        
            % - Recover of values per elements
            faces   = T.facebyele'; faces=faces(:);                         % Recovery of faces information
            eta_aux = reshape(Eta_n(:, faces, :),[4*d2, Nelts]);            % Reshape on the faces of Uhat
        
            % Update of energies
            lambda(N) = lambda_n;
        
            % Exit variable of parfor
            Eta(:, :, N) = Eta_n;
        
            % Reconstruction of solution per energy
            for i = 1:Nelts
        
                % Reconstruction of solution
                sol = A1(:,:,i)\A2(:,:,i) * eta_aux(:,i);
        
                % Fill of spaces
                temporal = (eye(4*d3, 4*d3) - lambda_n * U_w(:, :, i)) \ sol(:);
        
                Uh(:,i,N) = temporal(block3(4));
        
                Qxh(:,i,N) = sol(block3(1)) + lambda_n * Q_w(block3(1), block3(1), i) * Uh(:,i,N);
                Qyh(:,i,N) = sol(block3(2)) + lambda_n * Q_w(block3(2), block3(2), i) * Uh(:,i,N);
                Qzh(:,i,N) = sol(block3(3)) + lambda_n * Q_w(block3(3), block3(3), i) * Uh(:,i,N);
            end
        
        end

        % Generation of mat structure
        name_mat = 'qho_var_k_' + string(k) + '.mat';
    
        qho_var_k_struct.Uh     = Uh;
        qho_var_k_struct.Qh     = cat(3, Qxh(:, :, 1), Qyh(:, :, 1), Qzh(:, :, 1));
        qho_var_k_struct.lambda = lambda;
    
        save(name_mat, '-struct', 'qho_var_k_struct')

        % Addition of errors
        disp('   - Adición de errores')

        % Norms definitions
        normU = errorElem(Tmax,u,zeros(d3,size(Tmax.elements, 1)),k,formulas{1});
        
        normQ = errorElem(Tmax,qx,zeros(d3,size(Tmax.elements, 1)),k,formulas{1})...
               + errorElem(Tmax,qy,zeros(d3,size(Tmax.elements, 1)),k,formulas{1})...
               + errorElem(Tmax,qz,zeros(d3,size(Tmax.elements, 1)),k,formulas{1});
        
        normEta = errorFaces(Tmax,u,zeros(d2,size(Tmax.faces, 1)),k,formulas{4});
        
        % Errors with exact solution
        error_eta = errorFaces(T,  u, Eta(:, :, 1), k, formulas{4});
        
        error_q   = errorElem (T, qx, Qxh(:, :, 1), k, formulas{1})...
                  + errorElem (T, qy, Qyh(:, :, 1), k, formulas{1})...
                  + errorElem (T, qz, Qzh(:, :, 1), k, formulas{1});
        
        error_u   = errorElem (T,  u, Uh(:, :, 1), k, formulas{1});
        
        error_eta_ = errorFaces(T,  u, -Eta(:, :, 1), k, formulas{4});
        error_q_   = errorElem (T, qx, -Qxh(:, :, 1), k, formulas{1})...
                   + errorElem (T, qy, -Qyh(:, :, 1), k, formulas{1})...
                   + errorElem (T, qz, -Qzh(:, :, 1), k, formulas{1});
        error_u_   = errorElem (T,  u, - Uh(:, :, 1), k, formulas{1});
        
        error_eta = min([error_eta, error_eta_]);
        error_u   = min([  error_u,   error_u_]);
        error_q   = min([  error_q,   error_q_]);
        
        ErrorU=[ErrorU error_u/normU];
        ErrorUhat=[ErrorUhat error_eta/normEta];
        ErrorLambda = [ErrorLambda abs((hbar*omega^2*(0 + 3/2)) - lambda) / ((hbar*omega^2*(0 + 3/2)))];

        disp('   - El error reportado en lambda es:')
        disp(ErrorLambda(end))
    end

    disp('   - Creación de tablas de convergencia')

    % Rate of sucesion: Convergence verification
    rateU      = log2(ErrorU(1:end-1)./ErrorU(2:end));
    rateEta    = log2(ErrorUhat(1:end-1)./ErrorUhat(2:end));
    rateLambda = log2(ErrorLambda(1:end-1)./ErrorLambda(2:end));

    % Display sequence of values
    format bank
    disp('Rate |U-Uh| Rate |U-Uhat|_h Rate |\lambda - \lambda_h|');
    disp([rateU' rateEta' rateLambda']);
    format shortE
    disp('Error |U-Uh| Error |U-Uhat|_h Error |\lambda - \lambda_h|');
    disp([ErrorU' ErrorUhat' ErrorLambda']);

    % Saving of values
    disp('   - Guardado de tablas de ratas y errores')

    % Generation of mat structure
    name_mat = 'qho_rat_err_k_' + string(k) + '.mat';

    rat_err_k_struct.rateU = rateU;
    rat_err_k_struct.rateEta = rateEta;
    rat_err_k_struct.rateLambda = rateLambda;

    rat_err_k_struct.ErrorU = ErrorU;
    rat_err_k_struct.ErrorUhat = ErrorUhat;
    rat_err_k_struct.ErrorLambda = ErrorLambda;

    save(name_mat, '-struct', 'rat_err_k_struct')
    
end

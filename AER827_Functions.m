%Matrix Library for AER827 Composite Materials
%Written by Jim Sadden, Feb 2019
format shortG;

%Functions:
% getT3: inputs: theta of ply in degrees
%        outputs: 3x3 Transformation matrix T

% getQ3: inputs: E1, E2, v12, v21, G12 
%        outputs: 3x3 stiffness matrix Q

%getQbar3: inputs: 3x3 stiffness matrix Q, theta in degrees
%          outputs: 3x3 Qbar matrix

%stressLocalToGlobal3: inputs: 3x1 local stress vector, 3x3 transformation
                               %matrix T
%                      outputs: 3x1 global stress vector

%stressGlobalToLocal3: inputs: 3x1 global stress vector, 3x3 transformation
                               %matrix T
%                      outputs: 3x1 local stress vector

%getS3: inputs: E1, E2, v12, v21, G12
%       outputs: 3x3 compliance matrix S


%getA3: inputs: matrix of k Qbar matrices Q, organize as Q = [Qbar1, Qbar2, ...Qbark]
              % (k+1)x1 z vector of locations of plies, organize as z = [z0; z1; ...zk]
              % z vector must be one longer than number of Qbar matrices in
              % Q
%       outputs: Extended stiffness matrix A

%getB3: inputs: matrix of k Qbar matrices Q, organize as Q = [Qbar1, Qbar2, ...Qbark]
              % (k+1)x1 z vector of locations of plies, organize as z = [z0; z1; ...zk]
              % z vector must be one longer than number of Qbar matrices in
              % Q
%       outputs: Coupling stiffness matrix B

%getD3: inputs: matrix of k Qbar matrices Q, organize as Q = [Qbar1, Qbar2, ...Qbark]
              % (k+1)x1 z vector of locations of plies, organize as z = [z0; z1; ...zk]
              % z vector must be one longer than number of Qbar matrices in
              % Q
%       outputs: Bending stiffness matrix D

%getMidPlaneStrains: inputs: A, B, D, matrices and 6x1 Force/Moment vector
%                    outputs: 6x1 mid plane strain and curvature vector

%getForcesMoments: inputs: A, B, D, matrices and 6x1 mid-plane
                          %Strain/Curvature vector 
%                  outputs: 6x1 forces/moments vector

%getGlobalStrains: inputs: 6x1 mid plane strain/curvature vector, z
                          %location of ply in focus
%                  outputs: 3x1 global strain vector for ply in focus

%getStressFromStrain: inputs: 3x3 stiffness matrix Q or Qbar and 3x1 strain
                              %vector
%                     outputs: 3x1 stress vector

%getStrainFromStress: inputs: 3x3 comliance matrix S (inverse of Q) and 3x1
                              %stress vector
%                     outputs: 3x1 strain vector

%getQorS: inputs: Q or S to be inverted
%         outputs: If Q input, returns S, if S input, returns Q

%getLocalStrainFromGlobal: inputs: 3x1 global strain vector, transformation
                                  %matrix T
%                          outputs: 3x1 local strain vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%From properties to local strains and everything in between
%
%givePropsToTheLocals: inputs: E1, E2, v12, v21, G12, row vector(1xk) of k ply thetas,
                              %column vector((k+1)x1) of k+1 ply z locations, 6x1 forces/moments vector
%          displays/outputs: --0 degree Q matrix
%                            --Q bar matrices --> organized as 3x3s in one
%                              large 3x3k matrix
%                            --A, B, and D matrices
%                            --6x1 Midplane Strains/Curvatures vector
%                            --3x(k+1) Global strain vectors organized in
%                              columns
%`                           --3x(k*2) Global Stress vectors organized in
%                              columns
%                            --3x(k*2) Local Stress vectors organized in
%                              columns
%                            --3x(k*2) Local Strain vectors organized in
%                              columns
%
%NOTE: Also outputs ply number, ply angle, and position of reading. Each
%      should fall above their respective column of stresses/strains,
%      except for global strains, which only show the position


%%%%%%%%%%%%
%DEFINE INPUTS
%%%%%%%%%%%%
E1 = 1.81e11;
E2 = 1.03e10;
v12 = 0.28;
v21 = 0.01593;
G12 = 7.17e9;
theta = [-10, 30, -30, 45];
z = [-0.010; -0.005; 0.000; 0.005; 0.010];
NM = [10000;0;0;0;0;0];


%From properties to local strains and everything in between
givePropsToTheLocals(E1, E2, v12, v21, G12, theta, z, NM)





function givePropsToTheLocals(E1, E2, v12, v21, G12, thetas, z, NM)
    if (size(NM, 2) > 1 || size(NM, 1) > 6) 
        disp('Error: givePropsToTheLocals was expecting 6x1 vector of forces/moments');
        return;
    end
    
    if (size(thetas, 1) > 1)
        disp('Error: givePropsToTheLocals was expecting row vector of thetas, not column');
        return;
    end
    if (size(z, 2) > 1)
        disp('Error: givePropsToTheLocals was expecting column vector of z, not row');
        return;
    end
    if (size(thetas, 2) ~= size(z,1) - 1)
        disp('Error: givePropsToTheLocals requires thetas be k length and z be k+1 length ');
        return;
    end
        
    numPlies = length(thetas);  
    Qbars = [];
    GlobalStrains = [];
    LocalStrains = [];
    GlobalStresses = [];
    LocalStresses = [];
    Q0deg = getQ3(E1, E2, v12, v21, G12)

    for i = 1:numPlies
        Qbars = [Qbars, getQbar3(Q0deg, thetas(i))];
    end
    Qbars
    
    A = getA3(Qbars, z)
    B = getB3(Qbars, z)
    D = getD3(Qbars, z)
    
    MidplaneStrainCurv = getMidPlaneStrains(A,B,D,NM)
    
    for i = 1:length(z)
        GlobalStrains = [GlobalStrains, getGlobalStrains(MidplaneStrainCurv, z(i))];
    end
    Position = z'
    GlobalStrains
    
    Position = [];
    Ply = [];
    PlyNum = [];
    for i = 1:numPlies
        GlobalStresses = [GlobalStresses, getStressFromStrain(Qbars(:,i*3-2:i*3), GlobalStrains(:,i))];
        Position = [Position, z(i)];
        GlobalStresses = [GlobalStresses, getStressFromStrain(Qbars(:,i*3-2:i*3), GlobalStrains(:,i+1))];
        Position = [Position, z(i+1)];
        Ply = [Ply, thetas(i), thetas(i)];
        PlyNum = [PlyNum, i,i];
    end

    Ply_Angle_And_Position = [PlyNum; Ply; Position]
    GlobalStresses
    
    for i = 1:numPlies
        LocalStresses = [LocalStresses, stressGlobalToLocal3(GlobalStresses(:,i*2-1), getT3(thetas(i)))];
        LocalStresses = [LocalStresses, stressGlobalToLocal3(GlobalStresses(:,i*2), getT3(thetas(i)))];
    end
    Ply_Angle_And_Position = [PlyNum; Ply; Position]
    LocalStresses
    
    for i = 1:numPlies
        LocalStrains = [LocalStrains, getLocalStrainFromGlobal(GlobalStrains(:,i), getT3(thetas(i)))];
        LocalStrains = [LocalStrains, getLocalStrainFromGlobal(GlobalStrains(:,i+1), getT3(thetas(i)))];
    end
    Ply_Angle_And_Position = [PlyNum; Ply; Position]
    LocalStrains
end

function localStrain = getLocalStrainFromGlobal(globalStrain, T)
    R = [1,0,0;
        0,1,0;
        0,0,2];
    localStrain = R*T*inv(R)*globalStrain;
end

function matInv = getQorS(QorS)
    matInv = inv(QorS);
end

function strain = getStrainFromStress(S, stress)
    strain = S*stress;
end

function stress = getStressFromStrain(Q, strains)
    stress = Q*strains;
end

function globalStrains = getGlobalStrains(strainCurv, z)
    midplaneStrains = [strainCurv(1); strainCurv(2); strainCurv(3)];
    curvatures = [strainCurv(4); strainCurv(5); strainCurv(6)];
    globalStrains = midplaneStrains + z*curvatures;
end

function strainCurv = getMidPlaneStrains(A, B, D, NM)
    comb6 = [A, B;
             B, D];
    strainCurv = inv(comb6)*NM;
end

function forcesMoments = getForcesMoments(A,B,D,MP)
    comb6 = [A, B;
             B, D];
    forcesMoments = comb6*MP;
end

function A = getA3(Q, z)
    A = zeros(3,3);
    for i = 1:3
        for j = 1:3
            for k = 1:length(z)-1
                A(i,j) = A(i,j) + Q(i,j + (k-1)*3)*(z(k+1) - z(k));
            end
        end
    end
    
end


function B = getB3(Q, z)
    B = zeros(3,3);
    for i = 1:3
        for j = 1:3
            for k = 1:length(z)-1
                B(i,j) = B(i,j) + Q(i,j + (k-1)*3)*(z(k+1)^2 - z(k)^2);
            end
        end
    end
    B = 1/2 * B;
end

function D = getD3(Q,z)
    D = zeros(3,3);
    for i = 1:3
        for j = 1:3
            for k = 1:length(z)-1
                D(i,j) = D(i,j) + Q(i,j + (k-1)*3)*(z(k+1)^3 - z(k)^3);
            end
        end
    end
    D = 1/3 * D;
end


function localStress = stressGlobalToLocal3(globalStress, T)    
    localStress = T*globalStress;
end

function globalStress = stressLocalToGlobal3(localStress, T)
    globalStress = inv(T)*localStress;
end

function Qbar = getQbar3(Q, theta)
    m = cosd(theta);
    n = sind(theta);
    
    Qbar11 = Q(1,1)*m^4 + 2*(Q(1,2) + 2*Q(3,3))*n^2*m^2 + Q(2,2)*n^4;
    Qbar12 = (Q(1,1) + Q(2,2) - 4*Q(3,3))*n^2*m^2 + Q(1,2)*(n^4+m^4);
    Qbar22 = Q(1,1)*n^4 + 2*(Q(1,2) + 2*Q(3,3))*n^2*m^2 + Q(2,2)*m^4;
    Qbar16 = (Q(1,1) - Q(1,2) - 2*Q(3,3))*n*m^3 - (Q(2,2) - Q(1,2) - 2*Q(3,3))*n^3*m;
    Qbar26 = (Q(1,1) - Q(1,2) - 2*Q(3,3))*n^3*m - (Q(2,2) - Q(1,2) - 2*Q(3,3))*n*m^3;
    Qbar66 = (Q(1,1) + Q(2,2) - 2*Q(1,2) - 2*Q(3,3))*n^2*m^2 + Q(3,3)*(n^4+m^4);
    
    Qbar = [Qbar11, Qbar12, Qbar16;
            Qbar12, Qbar22, Qbar26;
            Qbar16, Qbar26, Qbar66];
end

function Q = getQ3(E1, E2, v12, v21, G12) 

Q = [E1/(1-v21*v12), v12*E2/(1-v21*v12) 0;
    v12*E2/(1-v21*v12), E2/(1-v21*v12), 0;
    0                   0               G12];
end

function S = getS3(E1, E2, v12, v21, G12) 
   S = [1/E1, -v12/E1, 0;
       -v21/E2, 1/E2, 0;
       0,       0,     1/G12];
end

function T = getT3(theta)
m = cosd(theta);
n = sind(theta);

T = [m^2, n^2, 2*m*n;
    n^2, m^2, -2*m*n;
    -m*n, m*n, m^2-n^2];
end
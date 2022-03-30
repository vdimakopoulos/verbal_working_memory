function [ x, y ] = Project_Cartesian_3D_On_2D( X, Y, Z )

[TH,PHI,R] = cart2sph(X,Y,Z);
PHI = pi/2-PHI;

Rsurface = PHI*mean(R);
x = Rsurface.*cos(TH);
y = Rsurface.*sin(TH);

end
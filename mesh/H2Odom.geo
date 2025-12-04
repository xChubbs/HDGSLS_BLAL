Lc = 0.125;

// Definition of points
Point(1) = {0, 0, 0, Lc};
Point(2) = {5, 0, 0, Lc};
Point(3) = {0, 5, 0, Lc};
Point(4) = {5, 5, 0, Lc};

// Definition of lines: First face
Line(1) = {1, 3};
Line(2) = {3, 4};
Line(3) = {4, 2};
Line(4) = {2, 1};

// Definition of control loops: (5) Positive Orientation (7) Negative Orientation
Line Loop(5) = {-3, -4, -1, -2};
Line Loop(7) = {3, 4, 1, 2};

Plane Surface(6) = {7};
Plane Surface(8) = {5};

// Extrude face 8 (Wrong orientation for loop but positive for the rest)
Extrude {0, 0, Pi} {
  Surface{8};
}

// Surface loop to create Volume
Surface Loop(31) = {6, 17, 21, 25, 29, 30};

// Physical properties
Physical Volume(40) = {1};

Physical Surface(50) = {6, 17, 21, 25, 29, 30};

// Definition of volume using Loop
// Volume(32) = {31};

// Physical Surface(44) = {15, 27, 28};
// Physical Surface(45) = {23, 6, 19};

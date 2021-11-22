SetFactory("OpenCASCADE");
// 
lc = 2.0;
//
Point(1) = {-1, -1, -1, lc};
Point(2) = {-1, 1, -1, lc};
Point(3) = {1, -1, 1, lc};
Point(4) = {1, 1, 1, lc};
//
Line(1) = {1, 2};
Line(2) = {4, 3};
Line(3) = {3, 1};
Line(4) = {2, 4};
//
Curve Loop(1) = {1, 4, 2, 3};
//
Plane Surface(1) = {1};
//
Physical Point("boundary", 1) = {1, 2, 3, 4};
Physical Curve("boundary", 2) = {1, 4, 2, 3};
Physical Surface("domain") = {1};
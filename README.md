# GeoAc
Infrasonic ray tracing code

GeoAc is a numerical package written in C++ which solves the equations governing acoustic 
propagation through the atmosphere in the geometric limit using a RK4 algorithm.  It contains
multiple instances of said equation system and is able to model propagation in an azimuthal 
plane using the effective motionless medium approximation as well as in three dimensions 
using an inhomogeneous moving background medium.  The three dimensional propagation scheme 
include methods to model propagation in a Cartesian coordinate system as well as a spherical 
coordinate system which incorporates the curvature of the earth. 

INSTALL:

1) To build locally, simply run "make all"

2) To build and install, check the install directory in the makefile (defaults to /usr/local/bin/) 
      and run "make install"

See manual for more information.

___________________________________________
___________________________________________

Copyright (c) 2014, Los Alamos National Security, LLC 
All rights reserved.

Copyright 2014. Los Alamos National Security, LLC. This software was produced under U.S. Government 
contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos 
National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, 
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC
MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If 
software is modified to produce derivative works, such modified software should be clearly marked, so 
as not to confuse it with the version available from LANL.

Additionally, permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to 
the following conditions:
 
The above copyright notice and this permission notice shall be included in all copies or substantial 
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO 
EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR 
THE USE OR OTHER DEALINGS IN THE SOFTWARE.

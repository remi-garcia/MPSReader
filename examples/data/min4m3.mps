NAME          min4m3
ROWS
 N  OBJ
 L  A1
 E  A2
 G  A3
COLUMNS
    x1        OBJ                 -2   A1                   1
    x1        A2                   1   A3                   1
    x2        OBJ                 -2   A1                   2
    x2        A3                   1   A2                  -1
    x3        A1                   1   OBJ                 -1
    x4        A2                   2   A3                  -1
    x3        A2                   1   A3                   1
    x4        A1                   1

RHS
    RHS       A1                  -5   A2                  -8
* A comment should be skip, just like the empty line before RHS section
    RHS       A3                   2
BOUNDS
 LO bnd       x2                   0
 FR bnd       x4
 FX bnd       x3                   3
 PL bnd       x4
 LO bnd       x4                   0
 LO bnd       x3                   0
 PL bnd       x3
ENDATA


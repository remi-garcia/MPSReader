NAME          min4m3
ROWS
 N  OBJ
 L  A1
 E  A2
 G  A3
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    x1        OBJ                 -2   A1                   1
    x1        A2                   1   A3                   1
    MARK0000  'MARKER'                 'INTEND'
    x2        OBJ                 -2   A1                   2
    x2        A3                   1   A2                  -1
    x3        A1                   1   OBJ                 -1
    x3        A2                   1   A3                   1
    MARK0001  'MARKER'                 'INTORG'
    x4        A2                   2   A3                  -1
    x4        A1                   1
    MARK0001  'MARKER'                 'INTEND'
RHS
    RHS       A1                  -5   A2                  -8
    RHS       A3                   2
BOUNDS
 UP bnd       x1                   1
 BV bnd       x2
 UI bnd       x3                   1
 UP bnd       x4                   1
ENDATA


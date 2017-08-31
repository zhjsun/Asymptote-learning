settings.outformat = "pdf";
import CJCBI;
import graph;
size(150);

real x = 2;
real y = 1;

real[][] A = {{x^x,sqrt(0.5)},{sqrt(0.5),y^y}};
int n = 2;
real[][] B = A;
real[][] V = {{0,0},{0,0}};
real eps = 1e-6;
int jt = 100;
int bo = cjcbi(B, n, V, eps, jt);
write(bo);
write(B);
draw(rotate(degrees(atan(V[1][0]/V[0][0])))*ellipse((0,0),sqrt(B[0][0]),sqrt(B[1][1])),red);
draw(box((-x,-y),(x,y)),blue);
xaxis("$x$",Arrow);
yaxis("$y$",Arrow);

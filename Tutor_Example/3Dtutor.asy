//settings.outformat = "pdf";
//settings.prc = false;
//settings.render = 0;
import graph;
//import graph3;
import three;
currentprojection= orthographic(5,0,10, up=Y);//obliqueX;//perspective(1,1,1);
real unit = 2.1cm;
real truecm = cm / unit;
unitsize(unit);

// Save some important numbers
real xmin = -0.1;
real xmax = 2;
real ymin = -0.1;
real ymax = 2;
real zmin = -0.1;
real zmax = 2;
real margin = 0.4 truecm;

// Construct the graph
real f(real x) { return sqrt(x); }
path s = graph(f, 0, 2, operator..);
path fillregion = s -- (xmax,0) -- cycle;
path3 p3 = path3(s);
surface solidsurface = surface(p3, c = O, axis = X);

// Draw the graph
draw(solidsurface, yellow);
draw(p3, black);

// Draw the plane
path planeoutline = box((xmin, ymin), (xmax+margin, ymax+margin));
draw(surface(planeoutline ^^ fillregion), surfacepen=lightgray, light=nolight);
// Fill the area under the graph
draw(surface(fillregion), surfacepen=gray(0.6), light=nolight);

// Draw the the axes
draw(xmin*X -- xmax*X, arrow=Arrow3(TeXHead2(normal=Z), emissive(blue)), blue);
draw(ymin*Y -- ymax*Y, arrow=Arrow3(TeXHead2(normal=Z), emissive(green)), green);
//draw(zmin*Z -- zmax*Z, L = Label("$z$", position=EndPoint), arrow=Arrow3, red);

//Direction of a point toward the camera.
triple cameradirection(triple pt, projection P=currentprojection) {
    if (P.infinity) {
        return unit(P.camera);
    } else {
        return unit(P.camera - pt);
    }
}

//Move a point closer to the camera.
triple towardcamera(triple pt, real distance=1, projection P=currentprojection) {
    return pt + distance * cameradirection(pt, P);
}

// Label the axes
label("$x$", align=E, position=towardcamera(xmax*X), p=blue);
label("$y$", align=N, position=towardcamera(ymax*Y), p=green);

settings.outformat = "pdf";
import CJCBI;
import graph;
size(350, 300, IgnoreAspect);
//unitsize(x = 0.0001cm, y = 0.0001cm);
locale("en_US.UTF-8");

int dim = 6,
    Xindex = 0,
    Yindex = 1;
real Mtemp[] = new real [dim],
     Ptemp[][] = new real [dim][dim],
     Mtx_2D[][] = new real [2][2],
     Vec_2D[][] = new real [2][2],
     eps = 1e-6,
     xscale = 1e-6,
     yscale = 1e-6;

////////////////////////////////////////////////////////////////////////////
// load files
file fin = input("GMM_Final.dat");
real[][] GMM_MP = fin.dimension(0,dim*(dim+1));
close(fin);
int nGMM = GMM_MP.length;

file fin = input("MC_Points_Final.dat");
real[][] MC_Pt = fin.dimension(0,dim);
close(fin);
int nMC = MC_Pt.length;

////////////////////////////////////////////////////////////////////////////
// Print GMM
path p_GMM[];
int jt = 100;
for(int i = 0; i < nGMM; ++i) {
    // Initialize
    for(int j = 0; j < dim; ++j) {
        Mtemp[j] = GMM_MP[i][j];
        for(int k = 0; k < dim; ++k) {
            Ptemp[j][k] = GMM_MP[i][dim+dim*j+k];
        }
    }
    pair cent = (Mtemp[Xindex], Mtemp[Yindex]);
    Mtx_2D[0][0] = Ptemp[Xindex][Xindex];
    Mtx_2D[0][1] = Ptemp[Xindex][Yindex];
    Mtx_2D[1][0] = Ptemp[Yindex][Xindex];
    Mtx_2D[1][1] = Ptemp[Yindex][Yindex];

    // Eigen value and vector
    int flag = cjcbi(Mtx_2D, 2, Vec_2D, eps, jt);

    // Plot the ellipse
    if(flag != 0) {
        path ellip = rotate(degrees(atan(Vec_2D[1][0]/Vec_2D[0][0])), cent) *
            ellipse(cent, 3*sqrt(Mtx_2D[0][0]), 3*sqrt(Mtx_2D[1][1]) );
        p_GMM = p_GMM ^^ ellip;
    }
}

////////////////////////////////////////////////////////////////////////////
// Print MC Points
path p_MC[];
for(int i = 0; i < nMC; ++i) {
    p_MC = p_MC ^^ (MC_Pt[i][Xindex], MC_Pt[i][Yindex]);
}

///////////////////////////////////////////////////////////////////////////////
// Define pens and markers
pen p1 = rgb(0.8, 0.4, 0.2);
pen p2 = rgb(1, 1, 1);
pen pmark = rgb(0.2, 0.4, 0.8);

////////////////////////////////////////////////////////////////////////////////
// Draw
draw(scale(xscale, yscale)*p_GMM, p1, Label(s="GMM ellipse", p=fontsize(8)) );
draw(scale(xscale, yscale)*p_MC, p2, Label(s="MC point", p=fontsize(8)), 
        marker(scale(0.5mm)*unitcircle, FillDraw(pmark)) );

// xaxis(L="$x/10^6$", axis=BottomTop, ticks=LeftTicks(N=2, n=5), p=fontsize(10));
// yaxis(L=shift(20N+50E)*Label("$y/10^6$",1), axis=LeftRight, ticks=RightTicks(N=2, n=5), p=fontsize(10));
// add(legend(linelength=30), 10SE, UnFill);

real xmin = 4.8,
     xmax = 6,
     ymin = -4.8,
     ymax = -3.2;
limits((xmin, ymin), (xmax, ymax), Crop);
draw(box((xmin,ymin),(xmax,ymax)));
xaxis(L="$x/10^6$", xmin=xmin, xmax=xmax, axis=BottomTop, 
    ticks=LeftTicks(N=2, n=5), p=fontsize(10));
yaxis(L=shift(20N+50E)*Label("$y/10^6$",1), ymin=ymin, ymax=ymax, axis=LeftRight, 
    ticks=RightTicks(N=2, n=6), p=fontsize(10));
add(legend(linelength=30), (xmin, ymax), 10SE, UnFill);

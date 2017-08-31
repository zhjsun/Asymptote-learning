settings.outformat = "pdf";
import CJCBI;
import graph;
size(350, 300, IgnoreAspect);
//unitsize(x = 0.0001cm, y = 0.0001cm);
locale("en_US.UTF-8");


int dim = 6,
    Xindex = 0,
    Yindex = 1;
int jt = 100;
real Mtemp[] = new real [dim],
     Ptemp[][] = new real [dim][dim],
     Mtx_2D[][] = new real [2][2],
     Vec_2D[][] = new real [2][2],
     eps = 1e-6,
     xscale = 1e-6,
     yscale = 1e-6;

////////////////////////////////////////////////////////////////////////////
// load files
file fin = input("MC_Points_Final.dat");
real[][] MC_Pt = fin.dimension(0,dim);
close(fin);
int nMC = MC_Pt.length;

file fin = input("MC_MP_Final.dat");
real[] MC_MP = fin;
close(fin);

file fin = input("DA_MP.dat");
real[][] DA_MP = fin.dimension(0,dim*(dim+1)+1);
close(fin);
int nDA = DA_MP.length;
nDA = 4;

/////////////////////////////////////////////////////////////////////////////
// Print MC Points
path p_MC[];
for(int i = 0; i < nMC; ++i) {
    p_MC = p_MC ^^ (MC_Pt[i][Xindex], MC_Pt[i][Yindex]);
}

////////////////////////////////////////////////////////////////////////////////
// Print MC Covariance
pair cent_MC;
path ellip_MC;
// Initialize
for(int j = 0; j < dim; ++j) {
    Mtemp[j] = MC_MP[j];
    for(int k = 0; k < dim; ++k) {
        Ptemp[j][k] = MC_MP[dim+dim*j+k];
    }
}
cent_MC = (Mtemp[Xindex], Mtemp[Yindex]);
Mtx_2D[0][0] = Ptemp[Xindex][Xindex];
Mtx_2D[0][1] = Ptemp[Xindex][Yindex];
Mtx_2D[1][0] = Ptemp[Yindex][Xindex];
Mtx_2D[1][1] = Ptemp[Yindex][Yindex];

// Eigen value and vector
int flag = cjcbi(Mtx_2D, 2, Vec_2D, eps, jt);

// Plot the ellipse
if(flag != 0) {
    ellip_MC = rotate(degrees(atan(Vec_2D[1][0]/Vec_2D[0][0])), cent_MC) * 
        ellipse(cent_MC, 3*sqrt(Mtx_2D[0][0]), 3*sqrt(Mtx_2D[1][1]));
}

////////////////////////////////////////////////////////////////////////////////
// Print DA Covariance
pair cent_DA[] = new pair [nDA];
path ellip_DA[] = new path [nDA];
for(int i = 0; i < nDA; ++i) {
    // Initialize
    for(int j = 0; j < dim; ++j) {
        Mtemp[j] = DA_MP[i][j+1];
        for(int k = 0; k < dim; ++k) {
            Ptemp[j][k] = DA_MP[i][dim+1+dim*j+k];
        }
    }
    cent_DA[i] = (Mtemp[Xindex], Mtemp[Yindex]);
    Mtx_2D[0][0] = Ptemp[Xindex][Xindex];
    Mtx_2D[0][1] = Ptemp[Xindex][Yindex];
    Mtx_2D[1][0] = Ptemp[Yindex][Xindex];
    Mtx_2D[1][1] = Ptemp[Yindex][Yindex];

    // Eigen value and vector
    int flag = cjcbi(Mtx_2D, 2, Vec_2D, eps, jt);

    // Plot the ellipse
    if(flag != 0) {
        ellip_DA[i] = rotate(degrees(atan(Vec_2D[1][0]/Vec_2D[0][0])), cent_DA[i]) * 
            ellipse(cent_DA[i], 3*sqrt(Mtx_2D[0][0]), 3*sqrt(Mtx_2D[1][1]));
    }
}

///////////////////////////////////////////////////////////////////////////////
// Define pens and markers
pen[] pen_cov = new pen [nDA+1];
pen_cov[0] = heavyred + dotted + 1.0;
pen_cov[1] = orange + dashed + 1.0;
pen_cov[2] = deepgreen + linetype(new real[] {18, 6}) + 1.0;
pen_cov[3] = heavyblue + dashdotted + 1.0;
pen_cov[4] = olive + solid + 0.8;

marker[] mark_mean = new marker [nDA+1];
mark_mean[0] = marker(scale(1.0 mm)*cross(3), heavyred, FillDraw(blue));
mark_mean[1] = marker(scale(1.0 mm)*cross(4), orange, FillDraw(blue));
mark_mean[2] = marker(scale(1.0 mm)*cross(5), deepgreen, FillDraw(blue));
mark_mean[3] = marker(scale(1.0 mm)*cross(6), heavyblue, FillDraw(blue));
mark_mean[4] = marker(scale(1.5 mm)*polygon(3));

pen pmark = rgb(0.2, 0.4, 0.8);

////////////////////////////////////////////////////////////////////////////////
// Draw
draw(scale(xscale, yscale)*p_MC, invisible, Label(s="Monte Carlo sample", p=fontsize(8)), 
    marker(scale(0.3mm)*unitcircle, linewidth(0.1), FillDraw(pmark)) );
draw(scale(xscale, yscale)*ellip_MC, pen_cov[nDA], Label("Monte Carlo cov", p=fontsize(8)) );
draw(scale(xscale, yscale)*cent_MC, invisible, Label("Monte Carlo mean", p=fontsize(8)), mark_mean[nDA] );
for(int i = 0; i < nDA; ++i) {
    draw(scale(xscale, yscale)*ellip_DA[i], pen_cov[i], Label(s=string(i+1) + "-order cov", p=fontsize(8)) );
    draw(scale(xscale, yscale)*cent_DA[i], invisible, 
        Label(s=string(i+1) + "-order mean", p=fontsize(8)), mark_mean[i] );
}

// xaxis(L="$x/10^6$", axis=BottomTop, ticks=LeftTicks(N=2, n=5), p=fontsize(10));
// yaxis(L=shift(20N+50E)*Label("$y/10^6$",1), axis=LeftRight, ticks=RightTicks(N=2, n=5), p=fontsize(10));

real xmin = 4.8,
     xmax = 6,
     ymin = -4.8,
     ymax = -3.2;
limits((xmin, ymin), (xmax, ymax), Crop);
draw(box((xmin,ymin),(xmax,ymax)));
xaxis(L="$x/10^6$", xmin=xmin, xmax=xmax, axis=BottomTop, 
    ticks=LeftTicks(N=2, n=5), p=fontsize(10));
yaxis(L=shift(20N+50E)*Label("$y/10^6$",1), ymin=ymin, ymax=ymax, axis=LeftRight, 
    ticks=RightTicks(N=2, n=5), p=fontsize(10));
add(legend(linelength=30), (xmin, ymax), 10SE, UnFill);

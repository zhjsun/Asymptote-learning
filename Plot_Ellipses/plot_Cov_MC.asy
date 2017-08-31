settings.outformat = "pdf";
import CJCBI;
import graph;
size(350, 300, IgnoreAspect);
//unitsize(x = 0.0001cm, y = 0.0001cm);
locale("en_US.UTF-8");


int dim = 6,
    Xindex = 1,
    Yindex = 0;
int nNode = 6;
real xscale = 1e-0,
     yscale = 1e-0;
////////////////////////////////////////////////////////////////////////////
// load files

// Monte Carlo points
file fin = input("MC_Points_1.dat");
real[][] MC_Pt_1 = fin.dimension(0,dim);
close(fin);
int nMC = MC_Pt_1.length;

file fin = input("MC_Points_2.dat");
real[][] MC_Pt_2 = fin.dimension(0,dim);
close(fin);

file fin = input("MC_Points_3.dat");
real[][] MC_Pt_3 = fin.dimension(0,dim);
close(fin);

file fin = input("MC_Points_4.dat");
real[][] MC_Pt_4 = fin.dimension(0,dim);
close(fin);

file fin = input("MC_Points_5.dat");
real[][] MC_Pt_5 = fin.dimension(0,dim);
close(fin);

file fin = input("MC_Points_6.dat");
real[][] MC_Pt_6 = fin.dimension(0,dim);
close(fin);

// Trajectory state
file fin = input("Traj_state.dat");
real[][] state = fin.dimension(0,dim+1);
close(fin);
int nState = state.length;

// Mean and covariance
file fin = input("Cov_analytical.dat");
real[][] MP_ana = fin.dimension(0,dim*(dim+1)+1);
close(fin);

file fin = input("Cov_MonteCarlo.dat");
real[][] MP_MC = fin.dimension(0,dim*(dim+1)+1);
close(fin);

////////////////////////////////////////////////////////////////////////////////
// Path of relative trajectory
path p_state;
for(int i = 0; i < nState; ++i) {
    p_state = p_state .. (state[i][Xindex+1], state[i][Yindex+1]);
}

////////////////////////////////////////////////////////////////////////////////
// Path of ellipses
path ellip_MC[] = new path [nNode];
path ellip_ana[] = new path [nNode];
pair cent_MC[] = new pair [nNode];
pair cent_ana[] = new pair [nNode];

real Mtemp[] = new real [dim],
     Ptemp[][] = new real [dim][dim],
     Mtx_2D[][] = new real [2][2],
     Vec_2D[][] = new real [2][2];

for(int i = 0; i < nNode; ++i) {
    // Monte Carlo ellipses
    // Initialize
    for(int j = 0; j < dim; ++j) {
        Mtemp[j] = MP_MC[i][j+1];
        for(int k = 0; k < dim; ++k) {
            Ptemp[j][k] = MP_MC[i][dim+1+dim*j+k];
        }
    }
    cent_MC[i] = (Mtemp[Xindex], Mtemp[Yindex]);
    ellip_MC[i] = Solve_Ellipse(Mtemp, Ptemp, Xindex, Yindex);

    // Analytical ellipses
    // Initialize
    for(int j = 0; j < dim; ++j) {
        Mtemp[j] = MP_ana[i][j+1];
        for(int k = 0; k < dim; ++k) {
            Ptemp[j][k] = MP_ana[i][dim+1+dim*j+k];
        }
    }
    cent_ana[i] = (Mtemp[Xindex], Mtemp[Yindex]);
    ellip_ana[i] = Solve_Ellipse(Mtemp, Ptemp, Xindex, Yindex);
}

////////////////////////////////////////////////////////////////////////////////
// Path of MC points
path p_MC_1[];
for(int i = 0; i < nMC; ++i) {
    p_MC_1 = p_MC_1 ^^ (MC_Pt_1[i][Xindex], MC_Pt_1[i][Yindex]);
}

path p_MC_2[];
for(int i = 0; i < nMC; ++i) {
    p_MC_2 = p_MC_2 ^^ (MC_Pt_2[i][Xindex], MC_Pt_2[i][Yindex]);
}

path p_MC_3[];
for(int i = 0; i < nMC; ++i) {
    p_MC_3 = p_MC_3 ^^ (MC_Pt_3[i][Xindex], MC_Pt_3[i][Yindex]);
}

path p_MC_4[];
for(int i = 0; i < nMC; ++i) {
    p_MC_4 = p_MC_4 ^^ (MC_Pt_4[i][Xindex], MC_Pt_4[i][Yindex]);
}

path p_MC_5[];
for(int i = 0; i < nMC; ++i) {
    p_MC_5 = p_MC_5 ^^ (MC_Pt_5[i][Xindex], MC_Pt_5[i][Yindex]);
}

path p_MC_6[];
for(int i = 0; i < nMC; ++i) {
    p_MC_6 = p_MC_6 ^^ (MC_Pt_6[i][Xindex], MC_Pt_6[i][Yindex]);
}

///////////////////////////////////////////////////////////////////////////////
// Define pens and markers
pen[] pen_cov = new pen [5];
pen_cov[0] = heavyred + dotted + 1.0;
pen_cov[1] = orange + dashed + 1.0;
pen_cov[2] = deepgreen + linetype(new real[] {18, 6}) + 1.0;
pen_cov[3] = heavyblue + dashdotted + 1.0;
pen_cov[4] = olive + solid + 0.8;

marker[] mark_mean = new marker [5];
mark_mean[0] = marker(scale(1.0 mm)*cross(3), heavyred, FillDraw(blue));
mark_mean[1] = marker(scale(1.0 mm)*cross(4), orange, FillDraw(blue));
mark_mean[2] = marker(scale(1.0 mm)*cross(5), deepgreen, FillDraw(blue));
mark_mean[3] = marker(scale(1.0 mm)*cross(6), heavyblue, FillDraw(blue));
mark_mean[4] = marker(scale(1.5 mm)*polygon(3));

pen pmark = rgb(0.2, 0.4, 0.8);

////////////////////////////////////////////////////////////////////////////////
// Draw
draw(scale(xscale, yscale)*p_state, pen_cov[0], Label("Trajectory", p=fontsize(8)) );

draw(scale(xscale, yscale)*p_MC_5, invisible, Label(s="Monte Carlo sample", p=fontsize(8)), 
    marker(scale(0.3mm)*unitcircle, linewidth(0.1), FillDraw(pmark)) );
draw(scale(xscale, yscale)*ellip_MC[4], pen_cov[4], Label("Monte Carlo cov", p=fontsize(8)) );
draw(scale(xscale, yscale)*ellip_ana[4], pen_cov[3], Label("Analytical cov", p=fontsize(8)) );
// draw(scale(xscale, yscale)*cent_MC, invisible, Label("Monte Carlo mean", p=fontsize(8)), mark_mean[4] );

draw(scale(xscale, yscale)*p_MC_4, invisible, marker(scale(0.3mm)*unitcircle, linewidth(0.1), FillDraw(pmark)) );
draw(scale(xscale, yscale)*ellip_MC[3], pen_cov[4]);
draw(scale(xscale, yscale)*ellip_ana[3], pen_cov[3]);

draw(scale(xscale, yscale)*p_MC_3, invisible, marker(scale(0.3mm)*unitcircle, linewidth(0.1), FillDraw(pmark)) );
draw(scale(xscale, yscale)*ellip_MC[2], pen_cov[4]);
draw(scale(xscale, yscale)*ellip_ana[2], pen_cov[3]);

draw(scale(xscale, yscale)*p_MC_2, invisible, marker(scale(0.3mm)*unitcircle, linewidth(0.1), FillDraw(pmark)) );
draw(scale(xscale, yscale)*ellip_MC[1], pen_cov[4]);
draw(scale(xscale, yscale)*ellip_ana[1], pen_cov[3]);

int flag = 1;
if(flag == 0) {
    xaxis(L="V-bar (m)", axis=BottomTop, ticks=LeftTicks(N=2, n=5), p=fontsize(10));
    yaxis(L="R-bar (m)", axis=LeftRight, ticks=RightTicks(N=2, n=5), p=fontsize(10));
} else {
    real xmin = -32000,
         xmax = 6000,
         ymin = -8000,
         ymax = 2000;
    limits((xmin, ymin), (xmax, ymax), Crop);
    xaxis(L="V-bar (m)", xmin=xmin, xmax=xmax, axis=BottomTop, 
        ticks=LeftTicks(N=2, n=5), p=fontsize(10));
    yaxis(L="R-bar (m)", ymin=ymin, ymax=ymax, axis=LeftRight, 
        ticks=RightTicks(N=2, n=5), p=fontsize(10));
}

add(legend(linelength=30),  point(NW), 10SE, UnFill);

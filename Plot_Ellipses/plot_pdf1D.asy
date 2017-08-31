settings.outformat = "pdf";
import CJCBI;
import graph;
size(350, 300, IgnoreAspect);
//unitsize(x = 0.0001cm, y = 0.0001cm);
locale("en_US.UTF-8");

int dim = 6,
    index = 2,
    r = 5,
    ns = 80;
real Mtemp[] = new real [dim],
     Ptemp[][] = new real [dim][dim],
     xscale = 1e-4,
     yscale = 1e5;

////////////////////////////////////////////////////////////////////////////
// load files
file fin = input("MC_Points_Final.dat");
real[][] MC_Pt = fin.dimension(0,dim);
close(fin);
int nMC = MC_Pt.length;

file fin = input("GMM_Final.dat");
real[][] GMM_MP = fin.dimension(0,dim*(dim+1));
close(fin);
int nGMM = GMM_MP.length;

real MC_vec[] = new real [nMC];
for(int i = 0; i < nMC; ++i) {
    MC_vec[i] = MC_Pt[i][index];
}
real Dmin = min(MC_vec),
     Dmax = max(MC_vec),
     Dstep = (Dmax - Dmin)/ns,
     Ds[] = uniform(Dmin, Dmax+r*Dstep, ns+r),
     f_MC[] = array(ns+r+1, 0);

/////////////////////////////////////////////////////////////////////////////
// Print GMM PDF
real w = 1/nGMM;

real pdf(real x, real[][] GMM_MP) {
    real f = 0;
    for(int i = 0; i < nGMM; ++i) {
        // Initialize
        for(int j = 0; j < dim; ++j) {
            Mtemp[j] = GMM_MP[i][j];
            for(int k = 0; k < dim; ++k) {
                Ptemp[j][k] = GMM_MP[i][dim+dim*j+k];
            }
        }
        real M_GMM = Mtemp[index],
             P_GMM = Ptemp[index][index];
        f = f + 1/sqrt(2*pi*P_GMM)*exp(-0.5*(x - M_GMM)^2/P_GMM);
    }
    f = f*w;
    return f;
}

real pdf_1(real x){
    return pdf(x, GMM_MP);
}

path s_1 = graph(pdf_1, Dmin, Dmax + 0.1*(Dmax-Dmin), Hermite(monotonic) );

/////////////////////////////////////////////////////////////////////////////
// Print MC PDF
for(int i = 0; i < nMC; ++i) {
    int n = floor( (MC_vec[i]-Dmin)/Dstep + 0.5);
    f_MC[n] += 1;
}
f_MC = f_MC/nMC/Dstep;
path p_MC;
for(int i = 0; i < ns+r+1; ++i) {
    p_MC = p_MC -- (Ds[i], f_MC[i]);
}

///////////////////////////////////////////////////////////////////////////////
// Define pens and markers
pen[] pens = new pen [5];
pens[0] = heavyred + dotted + 1.0;
pens[1] = orange + dashed + 1.0;
pens[2] = deepgreen + linetype(new real[] {18, 6}) + 1.0;
pens[3] = heavyblue + dashdotted + 1.0;
pens[4] = olive + solid + 0.8;

/////////////////////////////////////////////////////////////////////////////
// Draw
draw(scale(xscale, yscale)*s_1, pens[0], Label(s="DA\_GMM", p=fontsize(8)) );
draw(scale(xscale, yscale)*p_MC, pens[4], Label(s="Monte Carlo", p=fontsize(8)) );

// xaxis(L="$x$", axis=BottomTop, ticks=LeftTicks(N=2, n=5), p=fontsize(10));
// yaxis(L="Probability density", axis=LeftRight, ticks=RightTicks(N=2, n=5), p=fontsize(10));
// add(legend(linelength=30), point(NW), 10SE, UnFill);

real xmin = -2.4,
     xmax = 1.2,
     ymin = 0,
     ymax = 9;
limits((xmin, ymin), (xmax, ymax), Crop);
// draw(box((xmin,ymin),(xmax,ymax)));
xaxis(L="$z/10^4$ (m)", xmin=xmin, xmax=xmax, axis=BottomTop, 
        ticks=LeftTicks(N=3, n=4), p=fontsize(10));
yaxis(L="Probability density /$10^{-5}$", ymin=ymin, ymax=ymax, axis=LeftRight, 
        ticks=RightTicks(N=2, n=3), p=fontsize(10));
add(legend(linelength=30), point(NW), 10SE, UnFill);
settings.outformat = "pdf";
size(8cm, 0);
defaultpen(fontsize(10pt));
import graph;

// Define the command drawshifted
void drawshifted(path g, pair trueshift, picture pic = currentpicture, Label label="",
                pen pen=currentpen, arrowbar arrow=None, arrowbar bar=None,
                margin margin=NoMargin, marker marker=nomarker)
{
    pic.add(new void(frame f, transform t)
    {
        picture opic;
        draw(opic, L=label, shift(trueshift)*t*g, p=pen, arrow=arrow, bar=bar, margin=margin, marker=marker);
        add(f,opic.fit());
    });
    pic.addBox(min(g), max(g), trueshift+min(pen), trueshift+max(pen));
}

// Save some important numbers
real xmin = -0.1;
real xmax = 2;
real ymin = -0.1;
real ymax = 2;

// Draw the graph and fill the area under it
real f(real x) { return sqrt(x); }
path s = graph(f, 0, 2, operator..);
pen fillpen = mediumgray;
fill(s--(xmax,0)--cycle, fillpen);
draw(s, L=Label("$y=f(x)$", position=EndPoint));

// Fill the strip of width dx
real x = 1.4;
real dx = .05;
real t0 = times(s, x)[0];
real t1 = times(s, x+dx)[0];
path striptop = subpath(s, t0, t1);
filldraw((x,0)--striptop--(x+dx,0)--cycle, black);

// Draw the bars labeling the width dx
real barheight = f(x+dx);
pair barshifty = (0, 0.2cm);
Label dxlabel = Label("$dx$", position=MidPoint, align=2N);
drawshifted((x,barheight)--(x+dx,barheight), trueshift=barshifty, label=dxlabel, bar=Bars);

// Draw the arrows pointing inward toward the dx label
real myarrowlength = 0.3cm;
margin arrowmargin = DotMargin;
path leftarrow = shift(barshifty)*((-myarrowlength,0)--(0,0));
path rightarrow = shift(barshifty)*((myarrowlength,0)--(0,0));
draw((x,barheight), leftarrow, arrow=Arrow(), margin=arrowmargin);
draw((x+dx,barheight), rightarrow, arrow=Arrow(), margin=arrowmargin);

// Draw the bar labeling the height f(x)
real barx = x + dx;
pair barshiftx = (0.42cm,0);
Label fxlabel = Label("$f(x)$", align=(0,0), position=MidPoint, filltype=Fill(fillpen));
drawshifted((barx,0)--(barx,f(x)), trueshift=barshiftx, label=fxlabel, arrow=Arrows(), bar=Bars);

// Draw the axes on top of everything
arrowbar axisarrow = Arrow(TeXHead);
Label xlabel = Label("$x$", position=EndPoint);
draw((xmin,0)--(xmax,0), arrow=axisarrow, L=xlabel);
Label ylabel = Label("$y$", position=EndPoint);
draw((0,ymin)--(0,ymax), arrow=axisarrow, L=ylabel);

// Draw the tick mark on the x-axis
path tick = (0,0)--(0,-0.15cm);
Label ticklabel = Label("$x$", position=EndPoint);
draw((x,0), tick, L=ticklabel);


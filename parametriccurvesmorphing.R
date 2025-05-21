# Parametric Curves Morphing
# www.overfitting.net
# https://www.overfitting.net/2025/05/morphing-de-curvas-parametricas-con-r.html

library(png)
library(Cairo)


###################################
# Lorenz attractor

library(deSolve)  # ordinary differential equations solver (ode)
library(plotly)  # high quality atialiased output


# Parameters:  a       b        c
prm=list( sigma=10, rho=28, beta=8/3 )

# Init values
varini=c(X=1, Y=1, Z=1)

Lorenz=function (t, vars, prm) {
    with(as.list(vars), {
        dX=prm$sigma*(Y - X)
        dY=X*(prm$rho - Z) - Y
        dZ=X*Y - prm$beta*Z
        return(list(c(dX, dY, dZ)))
    })
}

N=100000 # number of points
times=seq(from=0, to=100, length.out=N)
# Ordinary differential equations solver (ode)
out=ode(y=varini, times=times, func=Lorenz, parms=prm)

# Point colour assignation
gfill=function (repArr, long) {
    rep(repArr, ceiling(long/length(repArr)))[1:long]
}

dout=as.data.frame(out)
dout$colour=gfill(gray.colors(nrow(dout)), nrow(dout))

# Plot Lorenz attractor using plotly
fig=plot_ly(
    data=dout, x=~X, y=~Y, z=~Z,
    type='scatter3d', mode='lines',
    opacity=1, line=list(width=2, color=~colour, reverscale=FALSE)
)

# Apply isometric view
fig=fig %>% layout(
    scene=list(
        camera=list(
            # eye=list(x=1.25, y=1.25, z=1.25),  # diagonal top-front-right view
            eye=list(x=0, y=10, z=0),  # classical Lorenz XZ view
            projection=list(type="orthographic")  # removes perspective distortion
        )
    )
)

fig


# HTML output version
htmlwidgets::saveWidget(fig, "lorenzattractor.html", selfcontained=TRUE)
browseURL("lorenzattractor.html")



###################################
# Parameters

N=100000 # number of points
NFRAMES=50  #283  # number of frames per transition -> 1981 total frames
NMORPH=7
# Create Nx3 matrix to hold the parametric curves
X=matrix(nrow=N, ncol=3)  # X[,1]=start, X[,2]=frame, X[,3]=end
Y=X  # Y[,1]=start, Y[,2]=frame, Y[,3]=end

t=seq(0, 1, length.out=N)


###################################
# Animation

# Full HD
DIMX=1920 # 512
DIMY=1080 # 512
GAP=ifelse(DIMX==1920, 5, 1)  # plot border
LWD=ifelse(DIMX==1920, 2, 1)  # lines thickness

# All images normalized to xlim=c(-1,1), ylim=c(-1,1)
fm=0  # general frame count
for (morph in 1:NMORPH) {
    if (morph==1) {  # circle -> spiral
        X[,1]=cos(2*pi*t)
        Y[,1]=sin(2*pi*t)
        
        X[,3]=cos(2*pi*50*t)*t
        Y[,3]=sin(2*pi*50*t)*t
    } else if (morph==2) {  # -> grid
        X[,1]=X[,3]
        Y[,1]=Y[,3]
        
        X[,3]=cos(2*pi*35*t)
        Y[,3]=sin(2*pi*40*t)
    } else if (morph==3) {  # -> flower1
        X[,1]=X[,3]
        Y[,1]=Y[,3]
        
        X[,3]=sin(8*pi*2*t)*cos(2*pi*t)
        Y[,3]=sin(8*pi*2*t)*sin(2*pi*t)
    } else if (morph==4) {  # -> flower2
        X[,1]=X[,3]
        Y[,1]=Y[,3]
        
        t2=(t*16-8)
        X[,3]=6 * sin(9.52*t2) * round(cos(cos(4.8*t2))^0.5)
        Y[,3]=6 * cos(9.52*t2)^4 * sin(sin(4.8*t2))
        
        X[,3]=(X[,3]-min(X[,3]))/(max(X[,3])-min(X[,3]))*2-1
        Y[,3]=(Y[,3]-min(Y[,3]))/(max(Y[,3])-min(Y[,3]))*2-1
    } else if(morph==5) {  # -> butterfly
        X[,1]=X[,3]
        Y[,1]=Y[,3]        
        
        t2=t*12*pi
        X[,3]=sin(t2) * ( exp(cos(t2)) - 2*cos(4*t2) - sin(t2/12)^5 )
        Y[,3]=cos(t2) * ( exp(cos(t2)) - 2*cos(4*t2) - sin(t2/12)^5 )
        
        X[,3]=(X[,3]-min(X[,3]))/(max(X[,3])-min(X[,3]))*2-1
        Y[,3]=(Y[,3]-min(Y[,3]))/(max(Y[,3])-min(Y[,3]))*2-1 
    } else if (morph==6) {  # -> lorenz
        X[,1]=X[,3]
        Y[,1]=Y[,3]
        
        X[,3]=(out[,2]-min(out[,2]))/(max(out[,2])-min(out[,2]))*2-1
        Y[,3]=(out[,4]-min(out[,4]))/(max(out[,4])-min(out[,4]))*2-1  
    } else if (morph==7) {  # -> circle (initial image)
        X[,1]=X[,3]
        Y[,1]=Y[,3]
        
        X[,3]=cos(2*pi*t)
        Y[,3]=sin(2*pi*t)   
    }
    
    for (f in 0:(NFRAMES-1)) {
        name=paste0("curves", ifelse(fm<10, "000", ifelse(fm<100, "00",
                              ifelse(fm<1000, "0", ""))), fm, ".png")
        print(paste0(fm+1, "/", NFRAMES*NMORPH, ": Writing '", name, "'..."))
        alpha=f/(NFRAMES-1)  # alpha loops through 0..1
        X[,2] = (1-alpha)*X[,1] + alpha*X[,3]  # linear interpolation
        Y[,2] = (1-alpha)*Y[,1] + alpha*Y[,3]
        CairoPNG(name, width=DIMX, height=DIMY, antialias="subpixel")
            par(mar=c(0,0,0,0), oma=c(GAP,GAP,GAP,GAP))
            plot(X[,2], Y[,2], xlim=c(-1,1), ylim=c(-1,1), type='l', lwd=LWD,
                 axes=FALSE, xlab="", ylab="", main="", xaxs="i", yaxs="i")
        dev.off()
        writePNG(1-readPNG(name)[,,1], name)  # invert and monochrome PNG
        fm=fm+1
    }
}


#############################
# Add Full HD title...

bkg=readPNG("background.png")
for (fm in 0:1980) {
    name=paste0("curves", ifelse(fm<10, "000", ifelse(fm<100, "00",
                          ifelse(fm<1000, "0", ""))), fm, ".png")
    img=readPNG(name)
    img[img<bkg]=bkg[img<bkg]  # add background
    
    name2=paste0("curves_final", ifelse(fm<10, "000", ifelse(fm<100, "00",
                                 ifelse(fm<1000, "0", ""))), fm, ".png")
    writePNG(img, name2)
    print(paste0(fm+1, "/", NFRAMES*NMORPH, ": Writing '", name2, "'..."))
}


#############################
# Create animated files:

# magick -delay 3 -loop 0 curves*.png parametriccurves.gif

# MP4 Video (MPEG-4 AVC/H.264):
# ffmpeg -loop 1 -framerate 24 -i curves_final%04d.png -i deutschland.wav
# -t 412.690 -c:v libx264 -crf 23 -pix_fmt yuv420p parametriccurves.mp4



##_______________________________________________________________________________________ Plotting Functions



## Plot multible "temperature" dependend curves

plot.multicurve.xy <- function (df,main='Main Title', xlab='XLAB', ylab='YLAB', legend.title='LEGEND.TITLE', legend.position="topright", xlim=NA, ylim=NA, range=NA, log10X=FALSE){
# Plot multible curves by 3 columns, where (curveID like temperature, x,y)
# range = often which temperature curves should be plotted
# logX, if TRUE then second column (x-values) are taken as log10
  
  if (is.na(xlim[1])){xlim <- c(min(df[2]), max(df[2]))}
  if (is.na(ylim[1])){ylim <- c(min(df[3]), max(df[3]))}
  
  if (is.na(range) == TRUE){
    curveID <- unique(df[[1]])
  }
  else {					# if a range is given, trim curve ID
    curveID <- df[ df[1] <= max(range),]   
    curveID <- df[ df[1] >= min(range),]
  }
  
  if (log10X == TRUE){df[[2]] <- log10(df[[2]])}

  colors <- c()
  for (i in seq(1, length(curveID),1)){
    colors <- append(colors, rgb(i/length(curveID),0,1-i/length(curveID)))
  }
  
  plot(subset(df, df[1] == min(curveID))[[2]]
      ,subset(df, df[1] == min(curveID))[[3]]
      ,type = 'n'
      ,xlim = xlim,    ylim = ylim
      ,main = main
      ,xlab = xlab
      ,ylab = ylab)

  for (i in seq(1, length(curveID),1)){
    lines(subset(df, df[1] == curveID[i])[[2]]
	  ,subset(df, df[1] == curveID[i])[[3]], type = 'o', pch = i
	  ,col = colors[i])
  }
  grid()
  legend(legend.position, ncol = 2, cex = 0.6, title = legend.title ,legend = curveID, col = colors, lty = 1, pch = seq(1, max(curveID),1))

}

##____________________________________________________________________________________________ Fiber tools



lee.comp.cte <- function (psi,E.m,E.f,nu.m,nu.f,cte.m,rho.m,rho.f,a.11) {

    a.22 <- (1-a.11)/2
    ar <- a.11/a.22

   # calculate lame constants
    lambda.m <- (nu.m * E.m)/((1+nu.m)*(1-2*nu.m))
    lambda.f <- (nu.f * E.f)/((1+nu.f)*(1-2*nu.f))
    mu.m <- G.m <- E.m /(2*(1+nu.m))
    mu.f <- G.f <- E.f /(2*(1+nu.f))

    #convert weight percentage to volume percentage
    phi <- 1/(1+((1-psi)/psi)*(rho.f/rho.m))
    g <- (ar /((ar^2 - 1)^(3/2)))* ((ar*(ar^2 - 1)^(0.5))- acosh(ar))

    #Eshelby Components        

    S.1111 <- (1 /(2*(1-nu.m))) * (1 - 2*nu.m + ((3*ar^2 -1 )/(ar^2-1)) - (1 - 2*nu.m + (3*ar^2)/(ar^2 -1)) * g   )
    S.2222 <- S.3333 <- ((3*ar^2)/( (8*(1-nu.m))*(ar^2 -1))) + (1/(4*(1-nu.m))) * (1-2*nu.m - (9/(4*(ar^2 -1))))*g
    S.2233 <- S.3322 <- (1/(4*(1-nu.m))) * ((ar^2)/(2*(ar^2-1)) - ( 1- 2*nu.m + (3/(4*(ar^2 -1))) )*g  )
    S.2211 <- S.3311 <- ((-1*ar^2)/(2*(1-nu.m)*(ar^2-1))) + (1/(4*(1-nu.m))) * ( ((3*ar^2)/(ar^2 -1)) - (1-2*nu.m) )*g
    S.1122 <- S.1133 <- -1/(2*(1-nu.m)) * (1-2*nu.m + 1/(ar^2 -1)) + (1/(2*(1-nu.m))) * (1-2*nu.m + 3/(2*(ar^2 -1))) * g

    S.2323 <- S.3232 <- 1/(4*(1-nu.m)) * ((ar^2)/(2*(ar^2-1)) + (1-2*nu.m - 3/(4*(ar^2-1))) *g  )

    S.1212 <- S.1313 <- (1/(4*(1-nu.m))) * (1-2*nu.m - (ar^2+1)/(ar^2-1) - 0.5 * (1-2*nu.m - (3*(ar^2 -1))/(ar^2-1)) * g   )

    
    D.1 <- 1+(2*(mu.f-mu.m))/(lambda.f - lambda.m)
    ## !! accoring lee2007 D.2 first factor starts with lambda.f!! 
    D.2 <- (lambda.m + 2* mu.m)/(lambda.f - lambda.m)
    D.3 <- lambda.m / (lambda.f - lambda.m)
    D.4 <- (-2 * (mu.f - mu.m)) / ((lambda.f - lambda.m)) - 3   ## -3 maybe inside bracket?

    B.1 <- phi*D.1 + D.2 +   (1-phi) *     (D.1*S.1111 + 2*S.2211)
    B.2 <- phi + D.3 +       (1-phi) *     (D.1 * S.1122 + S.2222 + S.2233)
    B.3 <- phi + D.3 +       (1-phi) *     (S.1111 + (1 + D.1)*S.2211)
    B.4 <- phi*D.1 + D.2 +   (1-phi) *     (S.1122 + D.1*S.2222 + S.2233)
    B.5 <- phi + D.3 +       (1-phi) *     (S.1122 + S.2222 + D.1 * S.2233)
    
    A.11.t <- (D.4 * (2*B.2 - B.4 - B.5)) / (2*B.2*B.3 - B.1*B.4 - B.1*B.5)
    A.33.t <- (D.4 * (B.3 - B.1)) / (2*B.2*B.3 - B.1*B.4 - B.1*B.5)
    
    a.n.fac <- (1 + phi * A.11.t)
    a.p.fac <- (1 + phi * A.33.t)
    cte.n <-  (1 + phi * A.11.t) * cte.m
    cte.p <-   (1 + phi * A.33.t) * cte.m
    
    return(list(a.n.fac=a.n.fac,a.p.fac=a.p.fac,cte.n=cte.n, cte.p=cte.p, phi=phi))


}






Halpin.Tsai.Micromechanics <- function(psi,E.m,E.f=73000,nu.m,nu.f=0.22,rho.m,rho.f=2600,ar) {
        G.m <- E.m /(2*(1+nu.m)) 
        G.f <- E.f /(2*(1+nu.f)) 

        phi <- 1/(1+((1-psi)/psi)*(rho.f/rho.m))
        df.fac <- data.frame(zeta=rep(NA,4),eta=rep(NA,4))
        
        df.fac$zeta[1] <- 2*ar
        df.fac$zeta[2] <- 2
        df.fac$zeta[3] <- 1
        df.fac$zeta[4] <- (3-4*nu.m+(G.m/G.f))/(4*(1-nu.m))
        df.fac$eta[1] <- ((E.f/E.m)-1)/((E.f/E.m)+df.fac$zeta[1])
        df.fac$eta[2] <- ((E.f/E.m)-1)/((E.f/E.m)+df.fac$zeta[2])
        df.fac$eta[3] <- ((G.f/G.m)-1)/((G.f/G.m)+df.fac$zeta[3])
        df.fac$eta[4] <- ((G.f/G.m)-1)/((G.f/G.m)+df.fac$zeta[4])        

        

        E.11 <- E.m * ((1+(2*ar)*df.fac$eta[1]*phi)/(1-df.fac$eta[1]*phi))
        E.22 <- E.33 <- E.m * ((1+2*df.fac$eta[2]*phi)/(1-df.fac$eta[2]*phi))
        G.12 <- G.13 <- G.m *((1+df.fac$eta[3]*phi)/(1-df.fac$eta[3]*phi))
 
        G.23 <- G.m * ((phi+df.fac$zeta[4] * (1-phi))/((G.m/G.f)*phi+df.fac$zeta[4]*(1-phi)))   ## in excel table first zeta in formula is eta o.O

        nu.12 <- nu.13 <- (E.22/E.11) * (nu.f * phi + nu.m *(1-phi))
        nu.23 <- 0.5 * (E.22/G.23) -1
        #nu.21=nu.31 is missing in stommel2011p52, therefore definition of nu.21 taken from https://de.wikipedia.org/wiki/Transversale_Isotropie (sec. Materialparameter)
        nu.21 <- nu.31 <- (nu.12/E.11)*E.22 

        return(data.frame(psi=psi,phi=phi,E.11=E.11,E.22=E.22,G.12=G.12,G.23=G.23,nu.12=nu.12,nu.23=nu.23,nu.13=nu.13,nu.21=nu.21,nu.31=nu.31))


}


E.m.and.orientation.from.composite <- function(E.t.comp,clte.n.comp,psi,E.m.start=2000,nu.m=0.35,rho.m=1100,cte.m=100e-6,E.f=73000,nu.f=0.22,rho.f=2660,cte.f=10e-6,ar=20,max.iter=40,micromechmodel='tandon.weng') {
        require('minpack.lm')
        ## (190915) Function to extract Matrix Modulus E.m from given composite Modulus E.t.comp of short fiber filled thermoplastics
        ## Requires functions: advani.tucker.tensile.homogenization(), tandon.weng()
        # E.t.comp: Measured Modulus in flow direction of short fiber filled thermoplastic [MPa]
        # psi: Fiber weight ratio [kg/kg]
        # E.m.limits: Vector with min,max Values of search region of possible matrix moduli [MPa]
        # nu.m: Matrix Poisson Ratio [-]
        # rho.m: Matrix Density [kg/m^3]
        # E.f: Fiber Modulus [MPa]
        # nu.f: Fiber Poisson Ratio [-]
        # rho.f: Fiber Density [kg/m^3]
        # ar: Fiber aspect ratio [mm/mm] "l/d"
        ## a.11, a22: first and second Eigenvalue of orientation tensor of tensile bar a.ij [0-1] [-]
        # a.33: Is calculated from a.11 and a.22
        # micromechmodel: either tandon.weng or halpin.tsai approach
        ##Out: E.m: Matrix modulus [MPa]
        
        
        #E.11.composite [MPa]:  Input value for fit, "which E.m is required to get this E.11 of composite"
        #psi.comp [kg/kg]:      Fiber weight ratio of composite
        #ar [mm/mm]:            fiber aspect ratio min-max-limit vector for fit 10-30, default 20


                E.m <- E.m.resi<- E.m.start
                a.11 <- a.11.resi <- 0.4
                iter <- 0
                #psi <- 0.35
                #E.t.comp <- 11500
                #clte.n.comp <- 15e-6
                
                if (micromechmodel == "tandon.weng") {
                        while ((abs(E.m.resi) > 10 | abs(a.11.resi) > 1e-7) & iter < max.iter){ 
                            iter <- iter+1
                            resiFun.alpha <- function(a.11,E.m){
                                                        C.t.ho <- advani.tucker.tensile.homogenization(a.11,(1-a.11)*0.5,
                                                                    E.11=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.11,
                                                                    E.22=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.22,
                                                                    G.12=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.12,
                                                                    G.23=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.23,
                                                                    nu.12=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.12,
                                                                    nu.23=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.23,
                                                                    nu.21=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.21)$C.t.ho

                                                                    #return(clte.n.comp - Rosenberg.Hashin(psi,E.m,nu.m,rho.m,cte.m,E.f,nu.f,rho.f,cte.f,C.t.ho)$clte.n)
                                                                    return (clte.n.comp - lee.comp.cte(psi=psi,E.m=E.m, E.f=E.f,nu.m=nu.m,nu.f=nu.f,cte.m=cte.m,rho.m=rho.m, rho.f=rho.f,a.11=a.11)$cte.p                                                                                     
                                                                        )

                                                                    

                            }
                            
                            alpha.fit <- nls.lm(par=a.11, fn=resiFun.alpha ,E.m=E.m,control = nls.lm.control(nprint=0), upper=c(0.9), lower=c(0.3333))
                            a.11.1 <- alpha.fit$par
                            
                            
                            resiFun <- function(E.m,modul,a.11) {
                                                                            modul - advani.tucker.tensile.homogenization(a.11,(1-a.11)*0.5,
                                                                            E.11=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.11,
                                                                            E.22=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.22,
                                                                            G.12=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.12,
                                                                            G.23=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.23,
                                                                            nu.12=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.12,
                                                                            nu.23=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.23,
                                                                            nu.21=tandon.weng(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.21)$E.ho.n
                            }   
                            tandon.weng.fit <- nls.lm(par=E.m, fn=resiFun ,modul=E.t.comp, ,a.11=a.11.1,control = nls.lm.control(nprint=0))

                            E.m.resi <- E.m - tandon.weng.fit$par
                            a.11.resi <- a.11.1 - nls.lm(par=a.11, fn=resiFun.alpha ,E.m=tandon.weng.fit$par,control = nls.lm.control(nprint=0), upper=c(0.9), lower=c(0.333))$par
                            
                            E.m <- tandon.weng.fit$par
                            a.11 <- nls.lm(par=a.11, fn=resiFun.alpha ,E.m=tandon.weng.fit$par,control = nls.lm.control(nprint=0), upper=c(0.9999), lower=c(0.3333))$par
                    }                                                                                                                                                         
                }
                
                if (micromechmodel == "halpin.tsai") {
                        while ((abs(E.m.resi) > 10 | abs(a.11.resi) > 1e-7) & iter < max.iter){ 
                            iter <- iter+1
                            resiFun.alpha <- function(a.11,E.m){
                                                        C.t.ho <- advani.tucker.tensile.homogenization(a.11,(1-a.11)*0.5,
                                                                    E.11=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.11,
                                                                    E.22=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.22,
                                                                    G.12=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.12,
                                                                    G.23=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.23,
                                                                    nu.12=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.12,
                                                                    nu.23=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.23,
                                                                    nu.21=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.21)$C.t.ho

                                                                    #return(clte.n.comp - Rosenberg.Hashin(psi,E.m,nu.m,rho.m,cte.m,E.f,nu.f,rho.f,cte.f,C.t.ho)$clte.n)
                                                                    return (clte.n.comp - lee.comp.cte(psi=psi,E.m=E.m, E.f=E.f,nu.m=nu.m,nu.f=nu.f,cte.m=cte.m,rho.m=rho.m, rho.f=rho.f,a.11=a.11)$cte.p                                                                                    
                                                                        )

                                                                    

                            }
                            
                            alpha.fit <- nls.lm(par=a.11, fn=resiFun.alpha ,E.m=E.m,control = nls.lm.control(nprint=0), upper=c(0.9), lower=c(0.3333))
                            a.11.1 <- alpha.fit$par
                            
                            
                            resiFun <- function(E.m,modul,a.11) {
                                                                            modul - advani.tucker.tensile.homogenization(a.11,(1-a.11)*0.5,
                                                                            E.11=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.11,
                                                                            E.22=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.22,
                                                                            G.12=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.12,
                                                                            G.23=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.23,
                                                                            nu.12=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.12,
                                                                            nu.23=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.23,
                                                                            nu.21=Halpin.Tsai.Micromechanics(psi,E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.21)$E.ho.n
                            }   
                            halpin.tsai.fit <- nls.lm(par=E.m, fn=resiFun ,modul=E.t.comp, ,a.11=a.11.1,control = nls.lm.control(nprint=0))

                            E.m.resi <- E.m - halpin.tsai.fit$par
                            a.11.resi <- a.11.1 - nls.lm(par=a.11, fn=resiFun.alpha ,E.m=halpin.tsai.fit$par,control = nls.lm.control(nprint=0), upper=c(0.9), lower=c(0.333))$par
                            
                            E.m <- halpin.tsai.fit$par
                            a.11 <- nls.lm(par=a.11, fn=resiFun.alpha ,E.m=halpin.tsai.fit$par,control = nls.lm.control(nprint=0), upper=c(0.9999), lower=c(0.3333))$par
                    }                                                                                                                                                         
                }
                
                return(list(E.m=E.m,a.11=a.11,iter=iter))
                  
}

# example Ultramid A3WG7 GF35:
#E.m.and.orientation.from.composite(11500,17.5e-6,0.35,ar=20,E.m.start=2000,cte.m=40e-6,max.iter=40,rho.m=1160,nu.m=0.25,micromechmodel='tandon.weng')





isotropic.compliance.matrix <- function(E,nu) {
    ## calculate isotropic compliance matrix from E,nu for isotropic materials (190918)

    out <- matrix(nrow=6,ncol=6)
    out[1,] <- c(1,-nu,-nu,0,0,0)
    out[2,] <- c(-nu,1,-nu,0,0,0)
    out[3,] <- c(-nu,-nu,1,0,0,0)
    out[4,] <- c(0,0,0,2*(1+nu),0,0)
    out[5,] <- c(0,0,0,0,2*(1+nu),0)
    out[6,] <- c(0,0,0,0,0,2*(1+nu))

    return(out / E)

}

Rosenberg.Hashin <- function(psi,E.m,nu.m,rho.m,cte.m,E.f,nu.f,rho.f,cte.f,C.t.ho) {
    ## Function to calculate CLTE (Coefficient of linear thermal expansion) for short fiber filled composites according kennedy2013p126 (190922)
    ## requires isotropic.compliance.matrix()
    # C.t.ho:  3x3 Homogenized Stiffness matrix for C_ij i,j=1-3 from function advani.tucker.tensile.homogenization()
    # psi:      Fiber weight percentage [kg/kg]
    # cte.m     Matrix CLTE [1/K]
    # cte.f     Fiber CLTE [1/K]    
    # E.m:      Matrix Modulus [MPa]
    # nu.m: Matrix Poisson Ratio [-]
    # rho.m: Matrix Density [kg/m^3]
    # E.f: Fiber Modulus [MPa]
    # nu.f: Fiber Poisson Ratio [-]
    # rho.f: Fiber Density [kg/m^3]
    
    phi.f <- 1/(1+((1-psi)/psi)*(rho.f/rho.m))
    phi.m <- 1-phi.f

    S.m <- isotropic.compliance.matrix(E.m,nu.m)[1:3,1:3]
    S.f <- isotropic.compliance.matrix(E.f,nu.f)[1:3,1:3]
    #S.comp <-solve(C.t.ho[1:3,1:3]) 
    S.comp <- C.t.ho[1:3,1:3]   ## from ar approach S.not to be converted


    cte.comp <- phi.f * cte.f*diag(3)+ phi.m *cte.m*diag(3)  + ((S.comp*diag(3) - phi.f * S.f*diag(3) - phi.m * S.m*diag(3))*solve(S.f*diag(3) - S.m*diag(3))) * (cte.f*diag(3) - cte.m*diag(3))
    #cte.comp <- phi.f * cte.f*matrix(1,nrow=6,ncol=1)+ phi.m *cte.m *matrix(1,nrow=6,ncol=1) + (((S.comp - (phi.f * S.f) - (phi.m * S.m))%*%solve(S.f - S.m))) %*%( (cte.f - cte.m)*matrix(1,nrow=6,ncol=1))

    #cte.comp <- cte.f *diag(3) + ((S.f*diag(3) - S.comp*diag(3))%*%(solve(S.m - S.f))) * (cte.f*diag(3)-cte.m*diag(3))
    #cte.n <- cte.f  + ((S.f[1,1] - S.comp[1,1])/(S.m[1,1] - S.f[1,1])) * (cte.f-cte.m)
    #cte.p <- cte.f  + ((S.f[2,2] - S.comp[2,2])/(S.m[2,2] - S.f[2,2])) * (cte.f-cte.m)

    return(
    list(clte.n=sum(cte.comp[1,]),clte.p=sum(cte.comp[2,]))
    #list(clte.n=cte.n,clte.p=cte.p)
    )
}












Rosenberg.Hashin <- function(psi,E.m,nu.m,rho.m,cte.m,E.f,nu.f,rho.f,cte.f,C.t.ho) {
    ## Function to calculate CLTE (Coefficient of linear thermal expansion) for short fiber filled composites according kennedy2013p126 (190922)
    ## requires isotropic.compliance.matrix()
    # C.t.ho:  3x3 Homogenized Stiffness matrix for C_ij i,j=1-3 from function advani.tucker.tensile.homogenization()
    # psi:      Fiber weight percentage [kg/kg]
    # cte.m     Matrix CLTE [1/K]
    # cte.f     Fiber CLTE [1/K]    
    # E.m:      Matrix Modulus [MPa]
    # nu.m: Matrix Poisson Ratio [-]
    # rho.m: Matrix Density [kg/m^3]
    # E.f: Fiber Modulus [MPa]
    # nu.f: Fiber Poisson Ratio [-]
    # rho.f: Fiber Density [kg/m^3]
    
    phi.f <- 1/(1+((1-psi)/psi)*(rho.f/rho.m))
    phi.m <- 1-phi.f

    K.m <- E.m /(1-2*nu.m)
    K.f <- E.f /(1-2*nu.f)
    
    
    #############  >((
    S.m <- isotropic.compliance.matrix(E.m,nu.m)[1:3,1:3]
    S.f <- isotropic.compliance.matrix(E.f,nu.f)[1:3,1:3]
    solve(isotropic.stiffness.matrix(E.m,nu.m))
    
    #S.comp <- S.ortho.ar[1:3,1:3]
    S.comp <- solve(C.t.ho)  
    
    cte.f <- cte.f*diag(3)
    cte.m <- cte.m*diag(3)
    
    #phi.f * cte.f + phi.m *cte.m + (cte.f - cte.m)/(1/K.f-1/K.m)* (3/C.t.ho[1,1] - ((phi.f/K.f + phi.m/K.m)))
    ## isotropic components approach
    #cte.comp <- phi.f * cte.f + phi.m *cte.m + (cte.f - cte.m)/(1/K.f-1/K.m) * (3*((S.comp)) - (phi.f/K.f + phi.m/K.m))
  

    #(solve(S.f + S.m)) %*% (S.f + S.m)
    
    #cte.comp <- phi.f * cte.f*diag(6)+ phi.m *cte.m*diag(6)  + (solve(S.f + S.m) %*% (3*(S.comp) - (nu.f * S.f + nu.m * S.m))  ) %*% (cte.f*diag(6) - cte.m*diag(6))

    cte.f <- rbind(cte.f*matrix(1,nrow=3,ncol=1))#,matrix(0,nrow=3,ncol=1))
    cte.m <- rbind(cte.m*matrix(1,nrow=3,ncol=1))#,matrix(0,nrow=3,ncol=1))

    cte.comp <- phi.f * cte.f+ phi.m *cte.m + (((S.comp - (phi.f * S.f) + (phi.m * S.m))%*%solve(S.f - S.m))) %*%( cte.f - cte.m)
    
    ##Rosenberg Hashin
    cte.comp <- phi.f * cte.f+ phi.m *cte.m + ((S.comp - phi.f * S.f - phi.m * S.m))%*%solve(S.f - S.m) %*%( cte.f - cte.m)


    ## Benveniste Dovrak
    #cte.comp <- cte.f + ((S.f - S.comp)%*%solve(S.m-S.f)) %*%( cte.f - cte.m)
    
    
    #cte.comp <- cte.f *diag(3) + ((S.f*diag(3) - S.comp*diag(3))%*%(solve(S.m - S.f))) * (cte.f*diag(3)-cte.m*diag(3))
    #cte.n <- cte.f  + ((S.f[1,1] - S.comp[1,1])/(S.m[1,1] - S.f[1,1])) * (cte.f-cte.m)
    #cte.p <- cte.f  + ((S.f[2,2] - S.comp[2,2])/(S.m[2,2] - S.f[2,2])) * (cte.f-cte.m)

    return(
    list(clte.n=sum(cte.comp[1,]),clte.p=sum(cte.comp[2,]))
    #list(clte.n=cte.n,clte.p=cte.p)
    )
}


schapery.cte.m.from.comp.mixing <- function(psi ,E.m, E.f, nu.m, nu.f,cte.n,cte.p,rho.m,rho.f) {
    ## Function to extract matrix cte from measured composite values by simple mixing rule  (191124)

    phi.f <- phi <- 1/(1+((1-psi)/psi)*(rho.f/rho.m))
    phi.m <- 1-phi.f

    #cte.m from formula of cte.n 
    cte.m.n <- (cte.n * ( E.m *(1-phi) + E.f*phi) - E.f*cte.f*phi) / (E.m*(1-phi))
    #cte.m from formula of cte.p
    cte.m.p <- (cte.p - (1+phi)*cte.f*phi - cte.n * nu.f *phi) / ((1+nu.m) * (1-phi) - cte.n *(1-phi))
    
    cte.n.fit<-cte.n
    resifun <- function(cte.n.fit) {((cte.n.fit * ( E.m *(1-phi) + E.f*phi) - E.f*cte.f*phi) / (E.m*(1-phi))) - 
                                       ((cte.p - (1+phi)*cte.f*phi - cte.n.fit * nu.f *phi) / ((1+nu.m) * (1-phi) - cte.n.fit *(1-phi))) }
                                       
    cte.n.lm <- nls.lm(par=cte.n.fit, fn=resifun ,control = nls.lm.control(nprint=0))    
    cte.n.corr <- cte.n.lm$par
    
    cte.n.calc <- (E.m*cte.m.p*(1-phi)+E.f*cte.f*phi)/(E.m*(1+phi)+E.f*phi)
    cte.p.calc <- (1+nu.m)*cte.m.p*(1-phi)+ (1+nu.f)*cte.f*phi-cte.n.calc*(cte.m.p*(1-phi)+nu.f*phi)
   
    return(list(cte.m.n=cte.m.n,cte.m.p=cte.m.p,cte.n.corr=cte.n.corr,cte.n.calc=cte.n.calc,cte.p.calc=cte.p.calc))

}




shapery.composite.cte <- function(E.comp,cte.m,cte.f,E.m,nu.m,E.f,nu.f,rho.f,rho.m,psi){
    ## Not used. try     
    phi.f <- 1/(1+((1-psi)/psi)*(rho.f/rho.m))
    phi.m <- 1-phi.f


    K.m <- E.m /(1-2*nu.m)
    K.f <- E.f /(1-2*nu.f)
    G.m <- G.from.E(E.m,nu.m)
    G.f <- G.from.E(E.f,nu.f)


    cte.av <- phi.f * cte.f + phi.m *cte.m 
    K.av <- phi.f * K.f + phi.m *K.m 
    G.av <- phi.f * G.f + phi.m *G.m 
    K.cte.av <- phi.f * K.f *cte.f + phi.m *K.m * cte.m

    K.cte.2.av <- phi.f * K.f *((cte.m^2+cte.f^2)) + phi.m *K.m * (cte.m^2+cte.f^2)

    #E.j <- C.t.ho[1,1] + C.t.ho[2,2]+ C.t.ho[3,3]

    E.u <- (1/3) * ((1/G.av) + (1/(3*K.av)))
    E.L <- phi.f * 1/E.f + phi.m *1/E.m 


    delta.alpha <- (((3*(1/E.comp)- E.u)^0.5 * (E.L- 1/E.comp)^0.5)/ (E.L - E.u) ) * 
                                                                    ((K.cte.2.av- K.cte.av^2) * (E.L-E.u) - 1/9 * (cte.av - (K.cte.av/K.av)) ^2)^0.5
    delta.alpha <- 0
    cte <- cte.av + ((K.cte.av/K.av) - cte.av) * ((E.L - (1/E.comp))/(E.L - E.u)) 

    (1+nu.f)*cte.f*phi.f + (1+nu.m)*cte.m*phi.m - (phi.f*nu.f+phi.m*nu.m)*41e-6

    return(cte)

}



tandon.weng <- function(psi,E.m,E.f=73000,nu.m,nu.f=0.22,rho.m,rho.f=2600,ar){
    ## Function to calculate stiffness components for unidirectional symmetrical elepsoid compound
    ## Transverse Isotropic for full orientation #stommel2011p49
    #psi: Fiber weight percentage [kg/kg]
    #E.m: Matrix Modulus [MPa]
    #E.f: Fiber Modulus [MPa]
    #nu.m, nu.f: Poisson Ratio for Matrix and Fiber [-]
    #rho.m, rho.f: Density for Matrix and Fiber [kg/m^3]
    #ar: Fiber aspect ratio
    
    # calculate lame constants
    lambda.m <- (nu.m * E.m)/((1+nu.m)*(1-2*nu.m))
    lambda.f <- (nu.f * E.f)/((1+nu.f)*(1-2*nu.f))
    mu.m <- G.m <- E.m /(2*(1+nu.m))
    mu.f <- G.f <- E.f /(2*(1+nu.f))

    #convert weight percentage to volume percentage
    phi <- 1/(1+((1-psi)/psi)*(rho.f/rho.m))
    df <- data.frame(phi=phi,psi=psi)
    
    g <- (ar /((ar^2 - 1)^(3/2)))* ((ar*(ar^2 - 1)^(0.5))- acosh(ar))

    #Eshelby Components        
    if(ar > 1) {
    S.1111 <- (1 /(2*(1-nu.m))) * (1 - 2*nu.m + ((3*ar^2 -1 )/(ar^2-1)) - (1 - 2*nu.m + (3*ar^2)/(ar^2 -1)) * g   )
    S.2222 <- S.3333 <- ((3*ar^2)/( (8*(1-nu.m))*(ar^2 -1))) + (1/(4*(1-nu.m))) * (1-2*nu.m - (9/(4*(ar^2 -1))))*g
    S.2233 <- S.3322 <- (1/(4*(1-nu.m))) * ((ar^2)/(2*(ar^2-1)) - ( 1- 2*nu.m + (3/(4*(ar^2 -1))) )*g  )
    S.2211 <- S.3311 <- ((-1*ar^2)/(2*(1-nu.m)*(ar^2-1))) + (1/(4*(1-nu.m))) * ( ((3*ar^2)/(ar^2 -1)) - (1-2*nu.m) )*g
    S.1122 <- S.1133 <- -1/(2*(1-nu.m)) * (1-2*nu.m + 1/(ar^2 -1)) + (1/(2*(1-nu.m))) * (1-2*nu.m + 3/(2*(ar^2 -1))) * g

    S.2323 <- S.3232 <- 1/(4*(1-nu.m)) * ((ar^2)/(2*(ar^2-1)) + (1-2*nu.m - 3/(4*(ar^2-1))) *g  )

    S.1212 <- S.1313 <- (1/(4*(1-nu.m))) * (1-2*nu.m - (ar^2+1)/(ar^2-1) - 0.5 * (1-2*nu.m - (3*(ar^2 -1))/(ar^2-1)) * g   )
    }
    #Eshelby Components for spheres
    if (ar == 1) {
        S.1111 <- S.2222 <- S.3333 <- (7-5*nu.m)/(15*(1-nu.m))
        S.1122 <- S.2233 <- S.3311 <- (5*nu.m - 1) / (15*(1-nu.m))
        S.1212 <- S.2323 <- S.3131 <- (4-5*nu.m)/ (15*(1-nu.m))
    }
    
    D.1 <- 1+(2*(mu.f-mu.m))/(lambda.f - lambda.m)
    ## !! accoring lee2007 D.2 first factor starts with lambda.f!! 
    D.2 <- (lambda.m + 2* mu.m)/(lambda.f - lambda.m)
    D.3 <- lambda.m / (lambda.f - lambda.m)
    #D.4 <- (-2 * (mu.f - mu.m)) / ((lambda.f - lambda.m)) - 3   ## -3 maybe inside bracket?

    df$B.1 <- phi*D.1 + D.2 +   (1-phi) *     (D.1*S.1111 + 2*S.2211)
    df$B.2 <- phi + D.3 +       (1-phi) *     (D.1 * S.1122 + S.2222 + S.2233)
    df$B.3 <- phi + D.3 +       (1-phi) *     (S.1111 + (1 + D.1)*S.2211)
    df$B.4 <- phi*D.1 + D.2 +   (1-phi) *     (S.1122 + D.1*S.2222 + S.2233)
    df$B.5 <- phi + D.3 +       (1-phi) *     (S.1122 + S.2222 + D.1 * S.2233)

    df$A.1 <- D.1 * (df$B.4+df$B.5) - 2* df$B.2
    df$A.2 <- (1+D.1)*df$B.2 - (df$B.4+df$B.5)
    df$A.3 <- df$B.1 - D.1 * df$B.3
    df$A.4 <- (1+D.1) * df$B.1 - 2*df$B.3
    df$A.5 <- (1-D.1)/(df$B.4-df$B.5)
    df$A.6 <- 2*df$B.2*df$B.3 - df$B.1*(df$B.4 + df$B.5)


    df$E.11 <- E.m/(1+((df$phi * (df$A.1+2*nu.m*df$A.2))/df$A.6))
    df$E.22 <- df$E.33 <- E.m/(1+(((df$phi*(2*nu.m*df$A.3+(1-nu.m)*df$A.4+(1+nu.m)*df$A.5*df$A.6)))/(2*df$A.6)))
    df$G.12 <- G.13 <- G.m * (1+(df$phi/((G.m/(G.f-G.m))+2*(1-df$phi)*S.1212)))
    df$G.23 <-  G.m * (1+(df$phi/((G.m/(G.f-G.m))+2*(1-df$phi)*S.2323)))
    df$nu.12 <- df$nu.13 <- (nu.m * df$A.6-df$phi*(df$A.3-nu.m*df$A.4))/(df$A.6+df$phi*(df$A.1+2*nu.m*df$A.2))
    
    #df$nu.23 <- (df$E.22/(2*df$G.23))-1
    K.m <- E.m/(3-6*nu.m)
    K.f <- E.m/(3-6*nu.f)    
    
    
    #from here modifications from moldflow help http://help.autodesk.com/view/MFIA/2016/ENU/
    
 
    
    zeta <- (K.m/G.m)/((K.m/G.m)+2)
    eta.3 <- ((G.f/G.m)-1) / ((G.f/G.m)+1)
    eta.4 <- ((G.f/G.m)-1) / ((G.f/G.m)+zeta)
    
    df$G.12 <- G.m * (1+ phi*eta.3)/(1- phi*eta.3)
    df$G.23 <- G.m * (1+ zeta*phi*eta.4)/(1- phi*eta.4)
    
    df$nu.12 <- nu.f*phi + (1-phi)*nu.m
    
    df$nu.23 <- (zeta*(3-4*(df$nu.12)^2)-1)/(zeta+1)
    ##nu.21=nu.31 is missing in stommel2011p52, therefore definition of nu.21 taken from https://de.wikipedia.org/wiki/Transversale_Isotropie (sec. Materialparameter)
    df$nu.21 <- df$nu.31 <- (df$nu.12/df$E.11)*df$E.22
    
    

    
    #Composite Bulk modulus from lusti2003p14 
    df$K.23 <- K.m * ((1+nu.m)*(1-2*nu.m))/(1-nu.m*(1+2*df$nu.21)+phi*(2*(df$nu.21-nu.m)*df$A.3 + (1-nu.m*(1+2*df$nu.21))*df$A.4)/(df$A.6))
    
    #Composite Bulk modulus from Halping Tsai Paper affdl1976p5
    #zeta.b <- (4*G.m)/(3*K.m)  #bulk modulus zeta
    #M.R <- K.f/K.m
    #eta.b <- (M.R-1)/(M.R+zeta.b)
    #df$K.23 <- (1 + zeta.b * eta.b * phi)/(1- eta.b *phi) *K.m
    
    #df$nu.23 <- 1 - df$nu.21 - (df$E.22/(2*df$K.23)) + (df$E.22*(1-df$nu.12))/(2*df$E.11)

    df$nu.21 <- df$nu.12

    
    return(data.frame(psi=df$psi,phi=df$phi,E.11=df$E.11,E.22=df$E.22,G.12=df$G.12,G.23=df$G.23,nu.12=df$nu.12,nu.23=df$nu.23,nu.13=df$nu.13,nu.21=df$nu.31,nu.31=df$nu.31))
}



advani.tucker.tensile.homogenization <- function(a.11,a.22,E.11,E.22,G.12,G.23,nu.12,nu.23,nu.21){
    ## Fiber homogenization for given Orientation eigenvalues a.ii and Stiffness values of unidirectional values
    ## according #stommel2011p52 and #han2002  (190915)
    # Stniffness and poisson values can be taken from function tandon.weng() or Halpin.Tsai.Micromechanics()
    # a.11, a22: first and second Eigenvalue of orientation tensor a.ij [0-1] [-]
    # a.33: Is calculated from a.11 and a.22
    # E.11: Stiffness entry for in fiber direction [MPa]
    # E.22: Stiffness for perpendicular to fiber direction [MPa]
    # G.12, G.23: Shear Moduli [MPa]
    # nu.xx : Poisson Ratios [-]
    ##Out: C.t.ho=Stniffness value in fiber direction,
    ##     E.ho.n=Youngs Modulus in Fiber direction
    ##     E.ho.p=Youngs Modulus transverse Fiber direction
    
    
    a.33 <- 1-(a.11+a.22)

    # Assembly Transversal istoropic C Stiffness matric according https://de.wikipedia.org/wiki/Transversale_Isotropie
    lambda <- ((nu.12*nu.21+nu.23)/((1-nu.23-2*nu.12*nu.21)*(1+nu.23)))*E.22
    C.ti <- matrix(nrow=6,ncol=6) #transversal isotropic
    C.ti[1,1] <- ((1-nu.23)/(1-nu.23-2*nu.12*nu.21))*E.11
    C.ti[2,2] <- C.ti[3,3] <- lambda + 2*G.23
    C.ti[2,3] <- lambda
    C.ti[4,4] <- G.23
    C.ti[5,5] <- C.ti[6,6] <- G.12
    C.ti[1,2] <- C.ti[1,3] <- 2*nu.12*(lambda+G.23)
    C.ti[1:3,4:6] <- C.ti[4,5:6] <- C.ti[5,6] <- 0      #create zero entries
    C.ti[lower.tri(C.ti)] <- t(C.ti)[lower.tri(C.ti)]  # make symmetric

    # Advani Tucker B Values #stommel2011p52
    B <- c( C.ti[1,1] + C.ti[2,2] - 2*C.ti[1,2] - 4*C.ti[4,4] ,
            C.ti[1,2] - C.ti[2,3],
            C.ti[4,4] + 0.5*(C.ti[2,3]-C.ti[2,2]),
            C.ti[2,3],
            0.5* (C.ti[2,2]-C.ti[2,3])
     )
     
     #forth-order colsure by #hahn2002p368
     r <- a.22/(a.11+3*a.22)
     s <- a.33/(1+2*a.33)
     b.1111 <- (1-r)*(1-s)*a.11
     b.2222 <- (1-s)*a.22-r*(1-s)*a.11
     b.3333 <- 3*s*a.33
     
     # homogenized 3x3 Tensile Stiffness matrix
        C.t.ho <- matrix(nrow=3,ncol=3) 
        C.t.ho[1,1] <- B[1]*b.1111+B[2]*(a.11+a.11)+B[3]*(a.11+a.11+a.11+a.11)+B[4]+B[5]*2
        C.t.ho[2,2] <- C.t.ho[3,3] <- B[1]*b.2222+B[2]*(a.22+a.22)+B[3]*(a.22+a.22+a.22+a.22)+B[4]+B[5]*2
        C.t.ho[2,3] <- B[1]*((a.22-b.2222)/2)+B[2]*(a.22+a.33)+B[3]*(0+0+0+0)+B[4]+B[5]*0
        C.t.ho[1,2] <- C.t.ho[1,3] <- B[1]*((a.11-b.1111)/2)+B[2]*(a.11+a.22)+B[3]*(0+0+0+0)+B[4]+B[5]*0

        C.t.ho[lower.tri(C.t.ho)] <- t(C.t.ho)[lower.tri(C.t.ho)]  # make symmetric
        
     # Measured tensile modulus in fiber direction
        #tensile.strain <- 0.00276                                                                  #example tensile strain value to calculate stress from C matrix
        #strain.tensor.n <- c(tensile.strain,-1*tensile.strain*nu.12,-1*tensile.strain*nu.12)     #calculate strain tensor from nu 
        #strain.tensor.p <- c(-1*tensile.strain*nu.21,tensile.strain,-1*tensile.strain*nu.21)     #calculate strain tensor from nu 

        #stress.n <- strain.tensor.n * c(C.t.ho[1,1],C.t.ho[1,2],C.t.ho[1,3])
        #stress.p <- strain.tensor.p * c(C.t.ho[2,1],C.t.ho[2,2],C.t.ho[2,3])
        #E.ho.n <- sum(stress.n)/tensile.strain
        #E.ho.p <- sum(stress.p)/tensile.strain
        
        return(list(C.t.ho=C.t.ho,E.ho.n=C.t.ho[1,1],E.ho.p=C.t.ho[2,2],C.ti=C.ti))
        
 }       




E.m.from.composite <- function(E.t.comp,psi,E.m.limits=c(500,5000),nu.m,rho.m,E.f=73000,nu.f=0.22,rho.f=2660,ar=20,a.11=0.8,a.22=0.1,micromechmodel='tandon.weng') {
        ## (190915) Function to extract Matrix Modulus E.m from given composite Modulus E.t.comp of short fiber filled thermoplastics
        ## Requires functions: advani.tucker.tensile.homogenization(), tandon.weng()
        # E.t.comp: Measured Modulus in flow direction of short fiber filled thermoplastic [MPa]
        # psi: Fiber weight ratio [kg/kg]
        # E.m.limits: Vector with min,max Values of search region of possible matrix moduli [MPa]
        # nu.m: Matrix Poisson Ratio [-]
        # rho.m: Matrix Density [kg/m^3]
        # E.f: Fiber Modulus [MPa]
        # nu.f: Fiber Poisson Ratio [-]
        # rho.f: Fiber Density [kg/m^3]
        # ar: Fiber aspect ratio [mm/mm] "l/d"
        ## a.11, a22: first and second Eigenvalue of orientation tensor of tensile bar a.ij [0-1] [-]
        # a.33: Is calculated from a.11 and a.22
        # micromechmodel: either Tandon.Weng or Halpin-Tsai approach
        ##Out: E.m: Matrix modulus [MPa]
        
        
        #E.11.composite [MPa]:  Input value for fit, "which E.m is required to get this E.11 of composite"
        #psi.comp [kg/kg]:      Fiber weight ratio of composite
        #ar [mm/mm]:            fiber aspect ratio min-max-limit vector for fit 10-30, default 20
        require('minpack.lm')
        
        if (micromechmodel == "tandon.weng") {
            resiFun <- function(para,modul,psi,a.11,a.22) {
                                                            modul - advani.tucker.tensile.homogenization(a.11,a.22,
                                                            E.11=tandon.weng(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.11,
                                                            E.22=tandon.weng(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.22,
                                                            G.12=tandon.weng(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.12,
                                                            G.23=tandon.weng(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.23,
                                                            nu.12=tandon.weng(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.12,
                                                            nu.23=tandon.weng(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.23,
                                                            nu.21=tandon.weng(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.21)$E.ho.n
            }  
        }
        if (micromechmodel == "halpin.tsai") {
            resiFun <- function(para,modul,psi,a.11,a.22) {
                                                            modul - advani.tucker.tensile.homogenization(a.11,a.22,
                                                            E.11=Halpin.Tsai.Micromechanics(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.11,
                                                            E.22=Halpin.Tsai.Micromechanics(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$E.22,
                                                            G.12=Halpin.Tsai.Micromechanics(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.12,
                                                            G.23=Halpin.Tsai.Micromechanics(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$G.23,
                                                            nu.12=Halpin.Tsai.Micromechanics(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.12,
                                                            nu.23=Halpin.Tsai.Micromechanics(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.23,
                                                            nu.21=Halpin.Tsai.Micromechanics(psi,para$E.m,E.f=E.f,nu.m=nu.m,nu.f=nu.f,rho.m=rho.m,rho.f=rho.f,ar=ar)$nu.21)$E.ho.n
            }  
        }
        para <-                 list(E.m=mean(E.m.limits))
        tandon.weng.fit <- nls.lm(par=para, fn=resiFun ,modul=E.t.comp, psi=psi,a.11=a.11,a.22=a.22,control = nls.lm.control(nprint=1), 
                                upper=c(max(E.m.limits)),
                                lower=c(min(E.m.limits)))

        return(tandon.weng.fit$par$E.m)

}















Phi <- function(Psi=0.3, rho.f = 2550, rho.m=1000){
## R script for calculating Fiber volume fraction Phi from fiber weight fraction Psi [0,1]
# fiber density rho.f: Default 2550 for glas (! according data tables 6500)
# polymer density rho.m

Psi / (Psi + ((1 - Psi) * (rho.f)/(rho.m)))

}


apprx.e.matrix.reinforced <- function(E.11, E.f=70000,  psi=seq(0,.5,0.1), ld=15, rho.m=1120, rho.f=2550){
## Approximation of matrix modulus of a short fiber filled thermoplastic
# E11: Measured Tensile Modulus [MPa], psi: fiber weight fraction [0,1], ld: length/diameter ratio [-]
# rho.m: matrix density, rho.f: fiber density
phi <- psi / (psi + (1 - psi) * (rho.f/rho.m))  # calc of fiber volume fraction phi [-]
print('Fiber weight ratio:')
print(psi)
zeta <- 20 * 2
B1 <- E.f * (1 + zeta * phi)
B2 <- E.11^2 * zeta^2 + 2 * E.11^2 * zeta * phi + 2 * E.11 * E.f * zeta - 2 * E.11 * E.f * zeta^2 * phi
B3 <- E.11^2 * phi^2 - 2 * E.11 * E.f  * phi + 2 * E.11 * E.f * zeta * phi^2 + E.f^2 +2* E.f^2 * zeta * phi + E.f^2 * zeta^2 * phi^2 - 8 * E.11 * zeta * phi

E.m <- 0.5 * ((B1 - (sqrt(B2 + B3)))/(zeta*(phi-1)))
print('Matrix Modulus:')
print(E.m)
}

tsai.wu.transversal.isotropic <- function(sigma.allw, m=1.3,n.p.ratio=0.5 ){
library(plot3D)
scale = 2

# R.t tension strenght
R.t.p <- sigma.allw              # p parallel
R.t.n <- R.t.p * n.p.ratio             ## get normal to fiber by reverse engineered matrix

# R compression compression strenght
R.c.p <-  R.t.p * m
R.c.n <- R.t.n * m

# R shear
R.s.pn <- 0.577 * R.t.p          # according stommel2011p175

A1 <- (1 / R.t.p) - (1 / R.c.p)             # A1  according ehrenstein2011p256 
A2 <- (1 / R.t.n) - (1 / R.c.n)             # A2
A3 <- (1/ (R.t.p * R.c.p))                # A3
A4 <- (1/ (R.t.n * R.c.n))                # A4
A5 <- 1 / (R.s.pn^2)                 
A6 <-0/ (sqrt(R.t.p * R.c.p * R.t.n * R.c.n))

F1 <- paste(A1, '* sx ', sep = ' ')
F2 <- paste(A2, '* sqrt(sy^2 + sz^2) ', sep = ' ')
F3 <- paste(A3, '* sx * sx ', sep = ' ')
F4 <- paste(A4, '* sqrt(sy^2 + sz^2) * sqrt(sy^2 + sz^2) ', sep = ' ')
F5 <- paste(A5, '* syz * syz ', sep = ' ')
F6 <- paste(A6, '* sx * sqrt(sy^2 + sz^2) ', sep = ' ')             ### what about sz?!?!?!

F.tot <- paste(F1 , F2, F3 , F4 , F5, F6, sep = ' + ') 

print(paste('Workbench Ansys User defined result ', ' ', sep = ' '))
print(F.tot)

# Plotting
x<-y<-z<-seq(-1 * R.t.p * scale ,R.t.p * scale ,length=60)


M <- mesh(x, y, z)
R <- with (M, A1 * x + A2 * y+ A3 * x^2 + A4 * y^2 + A5 * z^(2) - A6 * x * y)
 A2 * 1 + A4 * 1
isosurf3D(x, y, z, colvar = R, level = c(1), col = c("red"), bty = "g", alpha=0.9, ticktype = "detailed", theta = 33)


}

tsai.wu.orthotropic <- function(sigma.allw, m=1.14 ){
scale = 2
library(plot3D)

# R.t tension strenght (3)
R.t <- sigma.allw  ## sigma.allw.matrix is dummy 


# R compression compression strenght (3)
R.c <-  R.t * m

# R shear (3)
R.sh <- 0.577 * R.t          # according stommel2011p175

# R biaxial (3)
R.bi <- R.t          # no influence on biaxial behavior assumed

# Linear Factors
A1 <- (1 / R.t[1]) - (1 / R.c[1])             # A1  
A2 <- (1 / R.t[2]) - (1 / R.c[2])             # A1  
A3 <- (1 / R.t[3]) - (1 / R.c[3])             # A1  
A4 <- A5 <- A6 <- 0

# quadratic factors
A11 <- (1/ (R.t[1] * R.c[1]))   
A22 <- (1/ (R.t[2] * R.c[2])) 
A33 <- (1/ (R.t[3] * R.c[3])) 
A44 <- 1 / (R.sh[1]^2)   
A55 <- 1 / (R.sh[2]^2)  
A66 <- 1 / (R.sh[3]^2)  

# Bilinear factors
A12 <- (1 / (2 * R.bi[1]^2)) * (1 - R.bi[1] * (A1 + A2) - (R.bi[1]^2 * (A11 + A22)) )
A13 <- (1 / (2 * R.bi[2]^2)) * (1 - R.bi[2] * (A1 + A3) - (R.bi[2]^2 * (A11 + A33)) )
A23 <- (1 / (2 * R.bi[3]^2)) * (1 - R.bi[3] * (A2 + A3) - (R.bi[3]^2 * (A22 + A33)) )


# Plotting
x<-y<-z<-seq(-1 * R.t * scale ,R.t * scale ,length=60)


M <- mesh(x, y, z)
R <- with (M, A1 * x + A2 * y+ A3 * x^2 + A4 * y^2 + A5 * z^(2) - A6 * x * y)
 A2 * 1 + A4 * 1
isosurf3D(x, y, z, colvar = R, level = c(1), col = c("red"), bty = "g", alpha=0.9, ticktype = "detailed", theta = 33)


}

## _________________________________________________________________________________________________Isotropic Failure Criteria
conic.crit.fun.stommel <- function (R1=1,R2=1,R3=1,m=1) {
library(plot3D)
## According to stommel  ! faulty
scale <- max(R1,R2,R3) * 1.5


x<-seq(-1 *scale * R1 ,scale * R1 ,length=60)
y<-seq(-1 *scale * R2 ,scale * R2 ,length=60)
z<-seq(-1 *scale * R3 ,scale * R3 ,length=60)
r<-outer(x,y,function(x,y,z=0,m=1) ((m-1)/(2*m) * (x+y+z) + ((m+1)/(2*m*sqrt(2)))* sqrt((x+y+z)^2 + 1/(2*m) * ((x-y)^2 + (z-x)^2 + (y-z)^2)))-1)

par(mfrow = c(1, 2))
contour(x,y,r,levels = c(0), ylim=c(-1*scale,scale),  xlim=c(-1*scale,scale))
abline(h=0)
abline(v=0)
grid()



M <- mesh(x, y, z)
m <- 1.3
R <- with (M, ((m-1)/(2*m) * (x+y+z) + ((m+1)/(2*m*sqrt(2)))* sqrt((x+y+z)^2 + 1/(2*m) * ((x-y)^2 + (z-x)^2 + (y-z)^2)))-1)

#slice3D(x, y, z, colvar = p, xs = 0, ys = c(-4, 0, 4), zs = NULL)

isosurf3D(x, y, z, colvar = R, level = 0, col = "red", bty = "g", alpha=0.9, ticktype = "detailed")

F1 <- paste(A1, '* sx ', sep = ' ')
F2 <- paste(A2, '* sqrt(sy^2 + sz^2) ', sep = ' ')
F3 <- paste(A3, '* sx * sx ', sep = ' ')
F4 <- paste(A4, '* sqrt(sy^2 + sz^2) * sqrt(sy^2 + sz^2) ', sep = ' ')
F5 <- paste(A5, '* syz * syz ', sep = ' ')
F6 <- paste(A6, '* sx * sqrt(sy^2 + sz^2) ', sep = ' ')             ### what about sz?!?!?!

F.tot <- paste(F1 , F2, F3 , F4 , F5, F6, sep = ' + ') 
print(paste('Workbench Ansys User defined result ', ' ', sep = ' '))
print(F.tot)

}

conic.crit.fun.ehrenstein <- function (R1=1,R2=1,R3=1,m=1.3, scale=2) {
library(plot3D)
## Plotting vector of booleans first entry 2D, second entry 3D plot
## According to ehrenstein2011p251
scale <- max(R1,R2,R3) * scale


x<-seq(-1 *scale * R1 ,scale * R1 ,length=60)
y<-seq(-1 *scale * R2 ,scale * R2 ,length=60)
z<-seq(-1 *scale * R3 ,scale * R3 ,length=60)


r<-outer(x,y,function(x,y,z=0) (1/(2*m)) * ((m-1) * (x+y+z) + ((m+1)/(sqrt(2))) * sqrt(((x-y)^2 + (z-x)^2 + (y-z)^2)))  -1)
par(mfrow=c(1,2))

contour(x,y,r,levels = c(0), ylim=c(-1*scale,scale),  xlim=c(-1*scale,scale), main = 'conic ehrenstein')
abline(h=0)
abline(v=0)
grid()

print(paste('User Defined Result Workbench for conic criteion ehrenstein2011p251 for m = ', m))
print(paste((1/(2*m)),  ' * (', ((m-1)), ' * (s1+s2+s3) + ', ((m+1)/(sqrt(2))), ' * sqrt(((s1 - s2)^2 + (s3 - s1)^2 + (s2 - s3)^2)))'         , sep = ''))


M <- mesh(x, y, z)
R <- with (M, (1/(2*m)) * ((m-1) * (x+y+z) + ((m+1)/(sqrt(2))) * sqrt(((x-y)^2 + (z-x)^2 + (y-z)^2)))  -1      )


isosurf3D(x, y, z, colvar = R, level = 0, col = "red", bty = "g", alpha=0.9, ticktype = "detailed")



}

drucker.prager.crit.fun <- function (R1=1,R2=1,R3=1,m=1.3, scale=2) {
library(plot3D)
## Plotting vector of booleans first entry 2D, second entry 3D plot
## According to wikipedia
scale <- max(R1,R2,R3) * scale


x<-seq(-1 *scale * R1 ,scale * R1 ,length=60)
y<-seq(-1 *scale * R2 ,scale * R2 ,length=60)
z<-seq(-1 *scale * R3 ,scale * R3 ,length=60)


par(mfrow = c(1, 2))
r<-outer(x,y,function(x,y,z=0) ((m-1)/(2) * (x+y+z) + ((m+1)/2)* sqrt(1/2 * ((x-y)^2 + (z-x)^2 + (y-z)^2)))-1)
contour(x,y,r,levels = c(0), ylim=c(-1*scale,scale),  xlim=c(-1*scale,scale), main = 'drucker prageer wiki')
abline(h=0)
abline(v=0)
grid()


print('User Defined Result Workbench for Drucker_Prager Wiki')
print('has to be defined :(')

M <- mesh(x, y, z)
R <- with (M, ((m-1)/(2) * (x+y+z) + ((m+1)/2)* sqrt(1/2 * ((x-y)^2 + (z-x)^2 + (y-z)^2)))-1)

#slice3D(x, y, z, colvar = p, xs = 0, ys = c(-4, 0, 4), zs = NULL)
isosurf3D(x, y, z, colvar = R, level = 0, col = "red", bty = "g", alpha=0.9, ticktype = "detailed")

}

parabolic.crit.fun.ehrenstein <- function (R1=1,R2=1,R3=1,m=1.3, scale=1, plots=FALSE) {
library(plot3D)
## Plotting vector of booleans first entry 2D, second entry 3D plot
## According to ehrenstein2011p251
scale <- max(R1,R2,R3) * scale


x<-seq(-1 *scale * R1 ,scale * R1 ,length=60)
y<-seq(-1 *scale * R2 ,scale * R2 ,length=60)
z<-seq(-1 *scale * R3 ,scale * R3 ,length=60)


r<-outer(x,y,function(x,y,z=0) (((m-1)/(2*m)) * (x+y+z) + sqrt( (((m-1)^2)/(4*m^2)) * (x+y+z)^2 + (1/(2*m)) * ((x-y)^2 + (z-x)^2 + (y-z)^2)))  -1   )  

par(mfrow=c(1,2))
contour(x,y,r,levels = c(0), ylim=c(-1*scale,scale),  xlim=c(-1*scale,scale), main = 'parabolic ehrenstein')
abline(h=0)
abline(v=0)
grid()


#print(paste('User Defined Result Workbench for parabolic criteion ehrenstein2011p251 for m = ', m))
math.string <- paste(((m-1)/(2*m)), ' * (s1+s2+s3) + sqrt( ', (((m-1)^2)/(4*m^2)),  ' * (s1 + s2 + s3)^2 + ', (1/(2*m)), ' * ((s1 - s2)^2 + (s3 - s1)^2 + (s2 - s3)^2))'            , sep = '')


M <- mesh(x, y, z)
R <- with (M, (((m-1)/(2*m)) * (x+y+z) + sqrt( (((m-1)^2)/(4*m^2)) * (x+y+z)^2 + (1/(2*m)) * ((x-y)^2 + (z-x)^2 + (y-z)^2)))  -1      )

if (plots == TRUE) {
#isosurf3D(x, y, z, colvar = R, level = 0, col = "red", bty = "g", alpha=0.9, ticktype = "detailed")
}
return(list("math.string"=math.string))
}

critical.strain.menges <- function (polymertype, fiber.w.per, temperatures, T.g){
## Critical Strain by Menges, polymertypes semicr = 1, amorphous = 2, PS =3
  if (polymertype == 1) {
    crit.strain <- c()
    for (i in temperatures){
      if (i >= T.g){crit.strain <- append(crit.strain, 0.02)}
      else {crit.strain <- append(crit.strain, 0.008)}
    }
    of <- data.frame(temp=temperatures, crit.strain= crit.strain * (1 - 1.5 * (fiber.w.per/100)) )  
  }
  if (polymertype == 2) {of <- data.frame(temp=temperatures, crit.strain= 0.008 * (1 - 1.5 * (fiber.w.per/100)) )}
  if (polymertype == 3) {of <- data.frame(temp=temperatures, crit.strain= 0.002 * (1 - 1.5 * (fiber.w.per/100)) )}

  of
}

fatique.strength.temp.bayer <- function (strength, K.fac.dyn,temperatures=NULL, cycles.num=c(1,1e1,1e2,1e3,1e4,1e5,1e6,1e7)){
## Function to calculate strength over cycles for each given temperature
  func_stress_amplitude <- function(x,strength) strength * (1-(1-K.fac.dyn)*(log10(x)/7))
  
  if (is.null(temperatures)){temperatures <- c(0)}
  df <- data.frame(temp=temperatures, strength=strength)
  of <- data.frame()
  
  for (i in df$temp){   
      of <- rbind(of,
	    data.frame(temp=rep(i, length(cycles.num))
		      ,cycles=cycles.num
		      ,strength=func_stress_amplitude(cycles.num,df[df$temp == i, "strength"])))
  }
  if (is.null(temperatures)){of[,!colnames(of)=="temp"]}
  else {of}
}

sig.ak <- function(temp,sig.wk,Mean.stress.limit,stress.abbr.str=list(siga="siga",sigu="sigu",sigm="sigm",temp='BFE')){
    # function to calculate sig.ak "Part Limiting Stress Amplitude with Mean stress influence" according korte2018p231
    # sig.ak = sig.wk * K.ak
    # sig.wk = The long term durability limit for 10^7 (derived from short term durablity F * Aw)
    # K.ak   = Mean stress reduction factor, dependent on the type of overload situation F1-4
        df.wk.temp <- data.frame(temp=temp,sig.wk=sig.wk)
        sig.wk.fit <- fit.function.limits(df.wk.temp,degree=3,xname=stress.abbr.str$temp,str.nmbr.digits=5,x.col=1,y.col=2)
        sig.wk.over.temp.math.string <- sig.wk.fit$math.string
        sig.wk.temp.fun <- sig.wk.fit$max.limit.fun
        M.sig <- sig.wk / Mean.stress.limit

        K.ak.fun <-      list(  F1=function(sig.m,temp) 1 - M.sig[1] * sig.m / sig.wk.temp.fun(temp),                 #F1 sig.m const
                                F2=function(sig.m,sig.a) 1 / (1 + (M.sig[1] * sig.m / sig.a)),        #F2 R const
                                F3=function(sig.u,temp) (1 - (M.sig[1] * sig.u/sig.wk.temp.fun(temp)))/(1+M.sig[1]),       #F3 sig.u const
                                F4=function(sig.o,temp) (1 - (M.sig[1] * sig.o/sig.wk.temp.fun(temp)))/(1+M.sig[1])        #F4 sig.o const
                                )
        K.ak.string <- list (   F1=paste('1-(',M.sig[1],'*',stress.abbr.str$sigm,'/(',sig.wk.over.temp.math.string,'))'),
                                F2=paste('1 / (1 +(', M.sig[1],' *', stress.abbr.str$sigm ,'/', stress.abbr.str$siga,'))'),
                                F3=paste( '(1 - (',M.sig[1],' * ',stress.abbr.str$sigu,'/(',sig.wk.over.temp.math.string,')))/(1+',M.sig[1],')'),
                                F4=paste( '(1 - (',M.sig[1],' * ',stress.abbr.str$sigo,'/(',sig.wk.over.temp.math.string,')))/(1+',M.sig[1],')' )
                                )     
        string <- list ( F1= paste('(',K.ak.string$F1,')*(',sig.wk.over.temp.math.string,')',sep='') ,
                                F2= paste('(',K.ak.string$F2,')*(',sig.wk.over.temp.math.string,')',sep='') ,
                                F3= paste('(',K.ak.string$F3,')*(',sig.wk.over.temp.math.string,')',sep='') ,
                                F4= paste('(',K.ak.string$F4,')*(',sig.wk.over.temp.math.string,')',sep='') 
                                )
                                
        fun <- list (    F1=function(sig.m,temp) K.ak.fun$F1(sig.m,temp)* sig.wk.temp.fun(temp),
                                F2=function(sig.m,sig.a,temp) K.ak.fun$F2(sig.m,sig.a)* sig.wk.temp.fun(temp),
                                F3=function(sig.u,temp) K.ak.fun$F3(sig.u,sig.wk)* sig.wk.temp.fun(temp),
                                F4=function(sig.o,temp) K.ak.fun$F4(sig.o,sig.wk)* sig.wk.temp.fun(temp)
                                )
        return(list(K.ak.fun=K.ak.fun,K.ak.string=K.ak.string,string=string,fun=fun,sig.wk.fit=sig.wk.fit,M.sig=M.sig))                       
}


p.swt.to.N <- function(E,e.f,S.u,pswt) {
	# function to  convert a calculated pswt value to expected cycle number for metals
	# E= Modulus
	# e.f = strain at break
	# S.u = stress at break
	# example p.swt.to.N(E=82600,e.f=.305,S.u=240,pswt=551)
	
    sig.f <- (0.623*(S.u/E)^0.832)*E
    ep.f <- (0.0196*(log(1+e.f))^0.155) * (S.u/E)^(-0.53)  #true failure strain

    b <- -0.09
    c <- -0.56

    n <- b/c
    K <- sig.f/(ep.f^n)

    coeff <- list()

    coeff$strength.coeff <- sig.f
    coeff$strength.expo <- b
    coeff$ductility.coeff <- ep.f
    coeff$ductility.expo <- c
    coeff$cyc.strength.coeff <- K
    coeff$cyc.strain.hard.expo <- n

    P.swt.fun <- function(coeff,E,N,sig.f,ep.f){
        (sig.f^2 *(2*N)^(2*coeff$strength.expo) + sig.f * ep.f *E * (2*N)^(coeff$strength.expo+coeff$ductility.expo) )^0.5
    }

    pswt.to.n <- approxfun(P.swt.fun(coeff=coeff,E=E,N=10^(seq(-1,7,length.out=500)),sig.f=sig.f,ep.f=ep.f),10^(seq(-1,7,length.out=500)))
    return(pswt.to.n(pswt))  # export N
}


## _________________________________________________________________________________________________________________ Polymer stiffness shifting WLF



wlf.equation <- function(T.0 = 23,T.g=50, T=seq(-50,40,1), clte=0,base.time=1300,base.temp=110,plot=FALSE) {

    ## equation for wlf equation based on T0 = Tg
    ## base time in h or in general unitless
    
    # T.0 should be set to Tg of the plastic material    # T is the temperature region of interest    # beta is the coefficient of linear expansion CLTE# result is c(log10(a.t), k, vertical Factor)
    T.0 <- T.g+43 + 273.18 
    
    T.g <- 50+ 273.18
    T <- T + 273.18
    delta.T <- (T - T.0)
    C.1 <- 8.86
    C.2 <- 101.6
    log.a.t <- ((-1 * C.1 * delta.T)/(C.2 + delta.T))

    C.1.g <- 15
    C.2.g <- 50

    C.1 <- function (T.0,T.g) {(C.1.g * C.2.g)/(C.2.g + (T.0 -T.g ))}
    C.2 <- function (T.0,T.g) {C.2.g + (T.0 -T.g )}
    
    log.a.t.temp.modif <- ((-1 * C.1(T.0,T.g) * delta.T)/(C.2(T.0,T.g) + delta.T))

    
    k <- (log.a.t)/((1/T) - (1/T.0))  # calc of arrhenius value k, just for information
    df <- data.frame(T.kelvin=T,T=T-273,log.a.t=log.a.t,log.a.t.temp.modif=log.a.t.temp.modif,k=k)  #assembly of data.frame
    
    ## Calculation of scaled time.
    df$time.scaled <- 10^log.a.t *  (base.time /(10^approxfun(df$T,df$log.a.t)(base.temp)))
    df$time.scaled.modif <- 10^log.a.t.temp.modif * (base.time /(10^approxfun(df$T,df$time.scaled.modif)(base.temp)))

    
    
    if (plot == TRUE){
        require(ggplot2)
        # log.at plot
        gg.at <- ggplot(df) + theme_bw() + 
        geom_line(aes(x=T,y=log.a.t,color="const")) +
        geom_line(aes(x=T,y=log.a.t.temp.modif,color="dep")) +
        labs(y="log(at)",x="Temperature [degC]")+
            scale_x_continuous(name="Temperature [degC]",breaks=round(seq(min(df$T),max(df$T),5)),0)  + geom_vline(xintercept=base.temp ) + geom_hline(yintercept=0 )
        print(gg.at)
        
        #time plot
        plot.df <- df[df$T >= base.temp,]
        gg.time <- ggplot(plot.df) + theme_bw() + 
        geom_line(aes(x=T,y=time.scaled,color="C const.")) +
        geom_line(aes(x=T,y=time.scaled.modif,color="C(T)")) +
        labs(y="time scaled",x="Temperature [degC]") +
        scale_y_continuous(name="time scaled [h]", breaks=seq(0,base.time,50)) +
        scale_x_continuous(name="Temperature [degC]",breaks=round(seq(min(plot.df$T),max(plot.df$T),1)),0)  + geom_vline(xintercept=base.temp )
        ggsave(paste("ttsp-T0_",T.0-273.18,"--basetime-",base.time,".png",sep=""),plot=gg.time, width =40, height = 30, units="mm",scale=4)
        print(gg.time)
    }

   return(df)
}
#wlf.equation(T.g=50,T=seq(50,130,1),plot=TRUE,base.time=1375,base.temp=120)
#wlf.equation(T.g=50,T=seq(50,130,1),plot=TRUE,base.time=3280,base.temp=110)


## _____________________________________________________________________________________________________________________________  Material Utils

lin.exp.to.vol.exp <- function (lin.sh){
## Conversion of a linear expansion like CLTE[1/K] or linear shrinkage value
## to volumetric expansion/shrinkage
1-(1-(lin.sh))^3
}

vol.exp.to.lin.exp <- function (vol.sh){
## Conversion of a volumetric expansion like CLTE[1/K] or volumetric shrinkage value
## to linear expansion/shrinkage
1-(1-(vol.sh))^(1/3)
}

strain.stress.convert <- function (strain, stress, mode, temperatures=NULL) {
## Summary function for calculation between stress and strain types
## Modes eng->true
  of <- data.frame(strain,stress)
  if (mode == 'eng->true'){
    of[1] <- log(1 + strain)
    of[2] <- stress*(1 + strain)
  }
    
  if (is.null(temperatures) == FALSE) {
    of <- data.frame(temp=temperatures, strain=of[[1]], stress=of[[2]])
  }

  of
}

lin.el.extract <- function (total.strain, stress, temperatures=NULL){
## function to remove the elastic part of strains from stress-strain data
## ALSO gives back the corresponding lin. elastic E-Modulus for each temperature
  if (is.null(temperatures) == TRUE){temp = rep(0, length(total.strain))}			# if no additonal temperature column is given, create dummy temp vector with zeros
  else {temp <- temperatures}									# if temperatures are provided put into temp vector
  E.el <- rep(0, length(total.strain))

  of <- data.frame()										# initiate ouput dataframe 'of'
  df <- data.frame(temp=temp, tot.strain=total.strain, stress=stress, e.el=E.el)			# create working dataframe 'df' where operations are applied on

  for (i in unique(df$temp)){									# for each temperature in df
    df.temp <- subset(df, df$temp == i,)							# extract subset of strains/stresses for each temperature
    df.temp <- df.temp[order(df.temp$tot.strain),]						# order it by strains, just to get sure to calc right E_el in next step
    
    if (!(df.temp$tot.strain[1] <= 0)) {E <- df.temp$stress[1] / df.temp$tot.strain[1]}		# if first supplied strain is NOT zero or less create E_el from first point
    else {E <- df.temp$stress[2] / df.temp$tot.strain[2]}					# else take the second data point to create E_el
    df.temp$el.strain <- (df.temp$stress) / E							# e_pl = e_tot - e_el, where e_el is derived from trueStress/trueE
    df.temp$pl.strain <- df.temp$tot.strain - (df.temp$stress) / E							# e_pl = e_tot - e_el, where e_el is derived from trueStress/trueE
    df.temp$e.el <- E
    of <- rbind(of,df.temp)  									# rbind df.temp to output file
  }
 
  of <- of[!(of$pl.strain < 0),]								# output of file without negative plastic strain rows

  if (is.null(temperatures) == TRUE){of[,!colnames(of)=="temp"]}
  else {of}
}

secant.modulus <- function (strain, stress, temperatures=NULL){
## Function to extract secant modulus from stress/strain data
## returns data.frame((temp),stress,strain,e.sec)
  if (is.null(temperatures) == TRUE){temp = rep(0, length(strain))}				# if no additonal temperature column is given, create dummy temp vector with zeros
  else {temp <- temperatures}									# if temperatures are provided put into temp vector
  of <- data.frame()										# initiate ouput dataframe 'of'
  df <- data.frame(temp=temp, strain=strain, stress=stress)					# create working dataframe 'df' where operations are applied on
  df <- df[!(df$strain == 0),]									# Remove rows with zero strain to avoid inf (E.s = sigma/strain.short.temp.ref
  
  df$e.sec <- df$stress / df$strain
  
  if (is.null(temperatures) == TRUE){ df[,!(colnames(df)==c("temp","","",""))] }
  else { df[,!(colnames(df)==c("","","",""))] }
}


stress.strain.function.D1.D2.Menges <- function (E,yield.strain,strain.at.break,yield.stress,stress.at.break,create.plot=FALSE) {
    require('minpack.lm')
    ## Function to fit a stress-strain function from E, yield and breakpoint for crystalline Thermoplastics according menges2011p193
    #testers:
    #s.s.fun <- stress.strain.function.D1.D2.Menges(20000,0.01,0.018,90,97)
    #s.s.fun(0.01)
    s.s.points <- data.frame(strain=c(0,yield.strain,strain.at.break),stress=c(0,yield.stress,stress.at.break))


    para <- list(D1=1,D2=1)

    resiFun <- function(strain, para,stress,E) (stress - (E *strain *(1-para$D1*strain)/(1+para$D2 * strain)))


    s.s.fit <- nls.lm(par=para, fn = resiFun, strain = s.s.points$strain,stress = s.s.points$stress, E=E
                ,  control = nls.lm.control(nprint=1))

                
    D1 <- s.s.fit$par$D1
    D2 <- s.s.fit$par$D2

    stress.fun <- function(strain) (E * strain * (1-D1*strain)/(1+D2*strain))
                        
    if (create.plot==TRUE) {
        plot.df <- data.frame(strain=seq(0,strain.at.break,length.out=100),stress=stress.fun(seq(0,strain.at.break,length.out=100),E,D1,D2))
        ggplot(data=plot.df) + theme_bw() + geom_line(aes(x=strain,y=stress))
    }
    return(stress.fun)
}    


stress.strain.function.Schoeche <- function (E,yield.strain,strain.at.break,yield.stress,stress.at.break,create.plot=FALSE) {
    require('minpack.lm')
    ## Function to fit a stress-strain function from E, yield and breakpoint for crystalline Thermoplastics according menges2011p193
    #testers:
    #s.s.fun <- stress.strain.function.D1.D2.Menges(20000,0.01,0.018,90,97)
    #s.s.fun(0.01)
    s.s.points <- data.frame(strain=c(0,yield.strain,strain.at.break),stress=c(0,yield.stress,stress.at.break))


    para <- list(ep=1)

    resiFun <- function(strain, para,stress,E) (stress - (E *para$ep *(1-exp(-1*strain/para$ep))))


    s.s.fit <- nls.lm(par=para, fn = resiFun, strain = s.s.points$strain,stress = s.s.points$stress, E=E
                ,  control = nls.lm.control(nprint=1))

                
    ep <- s.s.fit$par$ep

    stress.fun <- function(strain) (E * ep * (1-exp(-1*strain/ep)))
                        
    if (create.plot==TRUE) {
        require(ggplot2)
        plot.df <- data.frame(strain=seq(0,strain.at.break,length.out=100),stress=stress.fun(seq(0,strain.at.break,length.out=100)))
        print(
            ggplot() + theme_bw() + geom_line(data=plot.df,aes(x=strain,y=stress))+ geom_point(data=s.s.points, aes(x=strain,y=stress))
        )
    }
    return(stress.fun)
}    

s.s.from.G.temp <- function (shearmod.file,stress.at.break.ref,strain.at.break.ref,ref.temp) {
   ## function to create generic stress-strain curves for thermoplastics based on shear modulus over temperature curves and one break point 

    g.temp <- read.table(shearmod.file, sep = ',', header = TRUE) #("type"= storage,loss,tand | "temp" [degC] | value = [MPa]/[ratio])

    colnames(g.temp) <- c("type","temp","value")



    df.mod.storage <- g.temp[g.temp$type == "storage",][,c(2,3)]
    colnames(df.mod.storage) <- c("temp","g.el")
    
    
    if (is.na(match(ref.temp,df.mod.storage$temp))) {  # if reftemp is not part of input, then interpolate
    df.mod.storage <- rbind(df.mod.storage,c(ref.temp,approxfun(df.mod.storage$temp,df.mod.storage$g.el)(23)))
    df.mod.storage <- df.mod.storage[order(df.mod.storage$temp),] # sort after additon
    }


    df.mod.storage$poisson <- 0.3 + 0.2 * (1- (df.mod.storage$g.el/max(df.mod.storage$g.el)))
    df.mod.storage$e.el <- 2 * df.mod.storage$g.el * (1+df.mod.storage$poisson)


    E.ref <- approxfun(df.mod.storage$temp,df.mod.storage$e.el)(ref.temp)
    

    df.s.s <- data.frame(temp=c(),strain=c(),stress=c())

    yield.stress <- stress.at.break.ref/2
    yield.strain <- yield.stress/E.ref
    s.s.fun <- stress.strain.function.Schoeche(E.ref,yield.strain,strain.at.break.ref,yield.stress,stress.at.break.ref,FALSE)

    for (i in 1:length(df.mod.storage$e.el)) {   # simple modulus scaling for other temps

    strain.temp <- seq(0,strain.at.break.ref*E.ref/df.mod.storage$e.el[i],length.out=20)
    stress.temp <- s.s.fun(strain.temp)*df.mod.storage$e.el[i]/E.ref

    df.s.s <- rbind(df.s.s, data.frame(temp=df.mod.storage$temp[i],strain=strain.temp,stress=stress.temp))

    }

    return(df.s.s)
}

#s.s.from.G.temp(shearmod.file,stress.at.break.ref,strain.at.break.ref,ref.temp)


fit.e.el.from.tan.mod <- function (s.s.short,strain.in.percent=TRUE) {
    ## Function to extract elastic modulus for temp. depended s.s curves by fit of tangential modulus at 0 strain
    ## This is done to avoid wrong elastic-modulus results from scattering multip-point test data
    
    df.e.tan.over.strain <- data.frame()
    e.el <- data.frame()
    strain.fac <- 1
    if (strain.in.percent) {strain.fac <- 100}
    for (temp in unique(s.s.short$temp)){
        
        df.e.tan.over.strain.temp <- data.frame(temp=temp,strain=head(s.s.short[s.s.short$temp == temp,]$strain,-1),e.tan=diff(s.s.short[s.s.short$temp == temp,]$stress)/(diff(s.s.short[s.s.short$temp == temp,]$strain /strain.fac))) #calculate e.tan points and store to df
        df.e.tan.over.strain.temp <- df.e.tan.over.strain.temp[df.e.tan.over.strain.temp$e.tan > 0,]  # remove negative e.tan values
        df.e.tan.over.strain.temp$fit.results <- fit.function.limits(df.e.tan.over.strain.temp[,c(2,3)],degree=3,xname='strain',str.nmbr.digits=20)$fit.results
        df.e.tan.over.strain <- rbind(df.e.tan.over.strain,df.e.tan.over.strain.temp)
        e.el.temp <- data.frame(temp=temp,e.el=df.e.tan.over.strain.temp$fit.results[1])
        e.el <- rbind(e.el,e.el.temp)

    }

    color.scale.temp <- rev(hcl(h = seq(1,230,length.out=length(unique(s.s.short$temp))), c = 100, l = 75))  # create temperature color scale for ggplots

    ##create e.tan plot with fit results from raw results
    gg.e.tan <- (ggplot(df.e.tan.over.strain) + theme_bw() + scale_colour_manual(name='temp [degC]',values = color.scale.temp) +
        geom_point(aes(x=strain,y=e.tan,group=temp,color=as.factor(temp))) +
        geom_line(aes(x=strain,y=fit.results,group=temp,color=as.factor(temp)))+
        scale_x_continuous(name="eng. strain [mm/mm]")+
        scale_y_continuous(name="Tangent Modulus [MPa]"))
        
    return(list(e.el=e.el,df.e.tan.over.strain=df.e.tan.over.strain,gg.e.tan=gg.e.tan))
}


##________________________________________________________________________________________ Electromagnetics

rel.cs.from.three.pts <- function(pt1, pt2, pt3) {
## Function to create relative local coordinate system in Maxwell by three
## Point Measurement in Rhino
## Do not forget to set the Reference CS in Maxwell to Global :)
  origin <- pt1
  x.axis <- pt2 - pt1
  y.point <- pt3 - pt1
  
  print('Origin:')
  print(paste(pt1[1],",",pt1[2],",",pt1[3]))
  print('X Axis:')
  print(paste(x.axis[1],",",x.axis[2],",",x.axis[3]))
  print('Y Point:')
  print(paste(y.point[1],",",y.point[2],",",y.point[3]))

}

pt1 <- c( -5.955,34.225,0.126)
pt2 <- c(-5.955,34.325,0.126 )
pt3 <- c(  -5.855,34.225,0.126)

rel.cs.from.three.pts(pt1, pt2, pt3)



## Located in criteria.lib.R::axial.step.fun.mag.rot


## Located in criteria.lib.R::axial.step.fun.mag.rot

axial.step.fun.mag.rot <- function(x.range,y.range,Time,k,plot3Drot=TRUE,plotResolution3D=200){
## Function for plotting an rotating half-axial (half disc) oriented magned

library('plot3D')
x <- seq(min(x.range),max(x.range),length.out=plotResolution3D)
y <- seq(min(y.range),max(y.range),length.out=plotResolution3D)
M <- mesh(x, y)

##Merged for Maxwell
z.comp <- function(x,y,Time,k){ ((2/(1+exp(-2*k*(x * cos(Time*0.01745329) + y * sin(Time*0.01745329))))))-1}
x.comp <- 1-abs(z.comp(x=x,y=0,Time=0,k=k) )

if (plot3Drot) {surf3D(M$x, M$y, z.comp(M$x,M$y,Time,k),bty="b2")}
#plot x-direction z.comp
plot(x,z.comp(x=x,y=0,Time=0,k=k),ylab='z.comp', xlab='x', type='l'); grid()
lines(x, x.comp, col='red')

}
#example
#fac <- 1/1000
#x.range <- y.range <- c(-4*fac,4*fac); Time <- seq(0,120,length.out=4); magnitude <- 242048.14 /2
#par(mfrow=c(1,1)) ; axial.step.fun.mag.rot(x.range,y.range,Time=0,k=500,plot3Drot=FALSE,plotResolution3D=2000)
#par(mfrow=c(2,2)); for (i in Time){axial.step.fun.mag.rot(x.range,y.range,i,2000,plot3Drot=TRUE,plotResolution3D=200)}






## _______________________________________________________________________________________  Mathematical Helpers
## Average of values by Integral  average.by.integral
average.by.integral <- function (table){
  tempfun <- approxfun(table[[1]],table[[2]])
  surf <- integrate(tempfun,min(table[[1]]),max(table[[1]]),subdivisions=1000)
  avg <- surf[["value"]] / (max(table[[1]]) - min(table[[1]]))
  plot(seq(min(table[[1]]),max(table[[1]]),length.out=200), tempfun(seq(min(table[[1]]),max(table[[1]]),length.out=200)), type='l',xlab= 'x', ylab ='y')
  points(table[[1]], table[[2]])
  polygon(seq(min(table[[1]]),max(table[[1]]),length.out=200), tempfun(seq(min(table[[1]]),max(table[[1]]),length.out=200)), border=NA, col="grey")
  abline(h=avg)
  grid()
  print(paste('Average value between x ', min(table[[1]]), ' and ', max(table[[1]]), ' is: ', avg ))
}


df.scalar.res.from.df.csv.list <- function(df.csv.list, exclude.df=c(), exclude.col.num=c()){
## Function to summarize some scalar results for each column in all df of df in a list
df.scalar.res <- maximum <- minimum <- p2p <- average <- df.name <- c()

for (csv.name in names(df.csv.list)){
  if (csv.name %in% exclude.df) next
  df.name <- append(df.name,rep(csv.name,length(names(df.csv.list[[csv.name]])[-1])))
  col.num <- 2
  for (col.name in names(df.csv.list[[csv.name]])[-1]){
    if (col.num %in% exclude.col.num) next
    cur.set <- df.csv.list[[csv.name]][[col.name]]
    df.scalar.res <- rbind(df.scalar.res, col.name)
    p2p <- rbind(p2p, max(cur.set) - min(cur.set))
    maximum <- rbind(maximum, max(cur.set))
    minimum <- rbind(minimum, min(cur.set))
    average <- rbind(average, mean(cur.set))
    col.num <- col.num + 1

}}
#df.scalar.res <- (data.frame(df=df.name,col.name=df.scalar.res[,1],max=maximum,min=minimum,mean=average,p2p=p2p))
return((data.frame(df=df.name,col.name=df.scalar.res[,1],max=maximum,min=minimum,mean=average,p2p=p2p)))
}


## p2p 
p2p <- function(x) {max(x)-min(x)}

## Integral from set of curves (curveID (mostly temperatures), x, y)
integral.multicurve <- function (df) {
  integral.values <- rep(0, length(unique(df[[1]])))
  i = 1
  for (curveID in unique(df[[1]])){
      integral.values[i] <- integrate(approxfun(subset(df, df[1] == curveID)[[2]], subset(df, df[1] == curveID)[[3]]),0,max(subset(df, df[1] == curveID)[[2]]))$value
      i = i + 1
  }
  integral.fun <- approxfun(unique(df[[1]]), integral.values)
  integral.table <- data.frame(temp=unique(df[[1]]), integral=integral.values)
  return(list(fun=integral.fun, table=integral.table))
}

## Convert deg Angle into rad angle
deg.to.rad <- function(deg){ (2*pi/360)*deg  }
rad.to.deg <- function(rad) {rad * 180 / pi}


## (Sensing) Estimate linearity error by multi point calibration
multi.pt.calibr <- function(x,y,x.grid=c(-1,0,1),resolution=10000){
      ## function to estimate error by linear sement-wise multi point calibration of ICs based on measured point data
      ## Point data x,y (eg. from simulation) is first fitted to polynominal function, this have to be veryfied by residuals
      ## x.grid is a vector of x values which are root positions for calibration
      x.plot <- seq(min(x.grid),max(x.grid),length.out=resolution)											#create high resolution x vector for output

      fit <- lm(y  ~  x + I(x ^2)+ I(x ^3)+ I(x ^4)+ I(x ^5)+ I(x ^6)+ I(x ^7))										#fit measured values to 7th order polynomial
      fit.fun.poly <- function(x) fit$coefficients[[1]] + x * fit$coefficients[[2]] + x^2 * fit$coefficients[[3]] + x^3 * fit$coefficients[[4]] +	#create function from fit
		x^4 * fit$coefficients[[5]] + x^5 * fit$coefficients[[6]] + x^6 * fit$coefficients[[7]] + x^7 * fit$coefficients[[8]]
      poly.sens.fun <- function(x) (fit$coefficients[[2]]) + x * 2 * fit$coefficients[[3]] + x^2* 3 * fit$coefficients[[4]] + x^3* 4 * fit$coefficients[[5]] + # create sensitivity function (derivative) for poly fit
				x^4* 5 * fit$coefficients[[6]] + x^5* 6 * fit$coefficients[[7]]  + x^6* 7 * fit$coefficients[[8]]  
				
      fit.fun.lin <- approxfun(x.grid, fit.fun.poly(x.grid))												#create multi-linear function by x grid roots
      lin.sens <- diff(fit.fun.lin(x.plot))/diff(x.plot)												#create sensitivity point for linear function
      
      #create output data frames 
      df.signal <- data.frame(x=x.plot, lin=fit.fun.lin(x.plot), poly=fit.fun.poly(x.plot), err=fit.fun.poly(x.plot)-fit.fun.lin(x.plot))
      df.sens <- data.frame(x=rowMeans(embed(x.plot, 2)), lin=diff(fit.fun.lin(x.plot))/diff(x.plot), poly=poly.sens.fun(rowMeans(embed(x.plot, 2))), err=poly.sens.fun(rowMeans(embed(x.plot, 2)))-(diff(fit.fun.lin(x.plot))/diff(x.plot)))

      return(list(signal=df.signal, fit.resi=data.frame(x=x, residual=fit$residuals ), sensitivity=df.sens, fun.poly=fit.fun.poly, fun.lin=fit.fun.lin))
}


## Convert a linear regression result with R^2 to a string to show in graphs
lm_eqn <- function(x,y){
    df <- data.frame(x=x,y=y)
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}


## ______________________________________________________________________________________ Rheology

pvt.tait <- function(b, bs, bm, T.degC, p.mpa){
## According kenedy2013p31
## second value for bs   'bs[2]' can be scaled by a factor to modify the spec. density at RT
## works for multipble temperatures but only one pressure at once
# zytel 70 G35 HSLRA4 BK267
#b <- c(493, 6.67e-8, 4.78e-5, 2.5e-2, 4.05e-9)
#bs <- c(7.28e-4, 9.36e-8, 2.72e8, 1.91e-3)
#bm <- c(7.76e-4, 1.56e-7, 1.79e8, 1.71e-3)

# Delrin 100 NC010
#b <- c(418, 7.67e-8, 9.09e-5, 4e-2, 5.15e-9)
#bs <- c(7.37e-4, 2.64e-7, 3.99e8, 1.1e-7)
#bm <- c(8.28e-4, 4.85e-7, 1.82e8, 4.4e-3)

  b <- append(rep(0,4), b)   # Create b vector
  C <- 0.0894					#Constant
  T <- T.degC + 273.15				# Unit to K
  p <- p.mpa * 1000000				# Unit to Pa

  T.trans <- b[5] + b[6] * p			
  
  df <- data.frame(T=T, V.0=0, B=0, V.t=0)
  
  df[df$T <= T.trans,"V.0"] <- bs[1] + bs[2] * (df[df$T <= T.trans,"T"] - b[5])			#for all Temps smaller than T.Trans 
  df[df$T >  T.trans,"V.0"] <- bm[1] + bm[2] * (df[df$T >  T.trans,"T"] - b[5])			#for all Temps bigger than T.Trans 
  
  
  df[df$T <= T.trans,"B"] <- bs[3] *  exp(-1 * bs[4] * (df[df$T <= T.trans,"T"] - b[5]))
  df[df$T >  T.trans,"B"] <- bm[3] *  exp(-1 * bm[4] * (df[df$T >  T.trans,"T"] - b[5]))
  
  df[df$T <= T.trans,"V.t"] <- b[7] * exp(b[8] * (df[df$T <= T.trans,"T"] - b[5]) - b[9] * p)
  
  df$V.0 *  (1 - C *  log(1 + (p)/(df$B))) + df$V.t  # output specific weight m^3/kg
}
# Example for plotting one Material with pvt.tait
#T <- 23:250
#plot(T, pvt.tait(b,bs,bm,T,0),type='l',xlab='T [degC]', ylab='nu [m^3/Kg]', ylim=c(min(pvt.tait(b,bs,bm,T,0))*0.9,max(pvt.tait(b,bs,bm,T,0))))
#lines(T, pvt.tait(b,bs,bm,T,50),type='l'); lines(T, pvt.tait(b,bs,bm,T,100),type='l');lines(T, pvt.tait(b,bs,bm,T,150),type='l');grid()

alpha.from.pvt <- function(b, bs, bm, T.range, p.mpa, alpha.ref=NULL){
## Helper function to extract approximation of CLTE (alpha) from a tait PVT diagam
## pvt.tait function must be sourced
## argument alpha.ref can be used to set fixed alpha at RT

  nu.ref <- pvt.tait(b,bs,bm,23,p.mpa)
  
  #delta.nu <- nu[2:length(x)] - nu[1:length(x)-1]    #just test -.-
  of <- data.frame(temp=min(T.range):max(T.range))
  #of <- data.frame(temp=of[!(of$temp == 23),])
  of$nu <- pvt.tait(b,bs,bm,of$temp,p.mpa)
  
  of.nu.tan <- c()
  for (i in 2:length(of$nu)){
    of.nu.tan <- append(of.nu.tan,((of$nu[i]/of$nu[i-1]) - 1 ))
  }
  of.nu.tan <- append(c(of.nu.tan[1]),of.nu.tan)
  of$beta.tan <- of.nu.tan
  of$alpha.tan <- vol.exp.to.lin.exp(of$beta.tan)
  
  of.nu.sec <- c()
  of$beta.sec <- ((of$nu/nu.ref - 1) * 1 / (min(T.range):max(T.range) - 23))
  of$alpha.sec <- vol.exp.to.lin.exp(of$beta.sec)
  
  if (!is.null(alpha.ref)) {
    fac <- alpha.ref / of[of$temp == 23,'alpha.tan']
    of$alpha.sec <- of$alpha.sec * fac
    of$alpha.tan <- of$alpha.tan * fac
  }
  par(mfrow=c(1,1))
  plot(of$temp,of$alpha.tan*1e6,xlab='temperature [degC]', ylab='alpha [1/10^-6]',type='l')
  lines(of$temp,of$alpha.sec*1e6,col='red')
  grid(); legend('topleft', ncol = 2, cex = 0.6, legend = c('tan.','sec'), col = c('black','red'), lty = 1, )
  
  of[,!(colnames(of)==c('','nu','beta.tan','','beta.sec',''))]

}

# Example for alpha.from.pvt
# Delrin 100 NC010
#b <- c(418, 7.67e-8, 9.09e-5, 4e-2, 5.15e-9)
#bs <- c(7.37e-4, 2.64e-7, 3.99e8, 1.1e-7)
#bm <- c(8.28e-4, 4.85e-7, 1.82e8, 4.4e-3)
#T.range <- c(-40,60)
#p.mpa <- 0

#alpha.from.pvt(b, bs, bm, T.range, p.mpa, alpha.ref=110e-6)
#points(c(-40,23,55),c(60,110,120))

## _______________________________________________________________________________________________ ggplot helpers

extract.extreme.from.molten <- function(mdf,min.and.max=FALSE) {
## Function to collapse a molten dataframe only to its most extreme (pos/neg max) values per variable, especially for barplots with ggplot
  colnames(mdf) <- c('id', 'variable', 'value')
  extr.err.df <- data.frame()
  for (id in unique(mdf$id)){
      mdf.i <- mdf[mdf$id==id,]
      i <- 0
      for (variable in unique(mdf.i$variable)){
	i <- i + 1
	mdf.i.i <- mdf.i[mdf.i$variable==variable,]
	
	if (!min.and.max){
	    if (abs(max(mdf.i.i$value)) >= abs(min(mdf.i.i$value))){    # if absolute val of smallest is lower than highest extract positive extreme
		extr.err.df <- rbind(extr.err.df, mdf.i.i[mdf.i.i$value==max(mdf.i.i$value),][1,])
	    }
	    if (abs(max(mdf.i.i$value)) < abs(min(mdf.i.i$value))){    # if absolute val of smallest is higher than highest extract negative extreme
		extr.err.df <- rbind(extr.err.df, mdf.i.i[mdf.i.i$value==min(mdf.i.i$value),][1,])
	    }
	}
	if (min.and.max) {
	    #index of midpoint (halfway of x values)
	    int.midpoint.x <- as.integer(length(mdf.i.i$value) - as.integer( length(mdf.i.i$value)/2))
	    extr.err.df <- rbind(extr.err.df, data.frame(id=id, variable=variable
			    ,mean=mean(mdf.i.i$value),min=min(mdf.i.i$value)[1],max=max(mdf.i.i$value)[1],midpoint.value=mdf.i.i$value[int.midpoint.x]))
	}
      }  
  }
  extr.err.df
}

gg.blank <- function() {
    # creates a blank white plot as space holder for grid.arrange() plots
    qplot() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),legend.position="none",
    panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank())
    }


grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)

}

##____________________________________________________________________________________ Data Handling Utils
import.temp.df.from.csv <- function(system='windows',custompath='',header=TRUE, sep=','){

if (system=='windows') {path='F:/Data/06_Bibliothek/temp.csv'}
if (system=='linux') {path='/media/qnap/tardis/Data/06_Bibliothek/temp.csv'}
if (custompath!='') {path=custompath}
input <- read.table(path, sep = sep, header = header)  

}

rbind.csvfiledir.to.df <- function(file.dir,use.filename.props=TRUE, id.name.trim.range=c(0,-4)){
## Function to read an join csv files row-wise and give extra colums for '_' seperated properties
## Properties, Number and names of columns in csv must be the same
  file.list <- paste(file.dir, list.files(file.dir), sep='')				# create full path for all files in the dir
  names <- c()	
  
  
  for (file.name in list.files(file.dir)){						# trim filenames to used names accord trim range
      names <- append(names,substring(file.name, first=id.name.trim.range[1], last=nchar(file.name)+id.name.trim.range[2])) }
      
  ## Check if all files have the same number of '_' seperated properties   
  if (use.filename.props == TRUE){
    numprop.prev <- length(strsplit(list.files(file.dir)[1], '_')[[1]]) 
    for (file in list.files(file.dir)){ 
      numprop <- length(strsplit(file, '_')[[1]])
      if (numprop.prev != numprop) {
	  print(paste('#######  Error: Number of properties do not match for file in csv rbind:', file)) ; stop() }
      numprop.prev <- numprop
    }
  } 
  
  ## Assemble output
  df <- data.frame()
  i = 1
  for (file in file.list){							# for each file in the directory
    df.i <- read.table(file, sep = ',', header = TRUE)				# read to new temporary df.i
    if (use.filename.props == TRUE){
      props <- strsplit(names[i], '_')[[1]]					# create temp vector of the properties by the current filename 
      df.i$id <- props[1]								# put first property in this vector to the $id column 
      j <- 2
      for (prop in props[-1]) {							# for each additional j property add new column to df.i
	df.i[[paste('prop', + j, sep = '')]] <- props[j]
	j <- j+1}
    }  
    df <- rbind(df, df.i)
    i = i+1
    }
 df
}  
  
#print(rbind.csvfiledir.to.df(file.dir)[c(1:5),])

cbind.csvfiledir.to.df <- function(file.dir, use.firstCol.of.first.as.x=TRUE, file.name.for.names=TRUE){
## Function to read an join csv files row-wise and give extra colums for '_' seperated properties
## Properties, Number and names of 	columns in csv must be the same
  i <- 1
  for (file.name in paste(file.dir,list.files(file.dir,pattern="*.csv"),sep='')){						# trim filenames to used names accord trim range
    df.temp <- read.csv(file.name)
    cur.name <- tools::file_path_sans_ext(list.files(file.dir)[i])
    names(df.temp) <- c(paste(cur.name,'id',sep='.'),paste(cur.name,seq(1,length(names(df.temp))-1),sep='.'))
    df.temp.y <- data.frame(df.temp[,-1])
    names(df.temp.y) <-paste(cur.name,seq(1,length(names(df.temp))-1),sep='.')
    names(df.temp.y)[1] <- cur.name
    
    if (i == 1){ df <- df.temp[1];  df <- cbind(df,df.temp.y)}
    if (i > 1 & use.firstCol.of.first.as.x){ df <- cbind(df,df.temp.y)}
    if (i > 1 & !use.firstCol.of.first.as.x){ df <- cbind(df,df.temp) }
    i <- i+1
    }
   names(df)
   return(df)
  }
  
csvfiledir.to.list <- function(file.dir,custom.colnames.per.df=c('x','y1','y2')){
## Function to read an join csv files in a list of data.frames 
## Column names per dataframe are modified by csv name
  i <- 1
  df.list <- list()
  for (file.name in paste(file.dir,list.files(file.dir,pattern="*.csv"),sep='')){						# trim filenames to used names accord trim range
    df.temp <- read.csv(file.name)
    cur.name <- tools::file_path_sans_ext(list.files(file.dir)[i])
    names(df.temp) <- paste(custom.colnames.per.df, cur.name, sep='_')  
    df.list[[cur.name]] <- df.temp
    i <- i +1
    }
   return(df.list)
  }



##__________________________________________________________________________________________________________________________ Engineering Helpers
inch_to_mm <- function(inch) {inch*25.4}
sqrft_to_mm2 <- function(sqrft) {sqrft / 0.0000107639104167}
lbs_to_n <- function(lbs) {lbs * 4.4482216282509}
feet_to_mm <- function(feet){feet*304.8}
ksi_to_mpa <- function(ksi){ksi * 6.89475908677537}

K.to.Tensile.Modulus <- function(K,nu) {
#Convert compression Modulus in Tensile modulus
3*K*(1-2*nu)
}

Tensile.Modulus.to.K <- function(E,nu) {
#Convert tensile Modulus in compression modulus
E / (3 - 6 * nu)
}

G.from.E <- function(E,nu) {
#Convert tensile Modulus in compression modulus
E / (2 * (1+nu))
}


orthotropic.compliance.matrix <- function(E.1,E.2,E.3,G.12,G.23,G.31,nu.12,nu.23,nu.13,nu.21,nu.31){
    ## calculate isotropic compliance matrix from E,nu for isotropic materials (191117) wiki

    out <- matrix(nrow=6,ncol=6)
    out[1,] <- c(1/E.1  ,       -(nu.21)/E.2    ,   -(nu.31)/E.3    ,   0,0,0)
    out[2,] <- c(-(nu.12)   ,   1/E.2           ,   -nu.23/E.3          ,0,0,0)
    out[3,] <- c(-(nu.13)/E.3,  -nu.23/E.3,         1/E.3           ,0,0,0)
    out[4,] <- c(0,0,0,1/G.23,0,0)
    out[5,] <- c(0,0,0,0,1/G.31,0)
    out[6,] <- c(0,0,0,0,0,1/G.12)

    return(out)


}

isotropic.compliance.matrix <- function(E,nu) {
    ## calculate isotropic compliance matrix from E,nu for isotropic materials (190918)

    out <- matrix(nrow=6,ncol=6)
    out[1,] <- c(1,-nu,-nu,0,0,0)
    out[2,] <- c(-nu,1,-nu,0,0,0)
    out[3,] <- c(-nu,-nu,1,0,0,0)
    out[4,] <- c(0,0,0,2*(1+nu),0,0)
    out[5,] <- c(0,0,0,0,2*(1+nu),0)
    out[6,] <- c(0,0,0,0,0,2*(1+nu))

    return(out / E)

}

isotropic.stiffness.matrix <-function(E,nu) {

lambda <- (nu/(1-2*nu)) * (1/(1+nu)) * E
G <- E / (2 * (1+nu))

    out <- matrix(nrow=6,ncol=6)
    out[1,] <- c(lambda + 2*G, lambda , lambda,0,0,0)
    out[2,] <- c(lambda, lambda + 2*G , lambda,0,0,0)
    out[3,] <- c(lambda , lambda , lambda + 2*G,0,0,0)
    out[4,] <- c(0,0,0,G,0,0)
    out[5,] <- c(0,0,0,0,G,0)
    out[6,] <- c(0,0,0,0,0,G)

    return(out)
}

##__________________________________________________________________________________________________________________________ Dynamics

time_to_frequency_octave_p_minute <- function(t,f1,oct.per.min=1){
# Convert time t [s] to corresponding frequency at a octave/min sweep
exp(log(2)*t/60*oct.per.min + log(f1))}

sine.sweep.time.by.oct.min <- function (f.min=10,f.max=1000,oct.per.min=0.33) { 
  # initial value from TAS Test specification reverse engineering
  # time according octave speed input, so oct.per.min gives minutes as result
  t <- 1/oct.per.min * (log((2*pi*f.max)/(2*pi*f.min)))/log(2)
  return(t)
}  


no.cycles.sine.sweep.counting <- function (oct.per.min=oct.per.min <- c(0.33,1,4),f.min=10,f.max=1000,no.sweeps=432)
{
# oct.per.min: usually 1 but also seen 0.33 in TAS Program
# no.sweep  total test time devided by time per sweep
# f.min: lowest frequency in observed range
# f.max: highest frequency in observed range

#number of sweeps * sweep-range * ni
no.cycles <- no.sweeps*(f.max-f.min)*(1/((oct.per.min/60)*log(2)))

return(no.cycles)
}

sine.sweep.oct.min.from.time <- function (f.min=10,f.max=1000,t=20) { 
  # octave speed according time unitinput, so min gives minutes as result
  vs <- 1/t * (log((2*pi*f.max)/(2*pi*f.min)))/log(2)
  return(vs)
}

##__________________________________________________________________________________________________________________________ DoE Helpers
modif.two.level.for.CCD <- function(factor.list, reduce.corners = TRUE){
## function which modifies the corners of a 2 level design in the way that 
## the star points do not exceed the intial corners. Usefull if only extreme simulateable
## values are known. 
# reduce.corners is default for modification, FALSE if star ist calculated from kept corners
  df <- as.data.frame(lapply(factor.list,min))
  df <- rbind(df,as.data.frame(lapply(factor.list,max)))
  df <- rbind(df,as.data.frame(lapply(factor.list,mean)))
  rownames(df) <- c('min','max','mean')

  alpha <- (2 ^ length(factor.list)) ^ 0.25   # alpha for rotateable

  df['star.min',] <- df['mean',] + ((df['mean',] - df['max',]) * alpha)
  df['star.max',] <- df['mean',] - ((df['mean',] - df['max',]) * alpha)

  df['two.lvl.min',] <- -1 * (((df['min',] - df['mean',]) / alpha) - df['mean',])
  df['two.lvl.max',] <- -1 * (((df['max',] - df['mean',]) / alpha) - df['mean',])


  if (reduce.corners == FALSE){return(as.list(df[c('star.min','star.max'),]))}
  if (reduce.corners == TRUE){return(as.list(df[c('two.lvl.min','two.lvl.max'),]))}
}


##____________________________________________________________________________________________________________________________ Critria Script functions
tabular.matrix.properties = function(
    fiber.l.d.ratio=30,
    m.density=1130,matrix.modulus=1500,
    df.matrix.alpha=data.frame(temp=c(-50,23,50,160,200),alpha=c(90,130,180,230,280))){
   
   
   #fiber.l.d.ratio = fiber l/D ratio [-]
   #m.density = Density of Matrix in kg/m^3
   #df.matrix.alpha = dataframe of alpha over temperature [temp in degC alpha in 10-6/K]
   #matrix.modulus = Matrix Modulus for Halpin-Tsai in Anisotropic APDL Export File
 
    ## linear fit of Matrix Alpha
 
  fit.m.alpha <- lm(df.matrix.alpha$alpha ~ df.matrix.alpha$temp +I(df.matrix.alpha$temp) )
 fit.m.alpha.fun <<- function (x){ fit.m.alpha[["coefficients"]][1] + fit.m.alpha[["coefficients"]][2] * x}
 fiber.l.d.ratio <<- fiber.l.d.ratio
 m.density <<-m.density
 # Glass Properties
 m.glass <- 73000
 poisson.glass <- 0.22
 glass.density <- 2600  #kg/m^3
 alpha.glass <- '5.4e-006'  # '' signs just to ensure that double zero e-006 is kept for right converse format  
 matrix.modulus <<- matrix.modulus
 colnames(s.s.short.matrix) <- c("temp", "strain", "stress")
 s.s.short.matrix$strain <- s.s.short.matrix$strain / 100
 s.s.short.matrix.true <- strain.stress.convert(s.s.short.matrix$strain,s.s.short.matrix$stress,'eng->true', s.s.short.matrix$temp)      # Calc True Strain, True stress
 temp.d.scalar.res.matrix <- lin.el.extract(s.s.short.matrix.true$strain, s.s.short.matrix.true$stress, s.s.short.matrix.true$temp)[,c(1,4)]
 temp.d.scalar.res.matrix <- temp.d.scalar.res.matrix[!duplicated(temp.d.scalar.res.matrix$e.el),]
 temp.d.scalar.res.matrix <- rbind(temp.d.scalar.res.matrix,c(temp=1000,tail(temp.d.scalar.res.matrix[2],1)))
 temp.d.scalar.res.matrix <<- rbind(c(temp=-270,head(temp.d.scalar.res.matrix[2],1)),temp.d.scalar.res.matrix)
 m.matrix.fun <<- approxfun(temp.d.scalar.res.matrix$temp,temp.d.scalar.res.matrix$e.el)
 return(list(temp.d.scalar.res.matrix=temp.d.scalar.res.matrix,m.matrix.fun=m.matrix.fun))
} 

func.s.s.isochr.input.treatment = function(s.s.isochr) {
  colnames(s.s.isochr) <- c("temp", "stress","time", "strain")
  s.s.isochr_org <<- s.s.isochr
  s.s.isochr <- s.s.isochr[s.s.isochr$stress > 0,]
  s.s.isochr$time <- s.s.isochr$time * 60 * 60  # Convert hours to seconds
  s.s.isochr$strain <- s.s.isochr$strain / 100  # Convert percent to mm/mm
  s.s.isochr$modulus <- s.s.isochr$stress / s.s.isochr$strain
  s.s.isochr <<- s.s.isochr
}

func.s.s.short.input.treatment <- function(s.s.short){
    #gives correct colnames "temp", "strain", "stress", convert strain in mm/mm, convert to true values
    s.s.short_org <- s.s.short 
    colnames(s.s.short) <- c("temp", "strain", "stress")
    s.s.short$strain <- s.s.short$strain / 100                          # convert strain unit from percent to mm/mm strain
    s.s.short <<- s.s.short
    s.s.short.true <- strain.stress.convert(s.s.short$strain,s.s.short$stress,'eng->true', s.s.short$temp) # Calc True Strain, True stress
    return(s.s.short.true)
}

fun.approx.s.s.hort <- function(s.s.df=s.s.short.true.A.stiff,strain=NA,stress=NA,temp) {
    #function to get either stress from strain or vise versa from approxfun of temp,stress,strain dataframe
    out <- c()
        if(is.na(stress[1])) {
            for (i in 1:length(strain)){     
                out <- c(out,approxfun(s.s.df[s.s.df$temp==temp[i],]$strain, s.s.df[s.s.df$temp==temp[i],]$stress)(strain[i]))
            }
            return(out)
        }
        if(is.na(strain[1])) {
            for (i in 1:length(stress)){     
                out <- c(out,approxfun(s.s.df[s.s.df$temp==temp[i],]$stress, s.s.df[s.s.df$temp==temp[i],]$strain)(stress[i]))
            }
            return(out)
        }

        
        if(is.na(stress[1]) && is.na(strain[1])) {print('give either strain or stress as X')}
}

get.room.temp <- function(s.s.short.true){
    # function to get closest temperature value of the measured data to 23, to ensure that this temperature exists in later processing
    return(unique(s.s.short.true$temp )[which(abs(unique(s.s.short.true$temp )-23)==min(abs(unique(s.s.short.true$temp )-23)))] )
}

Tg.from.eng.dens <- function(s.s.short.true){
    # function to get T.g by maximum value of eng.dens curve
    T.g <- optimize(integral.multicurve(s.s.short.true)$fun,interval=c(min(s.s.short.true$temp),max(s.s.short.true$temp)),maximum=TRUE)$maximum
    s.s.short.true$temp
    return(T.g)
}

A.w.s_factor_from_list <- function(fiber.w.perc,polymer.type){

    # function to get column wise list of reduction factors out of #korte2019 based on Fiber and polymer type
    load.types.min.max <- c("sig.one.time.max","sig.one.time.min","sig.multiple.max","sig.multiple.min","sig.relaxation.max","sig.relaxation.min",
    "sig.creep.max","sig.creep.min","sig.cyclic.max","sig.cyclic.min")
    load.types <- c("sig.one.time","sig.multiple","sig.relaxation","sig.creep","sig.cyclic")

    if (fiber.w.perc > 0 && polymer.type == 1) {
    A.w.s <- data.frame(min=c(0.55,0.45,0.35,0.45,0.25),max=c(0.7,0.55,0.45,0.55,0.25))}  # GF Semicrystalline

    if (fiber.w.perc == 0 && polymer.type == 1)  {
    A.w.s <- data.frame(min=c(0.8,0.6,0.5,0.6,0.2),max=c(1,0.8,0.6,0.7,0.3))}  # Semicrystalline

    if (fiber.w.perc > 0 && polymer.type == 2)  {
    A.w.s <- data.frame(min=c(0.55,0.45,0.35,0.45,0.17),max=c(0.7,0.55,0.45,0.55,0.17))}  # GF Amorphous

    if (fiber.w.perc == 0 && polymer.type == 2)  {
    A.w.s <- data.frame(min=c(0.65,0.5,0.4,0.5,0.17),max=c(0.8,0.65,0.5,0.6,0.2))}  # Amorphous

    if (polymer.type == 3)        {
    A.w.s <- data.frame(min=c(0.65,0.5,0.4,0.5,0.17),max=c(0.8,0.65,0.5,0.6,0.2))}  # PS (!No value given <- Value from Amorphous taken)

    
    A.w.s.min <- data.frame(t(A.w.s$min))
    A.w.s.max <- data.frame(t(A.w.s$max))
    colnames(A.w.s.min) <- colnames(A.w.s.max) <- load.types
    A.w.s <- data.frame(t(c(rbind(A.w.s$max,A.w.s$min)))) #create one alternating vector from A.ws max,min,max...
    colnames(A.w.s) <- load.types.min.max

    return(list(min.and.max=A.w.s,min=A.w.s.min,max=A.w.s.max))
}

safety.factors <- function(probability="high",severety="low",breaktype="ductile",anisofactor=1) {
    # function to retrieve safety factors according #korte2019p194
    if (severety == "high") {sev=1} ; if (severety == "mid") {sev=2}; if (severety == "low") {sev=3}
    if (probability == "high") {prob=1} ; if (probability == "mid") {prob=2}; if (probability == "low") {prob=3}

    df.ductile <- data.frame(high.severety=c(1.5,1.35),mid.severety=c(1.4,1.25),low.severety=c(1.3,1.2))
    df.brittle <- data.frame(high.severety=c(2,1.8),mid.severety=c(1.85,1.7),low.severety=c(1.75,1.6))

    rownames(df.brittle) <- rownames(df.ductile) <- c("high.probability","low.probability")
    out <- c()
    for (i in seq(1,length(breaktype),1)) {  # look for breaktype and assign corresponding safety factor
        if (breaktype[i] == "ductile"){out <- append(out, df.ductile[prob,sev])}
        if (breaktype[i] == "brittle"){out <- append(out, df.brittle[prob,sev])}
    }
    return(out)
}

eng.dens.el <- function(s.s.short.true)  {
    # function to extract only elastic part of engergy density of stress-strain curve
    # requires lin.el.extract function
    
    eng.dens.el <- data.frame(temp=lin.el.extract(s.s.short.true$strain, s.s.short.true$stress, s.s.short.true$temp)[1], eng.dens.el = 0.5* lin.el.extract(s.s.short.true$strain, s.s.short.true$stress, s.s.short.true$temp)[3] * lin.el.extract(s.s.short.true$strain, s.s.short.true$stress, s.s.short.true$temp)[5])
    # df of temperatures and strain * stress * 0.5 from lin.el.extract function
    
    out.eng.el <- c()
    for (temp in unique(eng.dens.el$temp)) { # for each temperature
        out.eng.el <- c(out.eng.el,max(eng.dens.el[eng.dens.el$temp == temp,][2])) #get rows for temperature and search for max energy
    }
    eng.dens.el <- data.frame(temp=unique(eng.dens.el$temp),eng.dens.el=out.eng.el)
    return(eng.dens.el)

}

ggplot_theme <-theme_bw() + theme(plot.background = element_blank() )+ theme(panel.border= element_blank())+  theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))

energy.density.over.temp <- function(s.s.short.true,brittle.smaller.than.ratio=2){
    #function to extract energy density data.frame
    
    eng.dens <- data.frame(integral.multicurve(s.s.short.true)$table,eng.dens.el(s.s.short.true)[2])
    eng.dens <- data.frame(eng.dens,eng.dens[2]-eng.dens[3])
    eng.dens <- data.frame(eng.dens,eng.dens[4]/eng.dens[3])
    colnames(eng.dens) <- c("temp","total","el","pl","ratio")
    eng.dens$breaktype = "ductile"
    if (min(eng.dens$ratio)<= brittle.smaller.than.ratio) {
    eng.dens[eng.dens$ratio <= brittle.smaller.than.ratio,]$breaktype <- "brittle"}
    #eng.dens[eng.dens$ratio > 3.3,]$breaktype <- "ductile"


    plot.df <- melt(eng.dens,id="temp")[melt(eng.dens,id="temp")[2] != "breaktype",]; plot.df$value <- as.numeric(plot.df$value)
    gg.eng.dens <- ggplot(data=plot.df,aes(x=temp,y=value,group=variable,color=variable)) + geom_line(aes(),size=2,linetype="solid") +
    #scale_color_manual(name='Energy type',labels=c('total','elastic','plastic','ratio el/pl'),values = rev(brewer.pal(4,"Spectral")))+ 
    scale_colour_manual(name='temp',values = rev(hcl(h = seq(1,230,length.out=4), c = 100, l = 75)))+
    geom_point(aes())+
    scale_y_continuous(breaks = seq(0, max(plot.df$value*1.5),0.25), limits= c(0,max(plot.df$value*1.1)))+ 
    scale_x_continuous(breaks = seq(0, max(s.s.short.true$temp*1.5),20), limits= c(0,max(s.s.short.true$temp*1.1)))+
    geom_hline(yintercept=brittle.smaller.than.ratio,linetype=4,color='red')

    gg.eng.dens <- gg.eng.dens + ggplot_theme +labs(title=paste('max. Energy Density over temperature'), x='Temperature [degC]', y='Energy Density [MPa]')
    print(gg.eng.dens)
    colnames(eng.dens) <- c('temp','w.total','w.el','w.pl','w.ratio','breaktype')
    return(eng.dens)
}


neuber.plot.energies <- function(s.s.short.true,room.temp,str.mat.name) {
    #function to plot stress strain curve with ggplot and overlay neuber curves
    
    #find closest temperature to 23

    # quadratic energy densities for min,rt,max temperatures ! only for plotting of neuber
    w.min.temp <-  max(s.s.short.true[s.s.short.true$temp == min(s.s.short.true$temp),]$stress * s.s.short.true[s.s.short.true$temp == min(s.s.short.true$temp),]$strain)
    w.rt <-        max(s.s.short.true[s.s.short.true$temp == room.temp,]$stress * s.s.short.true[s.s.short.true$temp == room.temp,]$strain)
    w.max.temp <-  max(s.s.short.true[s.s.short.true$temp == max(s.s.short.true$temp),]$stress * s.s.short.true[s.s.short.true$temp == max(s.s.short.true$temp),]$strain)

    # Plot Stress-Strain
    gg <- ggplot(data=s.s.short.true,aes(x=strain,y=stress,color=as.factor(temp))) + geom_point()+geom_line(aes(),size=1) + scale_colour_manual(name='temp',values = rev(hcl(h = seq(1,230,length.out=length(unique(s.s.short.true$temp))), c = 100, l = 75)))
    gg <- gg + ggplot_theme+ labs(title=paste(str.mat.name,'- Stress-Strain True Values'), x='true strain [mm/mm]', y='stress [MPa]') +
    #Add Neuber Lines
    geom_line(data=data.frame(x=seq(0,max(s.s.short.true$strain),length.out=50),y=w.min.temp/seq(0,max(s.s.short.true$strain),length.out=50)),aes(x=x,y=y,group='iso-engergy min'),color='black',linetype=5)+ 
    geom_line(data=data.frame(x=seq(0,max(s.s.short.true$strain),length.out=50),y=w.rt/seq(0,max(s.s.short.true$strain),length.out=50)),aes(x=x,y=y,group='iso-engergy mean'),color='black',linetype=6)+ 
    geom_line(data=data.frame(x=seq(0,max(s.s.short.true$strain),length.out=50),y=w.max.temp/seq(0,max(s.s.short.true$strain),length.out=50)),aes(x=x,y=y,group='iso-engergy max'),color='black',linetype=5)+ 
    scale_y_continuous(breaks = seq(0, max(s.s.short.true$stress*1.5),10), limits= c(0,max(s.s.short.true$stress*1.1)))+ 
    scale_x_continuous(breaks = seq(0, max(s.s.short.true$strain*1.5),0.01), limits= c(0,max(s.s.short.true$strain*1.1)))
    print(gg)
    gg.stress.strain <<- gg

}

sigma.max.UTS.over.temp <- function(s.s.short.true){
    #function to get maximum stress values per temperature
    value <- c()
    temp <- unique(s.s.short.true$temp)
    for (i in unique(s.s.short.true$temp)){
        value <- append(value,max(s.s.short.true[s.s.short.true$temp == i,]$stress))
    }

    return(data.frame(temp=temp,sigma.uts=value))
}

ep.from.sig.fun.list <- function(s.s.short.true){
    # function to get inverse approxfuns from stress-strain for all temperatures in a list with temp as names
    ep.from.sig.fun.list <- list()
    for (i in unique(s.s.short.true$temp)){
    strain <- s.s.short.true[s.s.short.true$temp == i,]$strain
    stress <- s.s.short.true[s.s.short.true$temp == i,]$stress
    ep.from.sig.fun.list <- append(ep.from.sig.fun.list,approxfun(stress,strain))
    }

    names(ep.from.sig.fun.list) <- unique(s.s.short.true$temp)
    return(ep.from.sig.fun.list)
}

stress.strain.limits.over.temp <- function(df.temp.sig.uts=df.eng.dens$sigma.uts,df.s.s.short, df.factors=df.reductionf.temp.load){

    # function to apply a matrix of reduction factors (loadtypes x temp) on a vector of stress UTS values to get allowable stress per load and temp
    # second function is to use the allowable stress to retrieve corresponding strain from stress-strain curves
    # dependencies:
    #   - ep.from.sig.fun.list() to retrieve approx functions for reverse-engineering
    
    df.sig.allow <- df.temp.sig.uts * df.factors[,-1] #create stress limits data.frame, remove temp from factors since temp column already in sig.uts
    ep.from.sig.fun.list <- ep.from.sig.fun.list(df.s.s.short)

    df.ep.allow <- data.frame(temp=unique(df.eng.dens$temp))
    for (load.type in colnames(df.sig.allow)){
        sig <- df.sig.allow[,load.type]
        ep.curr.load.type <- c()

            for (i in seq(1,length(unique(df.eng.dens$temp)))){
                ep.curr.load.type <- append(ep.curr.load.type, ep.from.sig.fun.list[[i]](sig[i])) }

        df.ep.allow <- cbind(df.ep.allow,ep.curr.load.type)
    }

    df.sig.allow <- cbind(temp=df.eng.dens$temp,df.sig.allow) 

    # create colnames for 
    load.types <- c("ep.one.time","ep.multiple","ep.relaxation","ep.creep","ep.cyclic")
    colnames(df.ep.allow) <- c("temp",load.types)

    return(list(df.ep.allow=df.ep.allow,df.sig.allow=df.sig.allow))
}

fit.function.limits <- function(df.allowable,degree=6,xname='x',str.nmbr.digits=5,x.col=1,y.col=2){
    
    # function to create a arbituary polynomial fit on the first two columns of a given data.frame
    # returns the function as mathematical string with xname for the x
    
    x.vals <- df.allowable[[x.col]]                                                     # store temperatures as x values
    y.vals <- df.allowable[[y.col]]                                                      # second values which indicates one.time.loading for y values
    fit.string <- 'x.vals'                                                          # initiate string which is used as polynomial formular in lm fit
    max.limit.function <- 'function(x){ max.limit.fit$coefficients[[1]]'            # initiate function fit
    if (degree==1) {fit.string <- paste(fit.string,'+I(x.vals^',1,')',sep='')}
    if (degree>1){
    for (i in seq(1,degree-1,1)){                                                   # assemble fit and function string according set degree
        fit.string <- paste(fit.string,'+I(x.vals^',i+1,')',sep='')
    } }
    fit.string <- paste('lm(y.vals ~ ',fit.string,")",sep='')                       # finialize fit string
    
    if (degree==1) {        max.limit.function <- paste(max.limit.function,
        '+max.limit.fit$coefficients[[',1+1,']]*x^',1,sep='')}
    
    if (degree>1){
    for (i in seq(1,degree,1)){  
        max.limit.function <- paste(max.limit.function,
        '+max.limit.fit$coefficients[[',i+1,']]*x^',i,sep='')
    }}

    max.limit.function <- paste(max.limit.function,'}')                             # finialize function string

    max.limit.fit <- eval(parse(text=fit.string))                                   # execute fit according fit.string        
    max.limit.fun <- eval(parse(text=max.limit.function))                           # assemble fit function
    fit.results <- max.limit.fun(x.vals)                                            # apply fitted function on x values
    
    coeff <- as.numeric(max.limit.fit$coefficients)                                 # store coefficient resutls of fit in short name vector
    math.string <- as.character(round(coeff[1],str.nmbr.digits))
   
   if (degree==1) {  cur.coeff <- round(coeff[1+1],str.nmbr.digits)
        if (cur.coeff >= 0){
        math.string <- paste(math.string,'+',cur.coeff,'*',xname,'^',1,sep='')}
        if (cur.coeff < 0){
        math.string <- paste(math.string,cur.coeff,'*',xname,'^',1,sep='')}}
   if (degree>1){
    for (i in seq(1,degree,1)){                                                   # assemble fit and function string according set degree
        cur.coeff <- round(coeff[i+1],str.nmbr.digits)
        if (cur.coeff >= 0){
        math.string <- paste(math.string,'+',cur.coeff,'*',xname,'^',i,sep='')}
        if (cur.coeff < 0){
        math.string <- paste(math.string,cur.coeff,'*',xname,'^',i,sep='')}
    } }
    
    return(list(max.limit.fit=max.limit.fit,max.limit.fun=max.limit.fun,fit.results=fit.results,math.string=math.string))

}

plot.limit.over.temp <- function(df.allowable, ylab='Allowable Strain [mm/mm]'){
plot.df <- melt(df.allowable,id = "temp")
ggplot(data=plot.df,aes()) + ggplot_theme +
geom_line(aes(x=temp,y=value,color=variable))+
geom_point(aes(x=temp,y=value,color=variable))+
labs(title=paste('Allowable values over temperature and load type'), x='Temperature [degC]', y=ylab)+
scale_color_discrete(name="Criterion")
}


plot.fatique.results <- function(temp.d.scalar.res,sig.bk.fun,sig.bk.gerber.fun,color.scale.temp,df.sig.allow,df.ep.allow){

    sig.m.plot <- c(0,seq(-1*max(df.sig.allow[["One Time"]]),max(df.sig.allow[["One Time"]]),length.out=100))


    df.bk.fat <- expand.grid(c(10^(0:7)),sig.m.plot,sort(unique(df.fatique$temp)))                      # create grid coordinates from N,sig.m and temperature
    colnames(df.bk.fat) <- c('N','sig.m','temp')                                                        # clean column names
    df.bk.fat$lin.bk.stress <- sig.bk.fun$F1(sig.m=df.bk.fat$sig.m,temp=df.bk.fat$temp,N=df.bk.fat$N)          # calculate (linear) sig.bk.fun, which spans from wk to uts(by ak)
    df.bk.fat$trim.bk.stress <- pmin(df.bk.fat$lin.bk.stress, sig.sk.temp.fun(df.bk.fat$temp)-abs(df.bk.fat$sig.m)) # trimm linear bk by the equal side triangle of sig.sk(one time limit)
    df.bk.fat$gerber.stress <- sig.bk.gerber.fun(sig.m=df.bk.fat$sig.m,temp=df.bk.fat$temp,N=df.bk.fat$N)      # calculate mean stress dependency by gerber approach
    df.bk.fat <- df.bk.fat[df.bk.fat$trim.bk.stress >= 0,]
    df.bk.fat$ep.m <- fun.approx.s.s.hort(stress=df.bk.fat$sig.m,temp=df.bk.fat$temp)
    df.bk.fat$trim.bk.strain <- fun.approx.s.s.hort(stress=df.bk.fat$lin.bk.stress,temp=df.bk.fat$temp)
    df.bk.fat$gerber.strain <- fun.approx.s.s.hort(stress=df.bk.fat$gerber.stress,temp=df.bk.fat$temp)


    gg.mean.stress <<- ggplot(data=df.bk.fat)+geom_vline(xintercept=0)+ theme_bw() +
    geom_line(aes(x=sig.m,y=trim.bk.stress,group=temp,color=as.factor(temp)))+
    facet_grid(~N)+ scale_colour_manual(name='temp [degC]',values = color.scale.temp)+
    geom_line(aes(x=sig.m,y=gerber.stress,group=temp,color=as.factor(temp)),alpha=0.2)+
    scale_y_continuous(limits=c(0,max(df.sig.allow[["One Time"]])*1.1),name="Allowed Amplitude Stress [MPa]")+
    scale_x_continuous(limits=c(0,max(df.sig.allow[["One Time"]])*1.1),name="Allowed Mean Stress [MPa]")

    gg.mean.strain <<- 
    ggplot(data=df.bk.fat)+geom_vline(xintercept=0)+ theme_bw() +
    geom_line(aes(x=ep.m,y=trim.bk.strain,group=temp,color=as.factor(temp)))+
    facet_grid(~N)+ scale_colour_manual(name='temp [degC]',values = color.scale.temp)+
    geom_line(aes(x=ep.m,y=gerber.strain,group=temp,color=as.factor(temp)),alpha=0.2)+
    #geom_line(aes(x=ep.m,y=gerber.strain.fit,group=temp,color=as.factor(temp)),alpha=0.2)+
    scale_y_continuous(name="Allowed Amplitude Strain [mm/mm]")+
    scale_x_continuous(name="Allowed Mean Strain [mm/mm]")+
    coord_cartesian(xlim=c(0,max(df.ep.allow[["One Time"]])*1.1),ylim=c(0,max(df.ep.allow[["One Time"]])*1.1))+
    geom_hline(yintercept = max(temp.d.scalar.res$crit.strain),linetype="dotted",color=color.scale.temp[length(color.scale.temp)])+
    geom_hline(yintercept = min(temp.d.scalar.res$crit.strain),linetype="dotted",color=color.scale.temp[1])
    #stat_function(fun = function(x) (max(df.ep.allow[["Cyclic 10^4"]])-1*x))



    gg.strength.fat <<- ggplot(data=df.bk.fat[df.bk.fat$sig.m==0,])+ theme_bw() + ggplot_theme+
    geom_line(aes(x=N,y=trim.bk.stress,group=temp,color=as.factor(temp)))+
    geom_point(aes(x=N,y=trim.bk.stress,group=temp,color=as.factor(temp)))+
    scale_colour_manual(name='temp [degC]',values = color.scale.temp)+
    scale_y_continuous(limits=c(0,max(df.sig.allow[["One Time"]])*1.1))+
    scale_x_continuous(breaks=c(10^(0:7)),labels=c(0:7),trans="log10")+
    labs(x="Number of Cycles [log10]",y="Strength [MPa]")

    gg.strain.fat <<- ggplot(data=df.bk.fat[df.bk.fat$sig.m==0,])+ theme_bw() + ggplot_theme+
    #geom_line(aes(x=N,y=trim.bk.strain,group=temp,color=as.factor(temp)))+
    geom_line(aes(x=N,y=trim.bk.strain,group=temp,color=as.factor(temp)))+
    geom_point(aes(x=N,y=trim.bk.strain,group=temp,color=as.factor(temp)))+
    #geom_line(aes(x=N,y=gerber.strain.fit,group=temp,color=as.factor(temp)))+
    scale_colour_manual(name='temp [degC]',values = color.scale.temp)+
    scale_y_continuous(breaks=seq(0,max(df.ep.allow[["One Time"]])*1.1,0.005),limits=c(0,max(df.ep.allow[["One Time"]])*1.1))+
    scale_x_continuous(breaks=c(10^(0:7)),labels=c(0:7),trans="log10")+
    labs(x="Number of Cycles [log10]",y="Allow. Strain BK [mm/mm]")+
    geom_hline(yintercept = max(temp.d.scalar.res$crit.strain),linetype="dotted",color=color.scale.temp[length(color.scale.temp)])+
    geom_hline(yintercept = min(temp.d.scalar.res$crit.strain),linetype="dotted",color=color.scale.temp[1])
    
    return(list(df.bk.fat=df.bk.fat,gg.strain.fat=gg.strain.fat,gg.mean.strain=gg.mean.strain,gg.mean.stress=gg.mean.stress))
}    

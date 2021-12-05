rm(list=ls(all=TRUE))

library('minpack.lm');library(reshape2);library(ggplot2);library('rgl');library(RColorBrewer)
require(gridExtra);require(grid);require(openxlsx)




## Things to add
# Interpolation between temperature curves for stress-strain & creep
# WB Command line and plot of conic criterion
# WB Command line and plot of tsai wu criterion
# Results in a Report
# Remove elastic strains from creep data
# Make all for multi temperatures
# create BISO fitting
# Change everything to Kelvin


## 0 Set working folder + Plot settings
setwd('')
source('criteria.lib_210917.R')


plot.fatique.curves <- FALSE



## 1. Material General information from Datasheet
grade.name                  <- "Material Grade name"
supplier.name               <- "Supplier Name"
polymer.abbr                <- "PA66"
E.Datasheet                 <- 7200 #MPa
strain.at.break.datasheet   <- 5    #%
stress.at.break.datasheet   <- 130 #MPa
k.RT                        <- 0.313 #W/m*K
c.RT                        <- 2036 #J/KG*K


## 1.1____Matermial-Name
str.mat.name <-  'xx_PA_Generic' # Name to be displayed in plots
polymer.type <- 1               # semicr = 1, amorphous = 2, PS =3, TPE/Rubbers=4

## 1.2__ Skalar Material Properties
T.g <- 39        # If not known, please insert 'NULL' to create an approximation by max. critical energy (sec 'Energy density under ref. Curve and Tg approx')
rho <- 1360        # Material Density kg/m^3
clte.flow.dir <- 27e-6
clte.perp.flow.dir <- 86e-6
fiber.w.perc <- 30             # Mass fraction of fiber in [%]
m.ratio <- 1.3                    # Ratio of Compression Strength to Tensile strength for isotropic criteria, default 1.3


## 2. Matrial Stress Strain Curves
## 2.2____Material Stress - Strain Curves over Tempereature

s.s.short.file.path <- 'Stress-strain.csv'  #("temp" [degC], "strain" [%], "stress" [MPa])#("temp" [degC], "strain" [%], "stress" [MPa])


## 2.3__ Creep Curve for Long term loads
create.creep.critera <- FALSE #TRUE/FALSE
evaluate.creep <- TRUE    #TRUE/FALSE
s.s.isochr.file.path <- 'isocr-stress-strain.csv'  #("temp" [degC], "stress" [MPa],"time" [h],"strain"[%])



######___________________ Misc Settings


### 3. Strenght Settings
## 3.1_________ Critical Strain Settings
brittle.smaller.than.ratio <- 2 # Ratio between w.el/w.pl [2-3.3] #korte2019p194 higher values means more brittle regions/ more conservative

## 4. ___ S - Safety Factor settings (high,mid,low) #korte2019p194
probability <- "high"
severety <- "low"


## 5. ____ Output of Converse Material settings
   #out.vec.conv.mat <- append(out.vec.conv.mat,(paste(temp.d.scalar.res$temp[i],temp.d.scalar.res$m.modulus[i],
    #temp.d.scalar.res$poisson[i],m.density/1000,alpha,m.glass,poisson.glass,glass.density/1000,

## 6. _____ Extrapolation of data settings

## 6.1 Short Term Stress Strain
s.s.point.no.reduced <- 20  # no of points per temperature to be used from stress-strain curve
temp.max <- 125             #max. temperature to evaluate from extrapolation
temp.min <- -40             #min. temperature to evaluate from extrapolation


## 6.2 Isochronous Stress Strain
s.s.iscrr.temp.to.extrapl <- 125  # only to higher temp. allowed so far

## 7 ___ Composite Matrix properties
ar <- 20                    # Fiber Aspect Ratio [-]  (10-30)
E.f <- 72000                # Fiber Modulus [MPa]
nu.f <- 0.22                # Fiber Poisson Ratio [-]
rho.f <- 2550               # Fiber Density [kg/m^3]
cte.f <- 5.1e-6             # Fiber Coefficient of thermal expansion [1/K]
k.f <- mean(c(1.2,1.35))    # Fiber Coefficient of thermal expansion [1/K]
c.f <- mean(c(800,805))    # Fiber Coefficient of thermal expansion [1/K]
fiber.orientation <- 0.8


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


df.general.info <- data.frame(name=NA,value=NA,unit=NA)
df.general.info <- rbind(df.general.info,c("Grade Name:",grade.name,""))
df.general.info <- df.general.info[-1,]
df.general.info <- rbind(df.general.info,c("Producer:",supplier.name,""))
if (fiber.w.perc == 0) {df.general.info <- rbind(df.general.info,c("Polymer Type:",paste(polymer.abbr,sep=""),"[%]"))} 
if (fiber.w.perc > 0) {df.general.info <- rbind(df.general.info,c("Polymer Type:",paste(polymer.abbr,"-GF",fiber.w.perc,sep=""),""))} 
df.general.info <- rbind(df.general.info,c("E-Modulus:",E.Datasheet,"[MPa]"))
df.general.info <- rbind(df.general.info,c("Density:",rho,"[kg/m^3]"))
df.general.info <- rbind(df.general.info,c("CLTE n/p:",paste(clte.flow.dir*10^6,"/",clte.perp.flow.dir*10^6,sep=""),"[10^-6/K]"))
df.general.info <- rbind(df.general.info,c("Strain at break RT:",strain.at.break.datasheet,"[%]"))
df.general.info <- rbind(df.general.info,c("Stress at break RT:",stress.at.break.datasheet,"[MPa]"))
    


# Create Output directory
dir.create(str.mat.name)
dir.create(paste(str.mat.name,"/Fiber-Files",sep=""))
dir.create(paste(str.mat.name,"/graphs",sep=""))


# Create copy of stress-strain and read
dir.create(paste(str.mat.name,"/raw-input-copy",sep=""))
file.copy(s.s.short.file.path,paste(str.mat.name,"/raw-input-copy/Stress-Strain.csv",sep=""),overwrite=TRUE)

s.s.short.file.path <- paste(str.mat.name,"/raw-input-copy/Stress-Strain.csv",sep="")
s.s.short <- read.table(s.s.short.file.path, sep = ',', header = TRUE) #("temp" [degC], "strain" [%], "stress" [MPa])#("temp" [degC], "strain" [%], "stress" [MPa])

# Create copy of isochr stress strain and read

if (evaluate.creep) {
    file.copy(s.s.isochr.file.path,paste(str.mat.name,"/raw-input-copy/Isochr-Stress-Strain.csv",sep=""),overwrite=TRUE)
    s.s.short.file.path <-paste(str.mat.name,"/raw-input-copy/Isochr-Stress-Strain.csv",sep="")
    s.s.isochr <- read.table(s.s.isochr.file.path, sep = ",", header = TRUE)
}

## ____ Reduce point in s.s.short


#reduce.s.s.points.raw <- function (s.s.short,s.s.point.no.reduced)

# preconditioning of s.s.short
colnames(s.s.short) <- c("temp","strain","stress")
s.s.short <- s.s.short[s.s.short$strain >=0,]  #remove negative strain values
s.s.short <- s.s.short[s.s.short$stress >=0,]  #remove negative stress values
s.s.short <- s.s.short[!(s.s.short$strain > 0 & s.s.short$strain <0.05e-2),]  #remove small strain values

s.s.short.red <- data.frame()

for (temp in unique(s.s.short$temp)) {  # for all temperatures
    #temp <- unique(s.s.short$temp)[1]

    if (length(s.s.short[s.s.short$temp == temp,][,1]) > s.s.point.no.reduced){
   
    s.s.short.temp <- s.s.short[s.s.short$temp == temp,]

    remove.rows <- as.integer(seq(1,length(s.s.short.temp[,1]),length.out=length(s.s.short.temp[,1])-s.s.point.no.reduced)) # create a sequence of the length of the items which are too much 
    if (remove.rows[1] == 1) {remove.rows[1] <- 2}  # avoid to remove first item because this is usually the 0,0 which is needed
    remove.rows <- sort(unique(remove.rows))        # delete duplicates form integer operation
    s.s.short.temp <- s.s.short.temp[-c(remove.rows),] # remove items
    s.s.short.red <- rbind(s.s.short.red,s.s.short.temp[s.s.short.temp$index==1,][,1:3]) # combine to overall df with all temperatures
    }
   
     if (length(s.s.short[s.s.short$temp == temp,][,1]) <= s.s.point.no.reduced){
   
    s.s.short.temp <- s.s.short[s.s.short$temp == temp,]
    s.s.short.red <- rbind(s.s.short.red,s.s.short.temp)
    }  


}


s.s.short <- s.s.short.red


##______ interpolate/extrapolate stress strain curves
#s.s.short <- s.s.short[s.s.short$strain > 0.3 | (s.s.short$strain == 0),]    ## remove all strain values below 0.5% since they are assumed to be part of linear viscoelastic region and may corrupt the modulus results
temp.interp.user <- c()  #extra temperatures to inter/extra polate -- not used yet
temp.interp <- c()

auto.fill.missing.temp <- TRUE
if (auto.fill.missing.temp == TRUE) {
    temp.add.extreme <- c()                     #vector of extreme temperatures to add
    temp.add <- sort(unique(s.s.short$temp))    #vector of intermediate temperatures to add if number less than 6

    #if supplied range not within the requested extreme temperatures, add extreme temperatures
    if (min(unique(s.s.short$temp)) > temp.min) {temp.add.extreme <- sort(c(temp.add.extreme,temp.min))}
    if (max(unique(s.s.short$temp)) < temp.max) {temp.add.extreme <- sort(c(temp.add.extreme,temp.max))}

    temp.interp <- temp.add.extreme             #helper vector for additional auto-interpolations

    #function to determine temperatures to interpolate by alternating adding a bi-section at the last two cold-hot temperature ranges
    add.temp <- function(cold.hot.indicator,temp.add) {
        if ((cold.hot.indicator %% 2) == 0){
            temp.out <- (temp.add[length(temp.add)] + temp.add[length(temp.add)-1])/2
        }
        if ((cold.hot.indicator %% 2) == 1){
            temp.out <- (temp.add[1] + temp.add[2])/2
        }
        return(temp.out)
    }

    # if number of temperatures is less than 6, including requested extreme temperatures then do addtional interpolation
    if (length(sort(unique(s.s.short$temp))) + length(temp.add.extreme) < 6) {
        for (i in 1:(6 - length(sort(unique(s.s.short$temp)))- length(temp.add.extreme))) { #for each missing temperature to 6
        temp.interp <- c(temp.interp,add.temp(i,temp.add))
        temp.add <- sort(c(temp.add,add.temp(i,temp.add)))
        }
    }


}

temp.interp <- sort(unique(c(temp.interp, temp.interp.user[!(temp.interp.user %in% unique(s.s.short$temp))])))     #combine automatically created temperatures with additonal user requests



##### ADD moduli to df.eng.dens
#if (!(all(diff(unique(lin.el.extract(s.s.short$strain/100, s.s.short$stress, s.s.short$temp)$e.el)<0)))) {stop("Not all E-Moduli are decreasing over temperature, check input stress strain curve and maybe adjust min. strain to use")}

e.el.over.temp.eng.fun <- approxfun(fit.e.el.from.tan.mod(s.s.short,strain.in.percent=FALSE)$e.el$temp,fit.e.el.from.tan.mod(s.s.short,strain.in.percent=TRUE)$e.el$e.el)


s.s.short.temp <- s.s.short
temp.int <- temp.interp
for (temp.int in temp.interp) {

    #interpolate
    if (temp.int < max(unique(s.s.short$temp)) && temp.int > min(unique(s.s.short$temp))) { e.int <- e.el.over.temp.eng.fun(temp.int)}
   
     #extrapolate
    if (temp.int < min(unique(s.s.short$temp)))  {  # to lower temp
        slope <- (e.el.over.temp.eng.fun(sort(unique(s.s.short$temp))[2]) - e.el.over.temp.eng.fun(sort(unique(s.s.short$temp))[1])) / (sort(unique(s.s.short$temp))[2] - sort(unique(s.s.short$temp))[1])
        delta.T <- temp.int - sort(unique(s.s.short$temp))[1]
        e.int <- slope * delta.T + e.el.over.temp.eng.fun(sort(unique(s.s.short$temp)))[1]
       
    }

    if (temp.int > max(unique(s.s.short$temp)))  {  # to higher temp
        slope <- (e.el.over.temp.eng.fun(rev(sort(unique(s.s.short$temp)))[2]) - e.el.over.temp.eng.fun(rev(sort(unique(s.s.short$temp)))[1])) / (rev(sort(unique(s.s.short$temp)))[2] - rev(sort(unique(s.s.short$temp)))[1])
        delta.T <- temp.int - rev(sort(unique(s.s.short$temp)))[1]
        e.int <- slope * delta.T + e.el.over.temp.eng.fun(rev(sort(unique(s.s.short$temp))))[1]
        if (e.int < 0) {e.int <- 1}
    }
   
   
    temp.base <- unique(s.s.short$temp)[which.min(abs(unique(s.s.short$temp) - temp.int))]
    e.base <- e.el.over.temp.eng.fun(temp.base)

    df.int <- s.s.short[s.s.short$temp == temp.base,]
    df.int$temp <- temp.int
   
    if (temp.base < temp.int) {  # softening
        ratio <- e.int/e.base   # softer/harder <1
        df.int$strain <- df.int$strain / ratio
        df.int$stress <- df.int$stress * ratio
    }
   
    if (temp.base > temp.int) {  #hardening
        ratio <- e.int/e.base     #harder/softer >1
        df.int$strain <- df.int$strain / ratio
        df.int$stress <- df.int$stress * ratio
    }
    s.s.short.temp <- rbind(s.s.short.temp,df.int)[order(rbind(s.s.short.temp,df.int)$temp),]
   
   
 
   
}


s.s.short <- s.s.short.temp

if (evaluate.creep) {
    ## s.s. isocr extrapolate
    e.el.over.temp.eng.fun <- approxfun(fit.e.el.from.tan.mod(s.s.short,strain.in.percent=FALSE)$e.el$temp,fit.e.el.from.tan.mod(s.s.short,strain.in.percent=TRUE)$e.el$e.el)  # update approx fun with new created temps


    s.s.ischr.extrap <- s.s.isochr[s.s.isochr$temp == max(unique(s.s.isochr$temp)),]
    s.s.ischr.extrap$strain <-  s.s.ischr.extrap$strain / (e.el.over.temp.eng.fun(s.s.iscrr.temp.to.extrapl) / e.el.over.temp.eng.fun(max(unique(s.s.isochr$temp))))
    s.s.ischr.extrap$temp <- s.s.iscrr.temp.to.extrapl

    s.s.isochr <- rbind(s.s.isochr,s.s.ischr.extrap)
}

## Tabular Matrix Properties
if (evaluate.creep == TRUE){func.s.s.isochr.input.treatment(s.s.isochr)}
s.s.short.true <- func.s.s.short.input.treatment(s.s.short)     #gives correct colnames "temp", "strain", "stress", convert strain in mm/mm, convert to true values

## Plot helpers
color.scale.temp <- rev(hcl(h = seq(1,230,length.out=length(unique(s.s.short$temp))), c = 100, l = 75))  # create temperature color scale for ggplots


##_________ Allowable Stress Settings
##___ Allowable Stress - K - Strength
K.fac <- 1         # Factorized reduction factor K


##___ Allowable Stress - A.stiff - Stiffness Reduction Factors
##_ A.stiff.F: Fiber Content / "Fertigungseinfluss"
A.stiff.F <- 1/1           # Reduction for homogenization of fiber filled material [1/A.stiff.F = 0.6 - 0.7 (1.428571)]
if (fiber.w.perc > 0){A.stiff.F <- 1/0.65}

##___ Allowable Stress - A.str - Strenght Reduction Factors
A.BN <- 1           # Reduction for welding line (2 - 1 for no welding line)
A.ex <- 1          # Reduction for uncertainty of experimental scalar.res
A.stiff <- A.stiff.F

{
##_________ BISO Settings
tan.mod.pivot <- 3       # 3 would mean that tangent modulus for biso is built by last stress-strain point and the 3rd last

# ANSYS allows max 6 temperatures to be used in BISO
# make a equal separation according 5 temperatures, add room temperature and then get rid of duplicates to ensure that room temp is inside
temp.index.to.use <-sort(unique(c(as.integer(seq(1,length(unique(s.s.short$temp)),length.out=5)),match(23,unique(s.s.short$temp)))))


indices <- as.integer(seq(1,length(unique(s.s.short$temp)),length.out=5))
rt. <- match(23,unique(s.s.short$temp))

# ensure that the RT was not matching the initial sequence of temperatures otherwise increase the separation to 6 and remove rt  from usual indices
if (!is.na(match(rt.,indices))){
    indices <- as.integer(seq(1,length(unique(s.s.short$temp)),length.out=6))
    item.remove <- which.min(abs(indices-rt.))
    indices <- indices[!indices==item.remove]
    temp.index.to.use <- sort(unique(c(indices,rt.)))
}

temp.index.to.use <- temp.index.to.use[-7]

###____________Viscoelastic Options
k_multiplicator <- 30               # Modify slope, NORMALLY 10 to reach normal values (hoeher desto flacher)  
E0.viscel.offset <- -2000                 # Global E-Modulus offset
ref_time <- 0.5                   # Usually about 0.5 s, since at this time 0.5 % as limit for linear viscoel is reached
temp_interval <- c(0,140)             # Intervall of termperatures which are used for fit (interval should not cross Tg)
use_extrapolated <- TRUE              # If True, than the Line viscoelastic props are extrapolated to very long times
                # Drawback -> Fitting lambdas have to be modified manually

###____________Time Hardening Options
time.hard.temp.range <- c(0,120)  # Deactivated so far !!!


ref_temp <- 23
create.perzynda.viscoplasticity <- FALSE}



###________________________________________________________________________  Energy density under ref. Curve and Tg approx

s.s.short.true.A.stiff <- s.s.short.true
s.s.short.true.A.stiff$stress <- s.s.short.true.A.stiff$stress * 1/A.stiff.F          # apply the global reduction factor on stress column of short term data
color.scale.temp <- rev(hcl(h = seq(1,230,length.out=length(unique(s.s.short.true.A.stiff$temp))), c = 100, l = 75))  # create temperature color scale for ggplots



df.eng.dens <- sigma.max.UTS.over.temp(s.s.short.true)
colnames(df.eng.dens)[2] <- "strength.datasheet"                                #function to get UTS stress over Temperature [temp,UTS]
df.eng.dens <- cbind(df.eng.dens,                                               #Add brittle/ductile property and break energies
               energy.density.over.temp(s.s.short.true,
               brittle.smaller.than.ratio=brittle.smaller.than.ratio)[,-1])
df.eng.dens$strength.A.stiff <- df.eng.dens$strength.datasheet * 1/A.stiff.F             





room.temp <- get.room.temp(s.s.short.true)                                      # Evaluate closest available temperature to 23degC
neuber.plot.energies(s.s.short.true,room.temp,str.mat.name)                     # Plot Stress-Strain Curve with neuber values from min-rt-max

ggsave(paste(str.mat.name,"/graphs/stress-strain.png",sep=""),plot=gg.stress.strain, width =40, height = 30, units="mm",scale=4)

###____________Assembling Factors
# add sefety.factor column
A.w.s <- A.w.s_factor_from_list(fiber.w.perc,polymer.type)$max                  # retrieve list of load case driven stress reduction factors
df.A.w.s <- cbind(temp=df.eng.dens$temp,A.w.s)                                  # expand A.w.s list to data.frame with row number of temperatures


df.eng.dens$safety.fac <- safety.factors(probability=probability,               # get safety factors according breaktype, probability and severety
                          severety=severety,breaktype=df.eng.dens$breaktype     # , anisofactor currently fixed to 1 according #korte2019
                          ,anisofactor=1)


df.reductionf.temp.load <- (df.A.w.s *(1/A.stiff.F))#(1/df.eng.dens$safety.fac) # Assembling factors to large df

df.ep.allow <- stress.strain.limits.over.temp(                                  # apply reductionfactors on UTS
                df.temp.sig.uts=df.eng.dens$strength.datasheet,
                df.s.s.short=s.s.short.true.A.stiff,
                df.factors=df.reductionf.temp.load)$df.ep.allow
df.sig.allow <- stress.strain.limits.over.temp(
                df.temp.sig.uts=df.eng.dens$strength.datasheet,
                df.s.s.short=s.s.short.true.A.stiff,
                df.factors=df.reductionf.temp.load)$df.sig.allow

 
#df.allowable <- df.sig.allow

df.ep.allow$fit <- fit.function.limits(df.ep.allow,degree=3,xname='BFE',str.nmbr.digits=20)$fit.results
df.sig.allow$fit <- fit.function.limits(df.sig.allow,degree=3,xname='BFE',str.nmbr.digits=20)$fit.results
colnames(df.ep.allow) <- colnames(df.sig.allow) <- c("temp","One Time", "Multiple Times", "Relaxation", "Creep", "Cyclic 10^7","Fit of 'One Time'")


sig.fit.math.string <- fit.function.limits(df.sig.allow,degree=3,xname='BFE',str.nmbr.digits=20)$math.string
sig.fit.math.string.base <- paste("(",fit.function.limits(df.sig.allow,degree=3,xname='BFE',str.nmbr.digits=20)$math.string , ")/",A.w.s$"sig.one.time")

ep.fit.math.string <- fit.function.limits(df.ep.allow,degree=3,xname='BFE',str.nmbr.digits=20)$math.string

sparab.string <- parabolic.crit.fun.ehrenstein(R1=1,R2=1,R3=1,m=m.ratio, scale=1, plots=FALSE)$math.string

workbench.SF <- paste("(",sig.fit.math.string , ")/(", sparab.string ,")",sep="")
workbench.index.datasheet <- paste("(",sparab.string , ")/(", sig.fit.math.string.base ,")",sep="")


A.str <- A.BN * A.ex * 1
A_S_relax <- 1
A_w <- A.s <- 0.17

gg.epallow.temp <- plot.limit.over.temp(df.ep.allow,ylab='Allowable Strain [mm/mm]')
gg.sigallow.temp <- plot.limit.over.temp(df.sig.allow,ylab='Allowable Stress [MPa]')


####__________________________________________________________________ Calculate E-Modulus over temperature and T.g

#df.eng.dens$e.el.a.stiff <- unique(lin.el.extract(s.s.short.true.A.stiff$strain, s.s.short.true.A.stiff$stress, s.s.short.true.A.stiff$temp)$e.el) #obsolete

df.eng.dens$e.el.a.stiff <- fit.e.el.from.tan.mod(s.s.short.true.A.stiff,strain.in.percent=FALSE)$e.el$e.el


#########
#df.eng.dens$e.el <- unique(lin.el.extract(s.s.short.true$strain, s.s.short.true$stress, s.s.short.true$temp)$e.el) #obsolete
df.eng.dens$e.el <- fit.e.el.from.tan.mod(s.s.short.true,strain.in.percent=FALSE)$e.el$e.el



#if (is.null(T.g)) {T.g <- Tg.from.eng.dens(s.s.short.true)}                    # If T.g is not provided estimate it from maximum Break Energy
                                                                                ## This method is deactivated since it searches for max break energy which is dependent
                                                                                ## on the quality of the measuring data about the strain at break.
if (is.null(T.g)) {
    df.eng.dens$e.el.slope <- c(diff(df.eng.dens$e.el)/diff(df.eng.dens$temp),c(diff(df.eng.dens$e.el)/diff(df.eng.dens$temp))[length(df.eng.dens$temp)-1])
                                                                                    # Calculation of slope of measured elastic modulus e.el        
    T.g <-df.eng.dens[df.eng.dens$e.el.slope == min(df.eng.dens$e.el.slope),]$temp  # T.g expected if the slope is at its minimum (max decrease of modulus per temperature.)
}

temp.d.scalar.res <- data.frame()
## here stopped to get rid of temp.d.scalar.res
temp.d.scalar.res <- data.frame(temp=df.eng.dens$temp,e.el=df.eng.dens$e.el)

####__________________________________________________________________ Calculate Poissions ratio over temperature

df.eng.dens$poisson <- 0.3 + 0.2 * (1- (df.eng.dens$e.el/max(df.eng.dens$e.el)))

####__________________________________________________________________ Calculate Mooney Rivlin

df.eng.dens$Sh.A <- (100*df.eng.dens$e.el -15.75) / (2.15+ df.eng.dens$e.el)
df.eng.dens$Mooney.R.2P.C1 <- df.eng.dens$e.el/(6*1.25)
df.eng.dens$Mooney.R.2P.C2 <- df.eng.dens$Mooney.R.2P.C1 * 0.25
df.eng.dens$Mooney.R.2P.D1 <- (2*(1-2*df.eng.dens$poisson))/(df.eng.dens$Mooney.R.2P.C1+df.eng.dens$Mooney.R.2P.C2)
#ggplot(data=df.eng.dens) + geom_line(aes(x=temp,y=Mooney.R.2P.C1,color='C1'))+ geom_line(aes(x=temp,y=Mooney.R.2P.C2,color='C2'))


### __________________________________________________________________ Exract Matrix Modulus



if (fiber.w.perc == 0 ){
    rho.m.rt <- rho
    df.eng.dens$cte.m <- NA
    df.eng.dens$e.m <- NA
    df.eng.dens$cte.m <- NA
    cte.m.rt <- mean(clte.perp.flow.dir,clte.flow.dir)
 
    df.eng.dens$e.m <- df.eng.dens$e.el
    e.m.tg <- df.eng.dens[df.eng.dens$temp <= T.g,]$e.m[length(df.eng.dens[df.eng.dens$temp <= T.g,]$e.m)]
    e.m.rt <-      approxfun(df.eng.dens$temp,df.eng.dens$e.m)(23)
 
    df.eng.dens$specvol.m.scale  <- (( df.eng.dens$e.m/e.m.tg) * ((T.g+273)/(df.eng.dens$temp+273))  -1)
       
            df.eng.dens$cte.m <- cte.m.rt - cte.m.rt * df.eng.dens$specvol.m.scale  #* df.eng.dens$cte.m.scale
            df.eng.dens[df.eng.dens$temp <= T.g, ]$cte.m <- cte.m.rt * e.m.rt/df.eng.dens[df.eng.dens$temp <= T.g, ]$e.m
            df.eng.dens[df.eng.dens$temp > 20 & df.eng.dens$temp < 24,]$cte.m <- cte.m.rt

            df.eng.dens$rho.m <- rho.m.rt/ (1+((df.eng.dens$temp-23)* lin.exp.to.vol.exp(df.eng.dens$cte.m)))
            
        gg.cte.temp <<- ggplot(data=df.eng.dens)+ theme_bw()  +  geom_vline(xintercept=23,linetype=5)+ geom_hline(yintercept=0)+ geom_vline(xintercept=0)+
        geom_line(aes(x=temp,y=cte.m,color="CTE"))+
        geom_point(data=data.frame(temp=c(23),cte=c(clte.flow.dir)),aes(x=23,y=cte,shape="data-sheet Flow Dir."),size=4)+

        scale_y_continuous(limits=c(0,max(df.eng.dens$cte.m)*1.1),name="CTE [10*-6/K]")+
        scale_x_continuous(name="Temperature [degC]")

        ggsave(paste(str.mat.name,"/graphs/cte-temp.png",sep=""),plot=gg.cte.temp, width =40, height = 30, units="mm",scale=4)
}


if (fiber.w.perc > 0 ){

    clte.datasheet.ratio <- clte.flow.dir / clte.perp.flow.dir

    #create interation loop for rho.m.rt and phi.f
    rho.m.rt <- 100
    i <- 1
    while (i <= 100){
        i <- i+1
        phi.f <- Phi(fiber.w.perc/100,rho.f=rho.f,rho.m=rho.m.rt)
        rho.m.rt    <- (rho - rho.f * phi.f) / (1-phi.f)
    }

    df.eng.dens$a.11 <- fiber.orientation
    df.eng.dens$cte.m <- NA
    df.eng.dens$e.m <- NA                                                       # initiating cte.m
    df.eng.dens$rho.m <- rho.m.rt
    df.eng.dens[df.eng.dens$temp > 20 & df.eng.dens$temp < 24,]$cte.m <-        # approximate cte.m.rt from data sheet perp. values
    cte.m.rt <-(clte.perp.flow.dir - phi.f * cte.f) / (1-phi.f)
    phi.f <- 1/(1+((1-(fiber.w.perc/100))/(fiber.w.perc/100))*(rho.f/rho.m.rt))
    phi.m <- 1-phi.f
    k.m <- (k.RT - phi.f * k.f) / (1-phi.f)
    c.m <- (c.RT - phi.f * c.f) / (1-phi.f)


    #df.eng.dens$cte.m <- cte.m.rt <-((clte.perp.flow.dir - phi.f * cte.f) / (1-phi.f) )    #obsolete method to derive cte.m at room temp from volume fraction

    #j <- -1
    #while (j <= 2){
    #j <- j+1

        ## Fitting Orientation degree from CLTE.n from datasheet and Datasheet compound modulus
        #df.eng.dens$a.11 <- E.m.and.orientation.from.composite(df.eng.dens[df.eng.dens$temp > 20 & df.eng.dens$temp < 24,"e.el"][1],
        #                                        clte.flow.dir,
        #                                        (fiber.w.perc+0)/100,cte.f=cte.f,E.f=E.f,nu.f=nu.f,rho.f=rho.f,
        #                                        ar=ar,E.m.start=1000,cte.m=df.eng.dens$cte.m[1],max.iter=50,rho.m=df.eng.dens$rho.m[1],
        #                                        nu.m=df.eng.dens$poisson[1],micromechmodel='tandon.weng')$a.11
        #                                        print(df.eng.dens$a.11)

    #       a.11.temp <- a.11.inc <- df.eng.dens$a.11
            #if fittet a.11 unreasonable use generic
    #        if (!(df.eng.dens$a.11 > 0.7 & df.eng.dens$a.11 < 0.9)[1]) {df.eng.dens$a.11 <- a.11.temp <- a.11.inc <- 0.8}

                for (i in 1:length(df.eng.dens$temp)){

                    ## e.m: Calculate Matrix modulus over temp from orientation                                                               
                    df.eng.dens$e.m[i] <- E.m.from.composite(df.eng.dens$e.el[i],
                                                                psi=fiber.w.perc/100,E.m.limits=c(100,10000),nu.m=df.eng.dens$poisson[i],
                                                                rho.m=df.eng.dens$rho.m[i],E.f=E.f,nu.f=nu.f,rho.f=rho.f,ar=ar,
                                                                a.11=df.eng.dens$a.11[1],
                                                                a.22=(1-df.eng.dens$a.11[1])*0.5,micromechmodel='tandon.weng')
                }


            ##cte.m: Approaches to extract Matrix CLTE over temperature
            e.m.tg <- df.eng.dens[df.eng.dens$temp <= T.g,]$e.m[length(df.eng.dens[df.eng.dens$temp <= T.g,]$e.m)]
            e.m.rt <- df.eng.dens[df.eng.dens$temp > 20 & df.eng.dens$temp < 24,]$e.m
       
            schapery.cte.m <- schapery.cte.m.from.comp.mixing(psi=fiber.w.perc/100, E.m=e.m.tg, E.f=E.f, nu.m=df.eng.dens[df.eng.dens$temp==23,]$poisson,
                                    nu.f=nu.f, cte.n=clte.flow.dir, cte.p=clte.perp.flow.dir, rho.m= df.eng.dens[df.eng.dens$temp==23,]$rho.m ,rho.f=rho.f)
            df.eng.dens$cte.m <- cte.m.rt <- schapery.cte.m$cte.m.p
            clte.flow.dir.corr <-  schapery.cte.m$cte.n.corr
            # TTSP vertical shift factor on cte.m > T.g, Normal e/e(t) factor for cte.m <= T.g
            df.eng.dens$specvol.m.scale  <- (( df.eng.dens$e.m/e.m.tg) * ((T.g+273)/(df.eng.dens$temp+273))  -1)
       
            df.eng.dens$cte.m <- cte.m.rt - cte.m.rt * df.eng.dens$specvol.m.scale  #* df.eng.dens$cte.m.scale
            df.eng.dens[df.eng.dens$temp <= T.g, ]$cte.m <- cte.m.rt * e.m.rt/df.eng.dens[df.eng.dens$temp <= T.g, ]$e.m

            df.eng.dens$rho.m <- rho.m.rt/ (1+((df.eng.dens$temp-23)* lin.exp.to.vol.exp(df.eng.dens$cte.m)))



    ## Create plots over orientation degree
    plot.df.orient <- data.frame()
    j <- 1
    i <- 1
    a.11.start <- 0.33334
    for (i in 1:length(df.eng.dens$temp)){
        nu.m.temp <- df.eng.dens$poisson[i]

        #ar <- 100
        #nu.m.temp <- 0.1
        rho.m.temp <- rho.m.rt / (1+((df.eng.dens$temp[i]-23)* lin.exp.to.vol.exp(df.eng.dens$cte.m[i])))
        df.tandon.weng <- tandon.weng(psi=fiber.w.perc/100,df.eng.dens$e.m[i],E.f=E.f,nu.m=nu.m.temp,nu.f=nu.f,rho.m=rho.m.temp,rho.f=rho.f,ar=ar)
        df.temp <- data.frame(a.11=sort(c(seq(a.11.start,0.999,length.out=100),0.8)))
        df.temp$temp <- df.eng.dens$temp[i]
        df.temp$e.m <- df.eng.dens$e.m[i]
        df.temp$nu.m <- nu.m.temp
        df.temp$e.el <- df.eng.dens$e.el[i]
        df.temp$cte.m <- df.eng.dens$cte.m[i]
        df.temp$E.n <- NA
        df.temp$E.p <- NA
        df.temp$clte.n <- NA
        df.temp$clte.p <- NA


        for (j in 1:length(df.temp$a.11)){
   
            df.temp$E.n[j] <- advani.tucker.tensile.homogenization(a.11=df.temp$a.11[j],a.22=(1-df.temp$a.11[j])*0.5,df.tandon.weng$E.11,df.tandon.weng$E.22,df.tandon.weng$G.12,df.tandon.weng$G.23,df.tandon.weng$nu.12,df.tandon.weng$nu.23,df.tandon.weng$nu.21)$E.ho.n
            df.temp$E.p[j] <- advani.tucker.tensile.homogenization(a.11=df.temp$a.11[j],a.22=(1-df.temp$a.11[j])*0.5,df.tandon.weng$E.11,df.tandon.weng$E.22,df.tandon.weng$G.12,df.tandon.weng$G.23,df.tandon.weng$nu.12,df.tandon.weng$nu.23,df.tandon.weng$nu.21)$E.ho.p
   
            cte.comp <- lee.comp.cte(psi=fiber.w.perc/100,
                        E.m=df.eng.dens$e.m[i],
                        E.f=E.f,
                        nu.m=df.eng.dens$poisson[i],
                        nu.f=nu.f,
                        cte.m=df.eng.dens$cte.m[i],
                        rho.m=df.eng.dens$rho.m[i],
                        rho.f=rho.f,
                        a.11=df.temp$a.11[j])
 
                   
        df.temp$clte.n[j] <- cte.comp$cte.n
        df.temp$clte.p[j] <- cte.comp$cte.p
        df.temp$a.p.fac[j] <- cte.comp$a.p.fac
        df.temp$a.n.fac[j] <- cte.comp$a.n.fac

            }

        plot.df.orient <- rbind(plot.df.orient,df.temp)
    }



    df.eng.dens$E.p <- plot.df.orient[plot.df.orient$a.11==df.eng.dens$a.11[1],]$E.p
    df.eng.dens$E.iso <- plot.df.orient[plot.df.orient$a.11==a.11.start,]$E.p
    df.eng.dens$inv.A.stiff <-  df.eng.dens$E.iso /df.eng.dens$e.el
    df.eng.dens$cte.iso <- plot.df.orient[plot.df.orient$a.11==min(plot.df.orient$a.11),]$clte.n
    df.eng.dens$cte.n <- plot.df.orient[plot.df.orient$a.11==max(plot.df.orient$a.11),]$clte.n
    df.eng.dens$cte.p <- plot.df.orient[plot.df.orient$a.11==max(plot.df.orient$a.11),]$clte.p

    ggplot(data=plot.df.orient)+ theme_bw() + geom_vline(xintercept=0.33)+ geom_hline(yintercept=0)+
        geom_line(aes(x=a.11,y=clte.n*10^6,color=as.factor(temp),linetype="parallel"))+
        geom_line(aes(x=a.11,y=clte.p*10^6,color=as.factor(temp),linetype="normal"))+
        geom_point(data=df.eng.dens, aes(x=rep(0.33,length(df.eng.dens$temp)),y=cte.m*10^6,color=as.factor(temp),shape="calc. Matrix CLTE"),size=4)+
        scale_y_continuous(limits=c(0,max(plot.df.orient$clte.p*10^6)*1.1),name="CLTE [10^-6/K]")+
        scale_x_continuous(limits=c(0.33,1),name="Max. Orientation Eigenvalue a.11 [-]")+
        geom_point(data=data.frame(temp=c(23,23),cte=c(clte.perp.flow.dir,clte.flow.dir)),aes(x=1,y=cte*1000000,shape="data-sheet",color=as.factor(temp)),size=4)+
        geom_point(data=data.frame(temp=c(23),cte=c(clte.flow.dir.corr)),aes(x=1,y=cte*1000000,shape="Schapery Correction",color=as.factor(temp)),size=4)+


        scale_colour_manual(name='temp [degC]',values = color.scale.temp)+
        scale_linetype_manual(name='Load Direction',values = c(1,5))


    ##example tandon.weng plot over fiber weight ratio

        gg.modulus.orientation <<- ggplot(data=plot.df.orient)+ theme_bw() + geom_vline(xintercept=0.33)+ geom_hline(yintercept=0)+
        geom_line(aes(x=a.11,y=E.p,color=as.factor(temp),linetype="parallel"))+
        geom_line(aes(x=a.11,y=E.n,color=as.factor(temp),linetype="normal"))+
        geom_point(aes(x=df.eng.dens$a.11[1],y=e.el,color=as.factor(temp)))+
        scale_y_continuous(limits=c(0,max(plot.df.orient$E.n)*1.1),name="Modulus [MPa]")+
        scale_x_continuous(limits=c(0.33,1),name="Max. Orientation Eigenvalue a.11 [-]")+
        scale_colour_manual(name='temp [degC]',values = color.scale.temp)+
        scale_linetype_manual(name='Load Direction',values = c(1,5))
   
   
        gg.inv.A.stiff.temp <<- ggplot(data=df.eng.dens)+ theme_bw() +  geom_vline(xintercept=23,linetype=5)+ geom_hline(yintercept=0)+
        geom_line(aes(x=temp,y=inv.A.stiff))+
        scale_y_continuous(limits=c(0.4,.75),name="1/A.stiff")+
        geom_hline(yintercept=0.65)+
        scale_x_continuous(name="Temperature [degC]")

   
        gg.e.m.temp <<- ggplot(data=df.eng.dens)+ theme_bw() +  geom_vline(xintercept=23,linetype=5)+ geom_hline(yintercept=0)+
        geom_line(aes(x=temp,y=e.m,color="Matrix Modulus"))+
        geom_line(aes(x=temp,y=E.p,,color="Tensile Bar perpendicular flow Modulus"))+
        geom_line(aes(x=temp,y=e.el,,color="Tensile Bar flow direction Modulus"))+
        geom_line(aes(x=temp,y=E.iso,,color="Quasi-Isotropic Modulus"))+
        scale_y_continuous(limits=c(0,max(df.eng.dens$e.el)*1.1),name="Elastic Modulus [MPa]")+
        scale_x_continuous(name="Temperature [degC]")
   
        gg.cte.temp <<- ggplot(data=df.eng.dens)+ theme_bw()  +  geom_vline(xintercept=23,linetype=5)+ geom_hline(yintercept=0)+ geom_vline(xintercept=0)+
        geom_line(aes(x=temp,y=cte.m*1000000,color="Matrix CTE"))+
        geom_line(aes(x=temp,y=cte.p*1000000,color="Perp. Flow CTE"))+
        geom_line(aes(x=temp,y=cte.n*1000000,color="In Flow CTE"))+
        geom_line(aes(x=temp,y=cte.iso*1000000,color="Isotropic CTE"))+
        geom_point(data=data.frame(temp=c(23),cte=c(clte.flow.dir)),aes(x=23,y=cte*1000000,shape="data-sheet Flow Dir."),size=4)+
        geom_point(data=data.frame(temp=c(23),cte=c(clte.perp.flow.dir)),aes(x=23,y=cte*1000000,shape="data-sheet Perp. Flow Dir."),size=4)+
        geom_point(data=data.frame(temp=c(23),cte=c(clte.flow.dir.corr)),aes(x=23,y=cte*1000000,shape="Perp. Flow Dir Schapery Correction"),size=4)+

        scale_y_continuous(limits=c(0,max(df.eng.dens$cte.p*1000000)*1.1),name="CTE [10*-6/K]")+
        scale_x_continuous(name="Temperature [degC]")

        grid.arrange(    gg.e.m.temp,gg.inv.A.stiff.temp,gg.modulus.orientation,layout_matrix = rbind(c(1,1),c(2,3)))
     

        ggsave(paste(str.mat.name,"/graphs/aniso-stifness.png",sep=""),plot=gg.modulus.orientation, width =50, height = 40, units="mm",scale=3)
        ggsave(paste(str.mat.name,"/graphs/cte-temp.png",sep=""),plot=gg.cte.temp, width =40, height = 30, units="mm",scale=4)

        ### Export scaled plasticity Stress-Strain cuves per temperature for 0 and 90 degree orientation
        
        df.s.s.scale <- plot.df.orient[plot.df.orient$a.11 == 0.999,]  #create df and select only highest orientations of given plot.df for orientations 
        df.s.s.scale$scale.n <- df.s.s.scale$E.n / df.s.s.scale$e.el   #calculate scale factor by Modulus in normal directon and datasheet modulus
        df.s.s.scale$scale.p <- df.s.s.scale$E.p / df.s.s.scale$e.el   #calculate scale factor by Modulus in perpendicular directon and datasheet modulus
        

        str.dir.pl.orient.folder.name <- paste(str.mat.name,"/Fiber-Files/","point_data-plStrain_Orient",sep="")
        dir.create(str.dir.pl.orient.folder.name)
        
        
}


## ___________________________________________________________________ Calculate Tangent Modulus for bilinear Hardening (BISO)
## BISO is calculated from un-homogenized data.sheet values --> only to be used in converse material output
df.eng.dens$tan.mod.biso <- NaN                   # craeate new colum
df.eng.dens$yield <- NaN
df.eng.dens$m.yield <- NaN

matrix.composite.ratio <- df.eng.dens$e.m/df.eng.dens$e.el



for (i in seq(1,length(unique(df.eng.dens$temp)))){               # for all existing temepratures in scalar res
  df.i <- s.s.short.true[s.s.short.true$temp == unique(s.s.short.true$temp)[i],]      # create filter for this temperature of s.s curves
  y2 <- tail(df.i,2)$stress[2]
  y1 <- tail(df.i,tan.mod.pivot)$stress[1]
  x2 <- tail(df.i,2)$strain[2]
  x1 <- tail(df.i,tan.mod.pivot)$strain[1]
  df.eng.dens[df.eng.dens$temp == unique(s.s.short.true$temp)[i],]$tan.mod.biso <- E.t <- (  # calculate tangent modulus from two points based on tan.mod.pivot
  (y2 - y1) / (x2 - x1))
  t0.biso <- y2-(E.t*x2)
  e.el.biso <- df.eng.dens$e.el[i]

  df.eng.dens[df.eng.dens$temp == unique(s.s.short.true$temp)[i],]$yield <- (t0.biso)/(df.eng.dens$e.el[i]-E.t) * df.eng.dens$e.el[i]
  df.eng.dens[df.eng.dens$temp == unique(s.s.short.true$temp)[i],]$m.yield <- ((t0.biso)/(df.eng.dens$e.el[i]-E.t) * df.eng.dens$e.el[i]) * matrix.composite.ratio[i]
  }
 
rm(df.i, x1,x2,y1,y2,t0.biso)



##___________________________________________________________________ Critical Strain reduction by fiber content
temp.d.scalar.res$crit.strain <- critical.strain.menges(polymer.type, fiber.w.perc,unique(s.s.short.true.A.stiff$temp),T.g)$crit.strain
df.ep.allow[["Critical Strain"]]    <- temp.d.scalar.res[["crit.strain"]]
df.sig.allow[["Critical Stress"]]   <- fun.approx.s.s.hort(strain=df.ep.allow[["Critical Strain"]],temp=df.ep.allow$temp)

####__________________________________________________________________ Calculate Poissions ratio over temperature #### obsolete for temp.d.scalar.res
temp.d.scalar.res$poisson <- 0.3 + 0.2 * (1- (temp.d.scalar.res$e.el/max(temp.d.scalar.res$e.el)))

#ggplÃ¶t
par(mfrow = c(1,1))
plot(temp.d.scalar.res$temp, temp.d.scalar.res$poisson
   ,main = 'Apprx. Poisson Ratio over Temperature'
   ,type = 'o', xlab = 'temperature [degC]', ylab = 'Poisson Ratio [-]'); grid()
par(mfrow = c(1,1))

### ____________________________________________________________________________ Modulus-Strenght correlation

df.eng.dens$strength.A.stiff


#grid plot of Temp. Modulus Strength correlation with linear regression
main.title <- "Modulus-Strength Corellation by Temperature"
gg.common <- theme_bw() + theme( axis.line = element_line( size = 1, linetype = "solid"))

gg.e.el.vs.temp <- ggplot(df.eng.dens,aes(x=temp, y=e.el)) +geom_line()+ labs(x='temperature[degC]',y='Elastic Modulus [MPa]') + gg.common
gg.strength.vs.temp <- ggplot(df.eng.dens,aes(y=temp, x=strength.A.stiff)) +geom_line() + labs(y='temperature[degC]',x='Strength [MPa]') + gg.common
gg.e.el.vs.strength <- ggplot(df.eng.dens,aes(x=strength.A.stiff, y=e.el))+geom_line() + labs(x='Strength [MPa]',y='Elastic Modulus [MPa]') + gg.common+
geom_smooth(method='lm',se = FALSE,size=5,linetype=3) +
geom_text(x = mean(df.eng.dens$strength.A.stiff), y = mean(df.eng.dens$e.el), label = lm_eqn(df.eng.dens$strength.A.stiff,df.eng.dens$e.el), parse = TRUE)


grid.arrange(gg.e.el.vs.temp,gg.e.el.vs.strength,gg.strength.vs.temp,fit.e.el.from.tan.mod(s.s.short,strain.in.percent=FALSE)$gg.e.tan,layout_matrix = rbind(c(1,2),c(4,3)),top=textGrob(main.title,gp=gpar(fontsize=20,font=3)))


###____________________________________________________________________ Plasticity short term MISO
s.s.short.true.pl <- lin.el.extract(s.s.short.true.A.stiff$strain, s.s.short.true.A.stiff$stress, s.s.short.true.A.stiff$temp)[, c("temp", "pl.strain", "stress")]
s.s.e.sec <- secant.modulus(s.s.short.true.A.stiff$strain, s.s.short.true.A.stiff$stress, s.s.short.true.A.stiff$temp)


#ggplÃ¶t
par(mfrow = c(2,2))      # Plot
plot.multicurve.xy(s.s.short, c('Stress-Strain: ', str.mat.name), 'strain [mm/mm]', 'stress [MPa]', 'Temp.')
plot.multicurve.xy(s.s.short.true.pl, c('Plastic Strain: ', str.mat.name), 'plastic true strain [mm/mm]', 'true stress [MPa]', 'Temp.')
plot.multicurve.xy(s.s.short.true.A.stiff, c('True Stress-Strain: ', str.mat.name), 'true strain [mm/mm]', 'true stress [MPa]', 'Temp.')
plot.multicurve.xy(s.s.e.sec[,!(colnames(s.s.e.sec)=="stress")], c('Tr. Secant Modulus: ', str.mat.name), 'true strain [mm/mm]', 'tr. Secant Modulus [MPa]', 'Temp.')
par(mfrow = c(1,1))

###____________________________________________________________________ Plasticity short term BISO


###____________________________________________________________________ Approximation of CLTE
temp.d.scalar.res$clte.appr <- 1 / temp.d.scalar.res$e.el
par(mfrow = c(1,1))

#ggplÃ¶t
#plot(temp.d.scalar.res$temp, temp.d.scalar.res$clte.appr * 10^(6), ylab='clte approx [10^-6/K] ', xlab='Temperature[degC]', main='Approximation of isotropic CLTE', type='o')
grid()


###________________________________________________________________________  Long term durability 10e7 with mean stress


#__________________________sig'_WK______________________________________________  Part Limiting Stress Amplitude 10^7 (de: Bauteilwechselfestigkeit)

df.fatique <- data.frame(temp=df.A.w.s$temp , A.w=A.w.s_factor_from_list(fiber.w.perc,polymer.type)$max$sig.cyclic
                        ,F=df.eng.dens$strength.datasheet)  #
df.fatique$sig.wk <- df.fatique$A.w * df.fatique$F / A.stiff                                              # Part Limiting Stress Amplitude 10^7 korte2019p231




#__________________________ sig'_AK = Slope of Goodman approach______________  Part Limiting Stress Amplitude 10^7 + (de: Bauteilgrenzspannungsamplitude)
#                                                                              with mean stress according Goodman Line

Mean.stress.limit <- df.fatique$F / A.stiff                                 # UTS reduced by fiber factor for Goodman line

sig.ak.res          <- sig.ak(temp=df.fatique$temp,sig.wk=df.fatique$sig.wk,
                                Mean.stress.limit=Mean.stress.limit)
sig.wk.temp.fun     <- sig.ak.res$sig.wk.fit$max.limit.fun                  # sig.wk(temp): math function long therm durability 10^7 over temp
df.fatique$M.sig    <- sig.ak.res$M.sig                                     # M.sig:        list Mean Stress sensitivity     
sig.ak.string       <- sig.ak.res$string                                    # sig.ak($F1-$F4): math string mean stress reduced long therm durability 10^7
sig.ak.fun          <- sig.ak.res$fun                                       # sig.ak($F1-$F4)(sig.m,sig.wk,temp): math function mean stress reduced long therm durability 10^7

#__________________________ sig'_BK _________________________________ Part Limiting Stress Amplitude N + (de: Bauteil-Zeitfestigkeit N<10^7)
#                                                                       with mean stress according Goodman Line to UTS

k.slope <-  (log10(10^7))/                                                                                      # slope of sig'bk over N
            (log10(A.w.s_factor_from_list(fiber.w.perc,polymer.type)$max$sig.one.time /
             A.w.s_factor_from_list(fiber.w.perc,polymer.type)$max$sig.cyclic))
sig.sk.temp.fun <- fit.function.limits(df.sig.allow,degree=3,xname='sparab',str.nmbr.digits=10)$max.limit.fun    # fit function of allowable stress over temperature

K.bk.fun <- function(N,k.slope,N.G=1e7) {(N.G/N)^(1/k.slope)}                                                   # "Increase Factor" for Load cycles <1e7

#__ sig'_BK _result assembly
sig.bk.fun <- list( F1=function(sig.m,temp,N)    {sig.ak.fun$F1(sig.m,temp) * K.bk.fun(N,k.slope=k.slope,N.G=1e7)},
                    F2=function(sig.m,sig.a,temp,N)     {sig.ak.fun$F2(sig.m, sig.a, temp) * K.bk.fun(N,k.slope=k.slope,N.G=1e7)},
                    F3=function(sig.u,temp,N)    {sig.ak.fun$F3(sig.u, temp)  * K.bk.fun(N,k.slope=k.slope,N.G=1e7)},
                    F4=function(sig.o,temp,N)    {sig.ak.fun$F4(sig.o, temp)  * K.bk.fun(N,k.slope=k.slope,N.G=1e7)}
)

sig.bk.string <-   list(F1=paste('(',sig.ak.string$F1 ,')*((1000000/N)^(1/',k.slope,'))' ),            
                        F2=paste('(',sig.ak.string$F2 ,')*((1000000/N)^(1/',k.slope,'))' ),
                        F3=paste('(',sig.ak.string$F3 ,')*((1000000/N)^(1/',k.slope,'))' ),
                        F4=paste('(',sig.ak.string$F4 ,')*((1000000/N)^(1/',k.slope,'))' )
)

#__________________________ sig'_BK_gerber_____________________________ Part Limiting Stress Amplitude N + (de: Bauteil-Zeitfestigkeit N<10^7)
#                                                                       with mean stress according Gerber parabolic function to short term limit

sig.sk.temp.string <- fit.function.limits(df.sig.allow,degree=3,xname='BFE',str.nmbr.digits=15)$math.string
sig.wk.temp.string <- fit.function.limits(df.fatique,degree=3,xname='BFE',str.nmbr.digits=15,x.col=1,y.col=4)$math.string
T.string <- paste('(',sig.wk.temp.string,')*(1-(sigm/(',sig.sk.temp.string,'))^2)')         # helpers for ger.string.N.from.sig
i.string <- paste('(1/',k.slope,')')

#__ sig'_BK_gerber_result assembly

sig.bk.gerber.fun <- function(sig.m,temp,N) { sig.wk.temp.fun(temp) * (1-(sig.m/sig.sk.temp.fun(temp))^2) * K.bk.fun(N,k.slope=k.slope,N.G=1e7)}

gerber.string.sig.from.N <- paste('=(',sig.wk.temp.string,')*((10000000/N)^(1/',k.slope,'))' ,'*(1-(sigm/(',sig.sk.temp.string,'))^2)')
#gerber.string.N.from.sig <- paste('10^7*(sparab/(',T.string,')^(-1/',i.string,'))')


gerber.string.N.from.sig <-   paste(
"10000000/(((",sparab.string,") / (",sig.wk.temp.string, " * (1-(sigm/(",sig.sk.temp.string,"))^2)))^(",k.slope,"))")

sigm.string <- paste("index_avg*(",sig.fit.math.string,")",sep=" ")



##_avg to be calculated within ansys user-defined result identifier from (minimum+maximum)/2 over time from this formula
##workbench.index.datasheet <- paste("(",sparab.string , ")/(", sig.fit.math.string ,")",sep="")
## sigm.string is only considered for the strength reduction, but not for the amplitude reduction to avoid negative values, this means a worst case
gerber.string.N.from.sig.sigm <-   paste(
"10000000/(((" , sparab.string , ") / (",sig.wk.temp.string, " * (1-((",sigm.string,")/(",sig.sk.temp.string,"))^2)))^(",k.slope,"))")

gerber.string.N.from.sig <-   paste(
"10000000/(((" , sparab.string , ") / (",sig.wk.temp.string, " * (1-((","0",")/(",sig.sk.temp.string,"))^2)))^(",k.slope,"))")

## Gerber string with reduction factor

gerber.reduction.A <- 2
gerber.string.N.from.sig.reduced <-   paste(
"10000000/(((" , sparab.string , ") / (((",sig.wk.temp.string, ")*1/",gerber.reduction.A,") * (1-((","0",")/(",sig.sk.temp.string,"))^2)))^(",k.slope,"))")



writeClipboard(gerber.string.N.from.sig.reduced)

#
#################
######## FUN for normalized N from sig stresses to include proper mean stress approach
###############

auslastung.string <- paste('(',sparab.string,")/(",sig.sk.temp.string,")",sep="")
## create user-defined restults with "auslastung.string" at time for low temp. Identifier: "lt"
## create user-defined restults with "auslastung.string" at time for high temp. Identifier: "ht"
writeClipboard(auslastung.string)

#mean.auslastung.string <- "abs(lt + ht)/2" #if the mean is to be assumed as average of high and low temp
mean.auslastung.string <- "rt"  # if mean is assumed after room temperature conditons


# --> simple fixed amplitude version
amplitude.string <- "abs(lt - ht)/2"
N.Cycles.normalized.string <- paste("10000000/(((",amplitude.string,") / (",(A.w.s_factor_from_list(fiber.w.perc,polymer.type)$max$sig.cyclic) / A.stiff ," * (1-(",mean.auslastung.string,")^2)))^(",k.slope,"))",sep="")



# --> simple fixed amplitude
amplitude.string <- "abs(lt - ht)/2"
N.Cycles.normalized.string <- paste("10000000/(((",amplitude.string,") / (",(A.w.s_factor_from_list(fiber.w.perc,polymer.type)$max$sig.cyclic) / A.stiff ," * (1-(",mean.auslastung.string,")^2)))^(",k.slope,"))",sep="")

## create user-defined restults with "N.Cycles.normalized.string"

#writeClipboard(N.Cycles)


#__________________________ ep'_BK_gerber______________________________________________________________________ Part Limiting Stress Amplitude N + (de: Bauteil-Zeitfestigkeit N<10^7)
#                                                                                                               with mean stress according Gerber parabolic function to short term limit
# creates df from functions as basis for plotting and strain fitting
plot.fat.res <- plot.fatique.results(temp.d.scalar.res,sig.bk.fun,sig.bk.gerber.fun,color.scale.temp,df.sig.allow,df.ep.allow)  # this also creates the gg grobs for cyclic plotting
df.bk.fat <- plot.fat.res$df.bk.fat

##fitting strain fatique # still free style formula **1 add strain life approach
para <- list(z=0,a=0,b=0,c=2,d=1,e=0)
resiFun <- function(N,para,strain,temp,ep.m) (strain - (    para$a * ((1/N)) ^para$b                                   
                                                * (temp*para$c) ^para$z  - (ep.m*para$d)^2 *(N^para$e )))
fatique.str.fit.df <- df.bk.fat[rowSums(is.na(df.bk.fat[,]))==0,]                                                           #only use rows where ep.m could be calculated

fat.strain.fit <- nls.lm(par=para, fn = resiFun, N = fatique.str.fit.df$N,
                        temp=fatique.str.fit.df$temp+273.15, strain = fatique.str.fit.df$gerber.strain,
                        ep.m=fatique.str.fit.df$ep.m, control = nls.lm.control(nprint=1))

para <- fat.strain.fit$par
fat.strain.fun <- function(ep.m,temp,N)(  para$a * ((1/N)) ^para$b *  (temp*para$c)^para$z - (ep.m*para$d)^2 *(N^para$e )  )
df.bk.fat$gerber.strain.fit <- fat.strain.fun(df.bk.fat$N,df.bk.fat$temp+273.15,df.bk.fat$ep.m)          # calculate function results for plotting

fat.strain.string <- paste(para$a ,"* ((1/N)) ^",para$b," *  (BFE*",para$c,")^",para$z," - (epmean*",para$d,")^2 *(N^",para$e," )")


## add tabular results

df.sig.allow[["Cyclic 10^4"]] <- sig.bk.gerber.fun(0,df.sig.allow$temp,10000)
df.ep.allow[["Cyclic 10^4"]]  <- fun.approx.s.s.hort(stress=df.sig.allow[["Cyclic 10^4"]],temp=df.sig.allow$temp)



### plottings cyclic

if (plot.fatique.curves == TRUE) {
    grid.arrange(gg.mean.stress,gg.strength.fat,nrow=2,ncol=1,top=textGrob(main.title,gp=gpar(fontsize=20,font=3)))
    grid.arrange(gg.mean.strain,gg.strain.fat,nrow=2,ncol=1,top=textGrob(main.title,gp=gpar(fontsize=20,font=3)))
}

ggsave(paste(str.mat.name,"/graphs/fatigue.png",sep=""),plot=gg.strength.fat, width =40, height = 30, units="mm",scale=4)

#grid_arrange_shared_legend(gg.mean.stress,gg.strength.fat, ncol = 1, nrow = 2,position="right")
#grid_arrange_shared_legend(plot.fat.res$gg.mean.strain,plot.fat.res$gg.strain.fat, ncol = 1, nrow = 2,position="right")




###________________________________________________________________________  Viscoelasticity

# Variable declaration
time_dec <- seq(-9,9,1)
time <- c()
modulus <- c()
strain_rate <- (50/60) / 100  # 50 mm/min * 60s / 100% [mm/mm*2]

s.s.e.sec$time <- s.s.e.sec$strain / strain_rate

# Arrhenius comparison-temperature is just minimum temperature in min(s.s.e.sec$temp), time(secantmodulus)
fit_comp_inv <- lm(s.s.e.sec[s.s.e.sec$temp == min(s.s.e.sec$temp),]$time ~ s.s.e.sec[s.s.e.sec$temp == min(s.s.e.sec$temp),]$e.sec
         +I(s.s.e.sec[s.s.e.sec$temp == min(s.s.e.sec$temp),]$e.sec ^2)+ I(s.s.e.sec[s.s.e.sec$temp == min(s.s.e.sec$temp),]$e.sec ^3))
fit_func_comp_inv <- function(x) fit_comp_inv$coefficients[[1]] + x * fit_comp_inv$coefficients[[2]] + x^2 * fit_comp_inv$coefficients[[3]] + x^3 * fit_comp_inv$coefficients[[4]]


##____ calc Arrhenius K factor from the times of reaching the same Modulus for two different temperatures
tref <- s.s.e.sec[s.s.e.sec$temp == ref_temp,]$time[1]        # Time at reference-Temperature to reach first secant modulus value
tcomp <- fit_func_comp_inv(s.s.e.sec[s.s.e.sec$temp == ref_temp,]$e.sec[1])    # Time at compasrison-Temperature to reach first secant modulus value
log.at <-

k <- (log10(tref / tcomp)) / (1/(ref_temp + 273.3) - 1/(min(s.s.e.sec$temp) + 273.3)) * k_multiplicator ##obsolete
#k_multiplicator <- 1

#k <- wlf.equation(T.0=T.g+50,T=23)$k * k_multiplicator


############################################# STopped here beim AufrÃ¤umen von time.stress.e.sec 30.1.17, darum neudeklaration damit es weirhin funktioniert
# time.stress.e.sec (1 temp,2 time,3 e.sec)
time.stress.e.sec <- data.frame( s.s.e.sec$temp, s.s.e.sec$time, s.s.e.sec$e.sec)

##____ Fit every secant modulus curve in the given range and evaluate Modulus at reference time
##---- then calculate shifted time by the k factor
for(temp in unique(s.s.e.sec$temp)){
 if (temp >= min(temp_interval) && temp <= max(temp_interval)){
  fit <- lm(time.stress.e.sec[time.stress.e.sec[1] == temp, 3] ~  time.stress.e.sec[time.stress.e.sec[1] == temp, 2]
       +I(time.stress.e.sec[time.stress.e.sec[1] == temp, 2] ^2)+ I(time.stress.e.sec[time.stress.e.sec[1] == temp, 2] ^3))
  fit_func <- function(x) fit$coefficients[[1]] + x * fit$coefficients[[2]] + x^2 * fit$coefficients[[3]] + x^3 * fit$coefficients[[4]]
 
  modulus <- append(modulus, fit_func(ref_time))
  time <- append(time, ref_time * 10^(k * (1/(ref_temp + 273.3)-1/(temp + 273.3))))
 }}


# Put scalar.res of shifted moduli in one package and fit for overall linear function by fitting for log10(t)
moduli <- data.frame(time, modulus)

fit_shifted <- lm(moduli$modulus ~ log10(moduli$time) +I(log10(moduli$time)^1))
fit_shifted_func <- function(x) fit_shifted$coefficients[[1]] + x * fit_shifted$coefficients[[2]]

# Extrapolate Moduli for all decades by fitted function
time_extra <- 1 * 10^(time_dec)
moduli_extra <- fit_shifted_func(time_dec)

moduli_extra <- data.frame(time_extra, moduli_extra)


if(use_extrapolated == T){modulus <- moduli_extra}
if(use_extrapolated == F){modulus <- moduli}


t <- modulus[[1]]

E <- modulus[[2]]
Estart <- 1/8
E_0 <- max(time.stress.e.sec[[3]]) + E0.viscel.offset

## Apply Prony fit on the extrapolated moduli
lambda <- c(1e-8, 1e-6, 1e-4, 1e-2, 1e0, 1e2, 1e4, 1e6, 1e9, 1e11)
para <- list(E1=Estart, E2=Estart,E3=Estart,E4=Estart,E5=Estart,E6=Estart,E7=Estart,E8=Estart,E9=Estart,E10=Estart)

resiFun <- function(t, para, E) (E - (E_0 * (para$E1 * exp(-1 * t/lambda[1])
                       +para$E2 * exp(-1 * t/lambda[2])
                       +para$E3 * exp(-1 * t/lambda[3])
                       +para$E4 * exp(-1 * t/lambda[4])
                       +para$E5 * exp(-1 * t/lambda[5])
                       +para$E6 * exp(-1 * t/lambda[6])
                       +para$E7 * exp(-1 * t/lambda[7])
                       +para$E8 * exp(-1 * t/lambda[8])
                       +para$E9 * exp(-1 * t/lambda[9])
                       +para$E10 * exp(-1 * t/lambda[10])
)
))

nls.out.visco <- nls.lm(par=para, fn = resiFun, E = E
            , t = t, control = nls.lm.control(nprint=1), lower = c(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
)

E_i <- unlist(nls.out.visco[[1]]) / max(cumsum(unlist(nls.out.visco[[1]]))) # Normalize E_i values for sum 1

fitted <- function(x) (E_i[1] * exp(-1 * x/lambda[1])
            +E_i[2] * exp(-1 * x/lambda[2])
            +E_i[3] * exp(-1 * x/lambda[3])
            +E_i[4] * exp(-1 * x/lambda[4])
            +E_i[5] * exp(-1 * x/lambda[5])
            +E_i[6] * exp(-1 * x/lambda[6])
            +E_i[7] * exp(-1 * x/lambda[7])
            +E_i[8] * exp(-1 * x/lambda[8])
            +E_i[9] * exp(-1 * x/lambda[9])
            +E_i[10] * exp(-1 * x/lambda[10])
      
      
)

par(mfrow = c(1,1))
plot.df <- data.frame(x=t, y=fitted(t)*E_0 )
arrh.points <- data.frame(x=moduli[[1]],y=moduli[[2]])
gg <- ggplot(data=plot.df,aes(x=x,y=y)) +geom_line() + scale_x_log10()
gg <- (gg + ggplot_theme +scale_colour_brewer(palette = "Dark2")
+ labs(title='Youngs Modulus over Time', x='time [s]', y='Modulus [MPa]')
        + geom_point(data=arrh.points,aes(x=x,y=y))
        #+ scale_x_continuous(labels = logplot.df$x,breaks = plot.df$x)
        #+ scale_y_continuous(breaks = seq(-200, 200,20))
        #+ scale_color_brewer(name='versions',palette='Set2')
        #+ facet_grid(variable~., scales='free')
        #+ xlim(-6,6)
        #+ ylim(-6,6)
 
    )
#Plot Points from creep curves at smallest stress
if (evaluate.creep == TRUE){
  plot.modulus <- data.frame(temp=c(),stress=c(), time=c(),strain=c(),modulus=c())
  for (temp in unique(s.s.isochr$temp)) {
    s.s.plot.temp <- s.s.isochr[s.s.isochr$temp == temp,]
    for (i in unique(s.s.plot.temp$time)) {  # filter for minimum stress per time
        plot.modulus <- rbind(plot.modulus,s.s.plot.temp[s.s.plot.temp$time == i,][1,])
    }
 
  }
  gg <- gg + geom_point(data=plot.modulus[plot.modulus$temp == ref_temp,],aes(x=time,y=modulus),color='green')
  gg <- gg + geom_point(data=plot.modulus[plot.modulus$temp == max(plot.modulus$temp),],aes(x=time,y=modulus),color='red')
  gg <- gg + geom_point(data=plot.modulus[plot.modulus$temp == min(plot.modulus$temp),],aes(x=time,y=modulus),color='blue')
}


print(gg)

df_prony <- data.frame(E_i, lambda)

temp.d.scalar.res$E0 <- E_0



if (evaluate.creep){
###____________________________________________________________________________________ Modified Time Hardening Fit
###____________________________________________Only made if creep modulus file is given, enlarged C4 - Higher spread at each temp
    s.s.isochr.org <- s.s.isochr
    TimeHardpara.out <- data.frame(temp=NA,C1=NA,C2=NA,C3=NA,C4=0)
    s.s.isochr.org$fit <- NA
    if (evaluate.creep == TRUE && is.element(ref_temp, unique(s.s.isochr.org$temp))){
    for (creep.temp in unique(s.s.isochr.org$temp)) {
        print(creep.temp)
        s.s.isochr <- s.s.isochr.org[ s.s.isochr.org$temp == creep.temp,]
        #s.s.isochr <- s.s.isochr[ s.s.isochr$temp >= min(time.hard.temp.range),]
    
        para <- list(C1=0,C2=0,C3=-0.5,C4=0)        # Initialize Parameters of Modified Time Hardening
        C4.lim <- c(0,0)
        #C4.lim <- c(-10000,20000)
    
        resiFun <- function(strain, stress, time, temp, para) (strain - (para$C1 * stress ^ para$C2
                        * time ^ (para$C3 + 1) * exp(-1 * para$C4 / temp) * 1 / (para$C3 + 1)))

        # FITTING witout C4 --> only one temperature considered so far, in workbench creep material at mm is OK
        nls.out.TimeHardening <- nls.lm(par=para, fn = resiFun, strain = s.s.isochr[[4]]
                , time = s.s.isochr[[3]], stress = s.s.isochr[[2]], temp = s.s.isochr[[1]]+273, control = nls.lm.control(nprint=1)
                , lower = c(0,0,-1,C4.lim[1]), upper = c(1000,1000,0,C4.lim[2]))

            TimeHardpara <- para <- list(C1=nls.out.TimeHardening$par$C1,C2=nls.out.TimeHardening$par$C2,C3=nls.out.TimeHardening$par$C3,C4=nls.out.TimeHardening$par$C4)
            TimeHard_fun <- function(time, temp, stress, para) (para$C1 * stress ^ para$C2 * time ^ (para$C3 + 1) * exp(-1 * para$C4 / temp)* 1 / (nls.out.TimeHardening$par$C3 + 1))
            TimeHard_fun_sig <- function(time, temp, strain,para) ((strain * (para$C3 + 1)) / (para$C1 * time ^ (para$C3 + 1) * exp(-1 * para$C4 / temp)))^(1/para$C2)
                    

            s.s.isochr.org[ s.s.isochr.org$temp == creep.temp,]$fit <-  TimeHard_fun(s.s.isochr$time,s.s.isochr$temp+273,s.s.isochr$stress,para)
            TimeHardpara.df <- data.frame(temp=creep.temp,C1=nls.out.TimeHardening$par$C1,C2=nls.out.TimeHardening$par$C2,C3=nls.out.TimeHardening$par$C3,C4=nls.out.TimeHardening$par$C4)
            TimeHardpara.out <- rbind(TimeHardpara.out, TimeHardpara.df)



        }
        gg.creep <- ggplot(data=s.s.isochr.org,group=time) +
        geom_point(aes(x=strain,y=stress,color=as.factor(time),group=time)) +
        geom_line(aes(x=fit,y=stress,color=as.factor(time),group=time))+
        scale_colour_brewer(palette="Dark2")+
        labs(name="isochr. stress-Strain with fitting",x="Strain [mm/mm]",y="Stress [MPa]")+
        facet_wrap(.~temp)+
        scale_color_manual(name="time",labels=c(paste(round(as.numeric(levels(as.factor(unique(s.s.isochr.org$time))))/3600,0),"h,",round(as.numeric(levels(as.factor(s.s.isochr.org$time)))/(3600*24),1),"d",sep="")
        ),values=rev(hcl(h = seq(1,230,length.out=length(unique(s.s.isochr.org$time))), c = 100, l = 75)))+
        theme( axis.line = element_line( size = 1, linetype = "solid"))
        TimeHardpara.out <- TimeHardpara.out[-1,]
        ggsave(paste(str.mat.name,"/graphs/stress-strain-isocr-fit.png",sep=""),plot=gg.creep, width =40, height = 30, units="mm",scale=4)
        gg.creep

    }
}
## _______________________ critera for creep
if (evaluate.creep == TRUE && create.creep.critera == TRUE){
  s.s.isochr_org$stress <- s.s.isochr_org$stress * 1/A.stiff       # apply the global reduction factor on stress colum of creep data
  s.s.isochr_org$strain <- s.s.isochr_org$strain / 100         # convert % to mm/mm for creep data
  for (i in seq(1, length(unique(s.s.short$temp)),1)){

  temp <- unique(s.s.short$temp)[i]
  times <- c(1,10,1000,10000,100000)            # Get unique times in creep measurement (10,100, ...)
  df.cr.criteria <- data.frame(times)
  sig.cr.criteria <- c()
  sig.allw.cr.temp <- (max(subset(s.s.short, s.s.short[1] == temp)[[3]]) * A_S_creep) / A.str   # allwoable strain for creep
  epsilon.allw.cr.temp <- approxfun(subset(s.s.short, s.s.short[1] == temp)[[3]],subset(s.s.short, s.s.short[1] == temp)[[2]])(sig.allw.cr.temp)

   plot(subset(s.s.short, s.s.short[1] == temp)[[2]]           # Plot short term data for this temperature
   ,subset(s.s.short, s.s.short[1] == temp)[[3]]
    ,xlim = c(0, epsilon.allw.cr.temp*1.5)
    ,ylim = c(0, sig.allw.cr.temp*1.5)
    ,main = c('Stress strain isochron: ', str.mat.name, temp)
    ,xlab = 'Strain strain [mm/mm]'
    ,ylab = 'Stress [MPa]'
    ,type = 'o')
  grid()
 
 
   for (j in seq(1 , length(times),1)){              # for each measured time in the current creep curve of this temp
  lines(seq(0,epsilon.allw.cr.temp*1.1,length.out=100), TimeHard_fun_sig(times[j],temp+273,seq(0,epsilon.allw.cr.temp*1.1,length.out=100)), type = 'l')
  sig.cr.criteria <- rbind(sig.cr.criteria,  TimeHard_fun_sig(times[j],temp+273,epsilon.allw.cr.temp )) # create function of current time and insert allow.stress
   }
   abline(v=epsilon.allw.cr.temp,col='red')
   abline(h=sig.allw.cr.temp,col='red')
 

   par(mfrow=c(1,1))
 
   df.cr.criteria <- data.frame(data.frame(times), sig.cr.criteria)
   plot(log10(df.cr.criteria[[1]]), df.cr.criteria[[2]], xlim = c(0,6), ylim = c(0, max(df.cr.criteria[[2]])*1.1), xlab = 'time [h]', ylab = 'creep criterion [MPa]',type='p', main=paste("Creep Criterion for ", temp, " degC"))
   fit.cr.criteria <- lm(df.cr.criteria[["sig.cr.criteria"]] ~ log10(df.cr.criteria[["times"]]) +I(log10(df.cr.criteria[["times"]])^1))
   cr.criteria.fun <- function (t) fit.cr.criteria[["coefficients"]][1] + fit.cr.criteria[["coefficients"]][2] * t

   df.cr.criteria = data.frame(seq(0,7,1), cr.criteria.fun(seq(0,7,1)), func_23_strain(cr.criteria.fun(seq(0,7,1))))
   print(i)
   lines(seq(0,7,1), cr.criteria.fun(seq(0,7,1)), col ='red')
   abline(v=log10(1128),col='red')

   grid()
  }
}

## TESTINGS Viscoplasticity

epto <- s.s.short.true.A.stiff[s.s.short.true.A.stiff[[1]] == ref_temp, 2]
sigma <- s.s.short.true.A.stiff[s.s.short.true.A.stiff[[1]] == ref_temp, 3]
sigma_dif <- append(sigma, 0, after = 0)
sigma_dif <- sigma_dif[1: length(sigma_dif) -1]
delta_sigma <- sigma - sigma_dif
E_el <- temp.d.scalar.res[temp.d.scalar.res[[1]] == ref_temp, 3]
dotep_pl <- strain_rate - (delta_sigma * strain_rate) / (epto * E_el)
j <- 0
remove <- c()
for (i in dotep_pl){
 j <- j + 1
 if (i <= 0 || is.nan(i)){
  remove <- append(remove, j)
 }
}
dotep_pl <- dotep_pl[remove * -1]
sigma <- sigma[remove * -1]
epto <- epto[remove * -1]

sigma_y <- 1

## For Perzynda Viscoplasticity
# TB,RATE,1,1,2,PERZYNA
# TBDATA,1,gamma,m
if (create.perzynda.viscoplasticity == TRUE){
para <- list('gamma' = 1, 'm' = 3 )
resiFun <- function(sigma, para) (dotep_pl - para$gamma * ((sigma / sigma_y) -1) ^ (1/para$m))

nls.out.viscopl <- nls.lm(par=para, fn = resiFun, sigma = sigma, control = nls.lm.control(nprint=1))

fitted <- function(sigma) (nls.out.viscopl$par$gamma * ((sigma / sigma_y) -1) ^ (1/nls.out.viscopl$par$m))

par(mfrow = c(1,1))
plot(sigma, dotep_pl, type = 'p', main = 'plastic strain rate', ylim = c(0, max(dotep_pl)), xlim = c(0, max(sigma)))
lines(sigma, fitted(sigma), lwd = 2, col = 'red' )
grid()
}
 
##_______________________________________________________________________________ Fiber critera
## Strenght value estimation

df.eng.dens.adv.script <- df.eng.dens


if (length(df.eng.dens.adv.script$temp) < 6){
    temp.index.to.use <- 1:6
    no.dummy <- 6 - length(df.eng.dens.adv.script$temp)
    for (k in 1:no.dummy) {
        df.eng.dens.adv.script <- rbind(df.eng.dens.adv.script,df.eng.dens.adv.script[length(df.eng.dens.adv.script[,1]),])
        df.eng.dens.adv.script[length(df.eng.dens.adv.script[,1]),]$temp <- df.eng.dens.adv.script[length(df.eng.dens.adv.script[,1]),]$temp + 0.1
    }

}
#parabolic.crit.fun.ehrenstein(sigma.allw.true, sigma.allw.true, sigma.allw.true, m=1.3)
out.vec.apdl <- c(paste('! Be carefull, this values are generated with a Fiber-Stiffness factor of:',A.stiff.F,sep=' '))

header <-(readLines(file(paste(getwd(),'/01_Hilfstabellen und Rechnungen/APDL/Anisotropic-BISO_Header.mac',sep=''),open="r")))
out.vec.apdl <- append(out.vec.apdl,header)



if (evaluate.creep) {
out.vec.apdl <- append(out.vec.apdl, paste('*DIM,creep_coeff,ARRAY,',length(TimeHardpara.out$temp),',5,1, , ,\n\n',sep=''))
out.vec.apdl <- append(out.vec.apdl, paste('num_temp_creep=',length(TimeHardpara.out$temp),'\t\t\t! number of creep temperatures\n\n',sep=''))

    for (i in 1:length(TimeHardpara.out$temp)){
    out.vec.apdl <- append(out.vec.apdl,(paste('*SET,creep_coeff(',i,',1,1) , ',TimeHardpara.out$temp[i],'\t\t\t!temp',sep='')))
    out.vec.apdl <- append(out.vec.apdl,(paste('*SET,creep_coeff(',i,',2,1) , ',TimeHardpara.out$C1[i],'\t\t\t!C1',sep='')))
    out.vec.apdl <- append(out.vec.apdl,(paste('*SET,creep_coeff(',i,',3,1) , ',TimeHardpara.out$C2[i],'\t\t\t!C2',sep='')))
    out.vec.apdl <- append(out.vec.apdl,(paste('*SET,creep_coeff(',i,',4,1) , ',TimeHardpara.out$C3[i],'\t\t\t!C3',sep='')))
    out.vec.apdl <- append(out.vec.apdl,(paste('*SET,creep_coeff(',i,',5,1) , ',TimeHardpara.out$C4[i],'\t\t\t!C4',sep='')))
    out.vec.apdl <- append(out.vec.apdl, c('!-------------------\n\n\n'))
    }


}
out.vec.apdl <- append(out.vec.apdl, paste('*DIM,ref_material,ARRAY,',length(temp.index.to.use),',7,1, , ,\n\n',sep=''))

# create output set per tem
j <- 1
for (i in temp.index.to.use){
   out.vec.apdl <- append(out.vec.apdl,(paste('*SET,ref_material(',j,',1,1) , ',df.eng.dens.adv.script$temp[i],'\t\t\t!temp',sep='')))
   out.vec.apdl <- append(out.vec.apdl,(paste('*SET,ref_material(',j,',2,1) , ',df.eng.dens.adv.script$e.el[i],'\t\t\t!E-Modulus-Datasheet MPa',sep='')))
   out.vec.apdl <- append(out.vec.apdl,(paste('*SET,ref_material(',j,',3,1) , ',df.eng.dens.adv.script$poisson[i],'\t\t\t!Poisson',sep='')))
   out.vec.apdl <- append(out.vec.apdl,(paste('*SET,ref_material(',j,',4,1) , ',df.eng.dens.adv.script$yield[i],'\t\t\t!yield stress with fiber from datasheet',sep='')))
   out.vec.apdl <- append(out.vec.apdl,(paste('*SET,ref_material(',j,',5,1) , ',df.eng.dens.adv.script$m.yield[i],'\t\t\t!matrix yield strength',sep='')))
   out.vec.apdl <- append(out.vec.apdl,(paste('*SET,ref_material(',j,',6,1) , ',df.eng.dens.adv.script$strength.datasheet[i],'\t\t\t!Composite ultimite strength, from critera',sep='')))
   out.vec.apdl <- append(out.vec.apdl,(paste('*SET,ref_material(',j,',7,1) , ',df.eng.dens.adv.script$tan.mod.biso[i]/df.eng.dens.adv.script$e.el[i],'\t\t\t!Datasheet Tangent Modulus ratio for BISO',sep='')))
   out.vec.apdl <- append(out.vec.apdl, c('!-------------------\n\n\n'))
   j <- j + 1
   }

footer <-(readLines(file(paste(getwd(),'/01_Hilfstabellen und Rechnungen/APDL/Anisotropic-BISO_Footer_210420.mac',sep=''),open="r")))
out.vec.apdl <- append(out.vec.apdl,footer)
 
## Export apdl File
filename <- paste(str.mat.name,"/Fiber-Files/",str.mat.name,"-adpd-adv-fiber.txt", sep="")
writeLines(as.character(out.vec.apdl), con=filename,sep='\n')


## Export converse material


out.vec.conv.mat <- c('1002,2','DEFAULT,0.3,30,1,0.7,0.083','23,2000,0.34,1.4,0.00011,73000,0.22,2.6,5.4e-006,0.1,2,0,0')
out.vec.conv.mat <- append(out.vec.conv.mat,paste(str.mat.name,fiber.w.perc/100,ar,length(temp.index.to.use),0.7,0.083,sep=','))


# store temperature lines in out.vec.conv.mat
for (i in temp.index.to.use){
   out.vec.conv.mat <- append(out.vec.conv.mat,(paste(df.eng.dens$temp[i],df.eng.dens$e.m[i],
    df.eng.dens$poisson[i],df.eng.dens$rho.m[i]/1000,df.eng.dens$cte.m[i],E.f,nu.f,rho.f/1000,cte.f,0.1,2,0,0,sep=',')))
}


## Export conversemat File
filename <- paste(paste(str.mat.name,"/Fiber-Files/",sep=""),str.mat.name,"_converse",".mdb", sep="")
writeLines(as.character(out.vec.conv.mat), con=filename,sep='\n')
temp.d.scalar.res
#print(fat.strain.cyclic.ep.max.limit.function)
print(df_prony)







## TTSP for thermal cycling loads
cte.iso.fun <- approxfun(df.eng.dens$temp,df.eng.dens$cte.iso)

#cte.iso.fun <- function (x){ cte.iso.fun(x) - df.eng.dens[df.eng.dens$temp == 23,]$cte.iso}



if (fiber.w.perc > 0){
    E.LT.extrapol <- (((diff(df.eng.dens$E.iso)[1] *-1) / diff(df.eng.dens$temp)[1]) * 60) +df.eng.dens$E.iso[1]
    cte.LT.extrapol <- (((diff(df.eng.dens$cte.iso)[1] *-1) / diff(df.eng.dens$temp)[1]) * 60) +df.eng.dens$cte.iso[1]
    cte.iso.fun <-approxfun(c(-100,df.eng.dens$temp),c(cte.LT.extrapol,df.eng.dens$cte.iso))
    e.iso.fun <- approxfun(c(-100,df.eng.dens$temp),c(E.LT.extrapol,df.eng.dens$E.iso))
}
if (fiber.w.perc == 0){
    E.LT.extrapol <- (((diff(df.eng.dens$e.m)[1] *-1) / diff(df.eng.dens$temp)[1]) * 60) +df.eng.dens$e.m[1]
    cte.LT.extrapol <- (((diff(df.eng.dens$cte.m)[1] *-1) / diff(df.eng.dens$temp)[1]) * 60) +df.eng.dens$cte.m[1]
    cte.iso.fun <-approxfun(c(-100,df.eng.dens$temp),c(cte.LT.extrapol,df.eng.dens$cte.m))
    e.iso.fun <- approxfun(c(-100,df.eng.dens$temp),c(E.LT.extrapol,df.eng.dens$e.m))
}



index.fun <- function (temp,cycle) {((cte.iso.fun(temp)- df.eng.dens[df.eng.dens$temp == 23,]$cte.iso)*e.iso.fun(temp))/sig.bk.gerber.fun(0,temp,cycle)}
index.fun <- function (temp,cycle) {abs(((cte.iso.fun(temp))*(abs(temp-23))*e.iso.fun(temp))/sig.bk.gerber.fun(0,temp,cycle))}

temp.from.index.fun <- function (index,cycle) {}

#nls(index.fun(temp,1500)-index.fun(85,3000))
#require(reconPlots)
#lm(0~index.fun(temp,1500)-index.fun(85,3000))

#install.packages("remotes")
#remotes::install_github("andrewheiss/reconPlots")
#library(reconPlots)


#curve_intersect(line1, line2)
start.temp <- 85
start.cycles <- 3000
speed.fac <- 2
temp.add.offs.pos <- 12
temp.add.offs.neg <- 15
temp.use <- seq(start.temp-temp.add.offs.neg,start.temp+temp.add.offs.pos,length.out=200)
#temp.use <- seq(-50,120,length.out=200)

cycles <- c(start.cycles,start.cycles/speed.fac,100,50)

df.thermal.index <- data.frame(cycle=NA,temp=NA,index=NA)
for (i in cycles) {
    for (j in temp.use) {
    df.thermal.index <- rbind(df.thermal.index,data.frame(cycle=i,temp=j,index=index.fun(j,i)))
    }
}

ggplot(df.thermal.index) + theme_bw() + geom_line(aes(x=temp,y=index,group=cycle,color=cycle)) +
geom_point(aes(x=start.temp,y=index.fun(start.temp,start.cycles))) + geom_hline(yintercept=index.fun(start.temp,start.cycles),linetype = "dashed") +
geom_point(aes(x=start.temp,y=index.fun(start.temp,100))) + geom_hline(yintercept=index.fun(start.temp,100),linetype = "dashed") +

scale_x_continuous(breaks = round(seq(start.temp-temp.add.offs.neg, start.temp+temp.add.offs.pos, by = 1),1))

#+  scale_y_continuous(breaks = round(seq(start.temp-5, start.temp+temp.add.offs.pos, by = 0.5),2)) + xlim(start.temp-5, start.temp+temp.add.offs.pos)



writeClipboard(gerber.string.N.from.sig.reduced)
writeClipboard(workbench.SF)
writeClipboard(workbench.index.datasheet)


##Mat Card exports
#Sig allow
write.table(df.sig.allow[,c(1,2,3,4,5,6,9)], "clipboard-16384", sep="\t", col.names=FALSE,row.names = FALSE)
#Ep allow
write.table(df.ep.allow[,c(1,2,3,4,5,6,9)], "clipboard-16384", sep="\t", col.names=FALSE,row.names = FALSE)
#crit values
write.table(cbind(df.ep.allow[,c(1,8)],df.sig.allow[,8]), "clipboard-16384", sep="\t", col.names=FALSE,row.names = FALSE)

##_____________________________________________________________________________________________________________________ Assemble output Excel file
OutsideBorders <-
  function(wb_,
           sheet_,
           rows_,
           cols_,
           border_col = "black",
           border_thickness = "medium") {
    left_col = min(cols_)
    right_col = max(cols_)
    top_row = min(rows_)
    bottom_row = max(rows_)
    
    sub_rows <- list(c(bottom_row:top_row),
                     c(bottom_row:top_row),
                     top_row,
                     bottom_row)
    
    sub_cols <- list(left_col,
                     right_col,
                     c(left_col:right_col),
                     c(left_col:right_col))
    
    directions <- list("Left", "Right", "Top", "Bottom")
    
    mapply(function(r_, c_, d) {
      temp_style <- createStyle(border = d,
                                borderColour = border_col,
                                borderStyle = border_thickness)
      addStyle(
        wb_,
        sheet_,
        style = temp_style,
        rows = r_,
        cols = c_,
        gridExpand = TRUE,
        stack = TRUE
      )
      
    }, sub_rows, sub_cols, directions)
  }



##_____________________________________________________________________________________________________________________ Assemble output Excel file

wb <- createWorkbook()
#white.bg <- createStyle(bgFill=rgb(1, 1, 1))

rgb(1, 1, 1)

addWorksheet(wb, "MaterialCard")
#addStyle(wb, "MaterialCard", white.bg, rows=1:200, cols=1:200,stack=TRUE,gridExpand = TRUE)


# create main outer frame
OutsideBorders(wb, "MaterialCard", rows_=2:15,cols_=2:15)

writeData(wb, "MaterialCard", df.general.info, startRow = 3, startCol = 3,borders = "all",colNames=FALSE)
mergeCells(wb, "MaterialCard", cols = 3:5, rows = 2)
writeData(wb, "MaterialCard", "General Info", startRow = 2, startCol = 3,borders = "all",colNames=FALSE)

header.style <- createStyle(  fontSize = 14, fontColour = rgb(0, 0, 0),  textDecoration = c("bold"),  
                    halign = "center", valign = "center", 
                    border =  c("top", "bottom", "left", "right"),borderStyle="thin")
         
addStyle(wb, "MaterialCard", header.style, rows=2, cols=3:5, gridExpand = FALSE, stack = TRUE)

writeData(wb,  "MaterialCard", df.sig.allow[,c(1,2,3,4,5,6,9)], startRow = 14, startCol = "F",borders = "all")
setColWidths(wb,  "MaterialCard", cols= 2:200,  widths = "auto",  ignoreMergedCells = FALSE)



# Frame for first Graph
#writeData(wb, "MaterialCard", matrix("",nrow=13,ncol=6), startRow = 2, startCol = 7,borders = "surrounding",colNames=FALSE)
insertImage(  wb,  "MaterialCard",  paste(str.mat.name,"/graphs/stress-strain.png",sep=""),  width =8.8, height = 6.6, units="cm",  startRow = 2,  startCol = "G",  dpi = 300)

insertImage(  wb,  "MaterialCard",  paste(str.mat.name,"/graphs/cte-temp.png",sep=""),  width =8.8, height = 6.6, units="cm",  startRow = 2,  startCol = "I",  dpi = 300)


addWorksheet(wb, "elastoplast")
if (fiber.w.perc > 0) { writeData(wb, "elastoplast", data.frame(temp=df.eng.dens$temp,E.iso=df.eng.dens$E.iso,poisson=df.eng.dens$poisson), startRow = 1, startCol = 1,borders = "all")}
if (fiber.w.perc == 0) { writeData(wb, "elastoplast", data.frame(temp=df.eng.dens$temp,E.iso=df.eng.dens$e.m,poisson=df.eng.dens$poisson), startRow = 1, startCol = 1,borders = "all")}

writeData(wb, "elastoplast", s.s.short.true.pl, startRow = 1, startCol = 5,borders = "all")

if (evaluate.creep) {
    addWorksheet(wb, "creep")
    writeData(wb, "creep", TimeHardpara.out, startRow = 1, startCol = 1,borders = "all")
}

addWorksheet(wb, "CLTE")
if (fiber.w.perc > 0) {writeData(wb, "CLTE", data.frame(temp=df.eng.dens$temp,cte.iso=df.eng.dens$cte.iso), startRow = 1, startCol = 1,borders = "all")}
if (fiber.w.perc == 0) {writeData(wb, "CLTE", data.frame(temp=df.eng.dens$temp,cte.iso=df.eng.dens$cte.m), startRow = 1, startCol = 1,borders = "all")}


addWorksheet(wb, "limits")
writeData(wb, "limits", df.sig.allow[,c(1,2,3,4,5,6,9)], startRow = 4, startCol = 1,borders = "all")
writeData(wb, "limits", df.ep.allow[,c(1,2,3,4,5,6,9)], startRow = 4, startCol = 9,borders = "all")
writeData(wb, "limits", cbind(df.ep.allow[,c(1,8)],df.sig.allow[,8]), startRow = 4, startCol = 18,borders = "all")
writeData(wb, "limits", c("Isotropic-SI Formula",workbench.index.datasheet), startRow = 1, startCol = 1,borders = "surrounding",colNames=FALSE)



addWorksheet(wb, "overview")
writeData(wb, "overview", df.eng.dens, startRow = 1, startCol = 1,borders = "all")

if (fiber.w.perc > 0) {

    addWorksheet(wb, "Lin.Matrix.prop")
    df.temp <- df.eng.dens[,c(1,19,12,18,20)]
    df.temp$k <- k.m
    df.temp$c <- c.m
    writeData(wb, "Lin.Matrix.prop", paste(str.mat.name,"_matrix",sep=""), startRow = 1, startCol = 1)

    writeData(wb, "Lin.Matrix.prop", df.temp, startRow = 3, startCol = 1,borders = "all")
    temps.string <- ""
    for (i in unique(df.s.s.scale$temp)){temps.string <- paste(temps.string,i,sep=",")}
    
    writeData(wb, "Lin.Matrix.prop", temps.string, startRow = 2, startCol = 1)

    #write orientation depended stress-strain
    addWorksheet(wb, "S-S-Orient")
    index <- 1
    for (temp in unique(df.s.s.scale$temp)){       
            
            # crate data.frame for 0 and 90 deg orientation and scale values from s.s.short.true according scale factors
            df.0 <- data.frame(strain=s.s.short.true[s.s.short.true$temp == temp,]$strain / df.s.s.scale[df.s.s.scale$temp == temp,]$scale.n,
                                stress=s.s.short.true[s.s.short.true$temp == temp,]$stress * df.s.s.scale[df.s.s.scale$temp == temp,]$scale.n)
            df.90 <- data.frame(strain=s.s.short.true[s.s.short.true$temp == temp,]$strain / df.s.s.scale[df.s.s.scale$temp == temp,]$scale.p,
                                stress=s.s.short.true[s.s.short.true$temp == temp,]$stress * df.s.s.scale[df.s.s.scale$temp == temp,]$scale.p)
            
            # Extract plastic portion for both directions
            df.0 <- lin.el.extract(df.0$strain, df.0$stress, temp)[, c("temp", "pl.strain", "stress")][-1]
            df.90 <- lin.el.extract(df.90$strain, df.90$stress, temp)[, c("temp", "pl.strain", "stress")][-1]

            # remove first Zero plastic strain row
            df.0 <- df.0[!df.0$pl.strain == 0,]
            df.90 <- df.90[!df.90$pl.strain == 0,]
            
            df.0[1,1] <- 0
            df.90[1,1] <- 0
            
            # export to csv
            
            mergeCells(wb, "S-S-Orient", cols = (index):(index+1), rows = 1)  # mere header cells
            writeData(wb, "S-S-Orient", paste("0deg_",temp,"C",sep=""), startRow = 1, startCol = (index),borders = "all") # add title
            writeData(wb, "S-S-Orient", df.0, startRow = 2, startCol = (index),borders = "all")
            index <- index + 3
            mergeCells(wb, "S-S-Orient", cols = (index):(index+1), rows = 1)  # mere header cells
            writeData(wb, "S-S-Orient", paste("90deg_",temp,"C",sep=""), startRow = 1, startCol = (index),borders = "all") # add title
            writeData(wb, "S-S-Orient", df.90, startRow = 2, startCol = (index),borders = "all")
       
             index <- index + 3
        }
    
    }




#insertImage(  wb,  "elasticity",  "stress-strain.png",  width =40, height = 30, units="cm",  startRow = 2,  startCol = 5,  dpi = 300)
saveWorkbook(wb, file = paste(str.mat.name,"/",str.mat.name,"--Output.xlsx",sep=""), overwrite = TRUE)

## Sine Vibration Test Input Plotter / Analyser for FEA

library(reshape2)
library(ggplot2)
library('minpack.lm')

rm(list=ls(all=TRUE))

# For Windows Systems
setwd('')

source('criteria.lib_210917.R')

## _______________________________________________________________________________________  General Data
oct.per.min <- 1
f1 <- 50					# Hz, Starting Frequency of range
t <- seq(0,4*60,1)				# s



f <- exp(log(2)*t/60*oct.per.min + log(f1))	# Hz, resulting frequency after time t

#Plotting Frequency increase over Time
{
  plot.df <- data.frame(x=t,freq=f)

  #x.column <- 'temp.degc'
  x.column <- names(plot.df)[1]
  y.exclude <- c('Tg.ratio')  # columns which should be not plotted
  plot.df <- plot.df[!(names(plot.df) %in% y.exclude)]
  plot.df <- melt(plot.df, id=x.column); names(plot.df)[1] <- 'x'


  gg <- ggplot(data=plot.df,aes(group=variable,color=variable)) +geom_line(aes(x=x,y=value),size=0)
  gg <- (gg + theme_bw()
	+ labs(title='Frequency over Time', x='t [s]', y='Frequency [Hz]')
	#+ scale_x_continuous(breaks = seq(-200, 200,20))
	#+ scale_y_continuous(breaks = seq(-200, 200,20))
	#+ scale_color_brewer(name='versions',palette='Set2')
	#+ facet_grid(variable~., scales='free')
	#+ xlim(-6,6)
	#+ ylim(-6,6)
  )
  print(gg)
}

sine.sweep.oct.min.from.time <- function (f.min=10,f.max=1000,t=20) {
  # octave speed according time unitinput, so min gives minutes as result
  vs <- 1/t * (log((f.max)/(f.min)))/log(2)
  return(vs)
}

oct.p.min <- sine.sweep.oct.min.from.time(f.min=10,f.max=1000,t=20)
nocycles <- no.cycles.sine.sweep.counting(oct.per.min=oct.p.min,f.min=0,f.max=seq(0,300,10),no.sweeps=144*(60/20))

log10(nocycles)

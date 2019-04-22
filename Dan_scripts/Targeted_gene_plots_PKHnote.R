load('Sleuth_exp.RData')
# Load the data from my computer:
load("C:/Users/Po-Kai/Box Sync/Mimulus_swc_timecourse_2018/RNAseq/limma/Dan_scripts/Sleuth_exp.RData")

library(cowplot)
library(reshape)


# Select the data (Col-FRI only) + drop the unused levels
sub_sample_info = droplevels(subset(sub_sample_info,Genotype == 'ColFRI' & Sampling.Day == '2014-02-26'))
# The samples
exp = exp[,sub_sample_info$sample]


# I don't have this file.
temps = read.delim('~/Documents/Arabidopsis/Lianas_experiment/Bolting/Data/ACE_temp_profile.txt')
ACE_22VarLD = subset(temps,Treatment == '22VarLD')[,-1]
# Create a new empty column: "Hour" 
ACE_22VarLD$Hour = 0
# Add a row: 12, 0, 0; row name: "Temp", "Duration", "Hour"
ACE_22VarLD = rbind(data.frame(Temp=12,Duration=0,Hour=0),ACE_22VarLD)
# Fill the "Hour" column by adding up the "Duration"
ACE_22VarLD$Hour = cumsum(ACE_22VarLD$Duration)
# Delete the first row
ACE_22VarLD_shift = ACE_22VarLD[-1,]
# Set the "Hour" column in the new data.frame: delete the last row of the original data.frame
ACE_22VarLD_shift$Hour = ACE_22VarLD$Hour[-nrow(ACE_22VarLD)]
# rbind the two data.frame
ACE_22VarLD = rbind(ACE_22VarLD,ACE_22VarLD_shift)
ACE_22VarLD = ACE_22VarLD[order(ACE_22VarLD$Hour),]


# Create a function called "se" for calculating the standard error
se = function(x) {
	x = na.omit(x)
	sd(x)/sqrt(length(x))
}
# These two functions are prepared for transformation (but Dan didn't use them eventually). 
transform_y = function(y,fun) fun(y)
itransform_y = function(y,fun) fun(y)


# The function for the plot?!
timecourse_plot = function(gene,transcript,cols = c('red','blue'),fun=function(x) x,ifun = function(x) x){
	# recover()
	# Extract the expression level from the data set
  y = exp[transcript,]
  
  # Working space in sub_sample_info:
	sub_sample_info$y = y
	
	# Extract the data at ZT==0 and duplicate them. The duplicates are called ZT=24
	sub_sample_info$ZT = as.numeric(as.character(sub_sample_info$ZT))
	sub_sample_info_ZT24 = sub_sample_info[sub_sample_info$ZT==0,]
	sub_sample_info_ZT24$ZT = 24
	sub_sample_info = rbind(sub_sample_info,sub_sample_info_ZT24)
	
	# use the "fun" set up previously, here we can call the "y" in the "sub_sample_info" data.frame
	# Calculate the mean and se for the samples (HOW does this loop???)
	means = tapply(fun(sub_sample_info$y),list(sub_sample_info$ZT,sub_sample_info$Treatment),mean)
	ses = tapply(fun(sub_sample_info$y),list(sub_sample_info$ZT,sub_sample_info$Treatment),se)
	means = melt(means)
	colnames(means) = c('ZT','Treatment','mean')
	ses = melt(ses)
	colnames(ses) = c('ZT','Treatment','se')
	mean_data = data.frame(means,ymin = means$mean-2*ses$se,ymax = means$mean+2*ses$se,stringsAsFactors=F)
	mean_data$Treatment = as.character(mean_data$Treatment)
	mean_data$mean = ifun(mean_data$mean)
	mean_data$ymin = pmax(0,ifun(mean_data$ymin))
	mean_data$ymax = ifun(mean_data$ymax)
	# mean_data$ZT = factor(mean_data$ZT)
	
	# recover()
	
	# For the expression of Light and Dark:
	light = data.frame(xmin = c(0,16),xmax = c(16,23.9),ymin=c(0,0),ymax = rep(max(sub_sample_info$y,mean_data$ymax),2)-.1,light=factor(c(1,0)),Hour=1,Temp=1)
	
	
	mean_data = subset(mean_data,Treatment %in% unique(Treatment)[!is.na(cols)])
	
	p = ggplot(sub_sample_info) + ggtitle(paste(gene))
	p = p + geom_rect(data=light,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=light)) + scale_fill_manual(values=c('0'='grey50','1'=NA))
	p = p + geom_point(aes(x = ZT,y=y,color = Treatment)) + scale_color_manual(values=cols)
	p = p + geom_ribbon(data=mean_data,aes(x=ZT,y=mean,ymin = ymin,ymax=ymax, linetype=NA,group = Treatment),alpha = .2)
	p = p + geom_line(data=mean_data,aes(x=ZT,y=mean,color=Treatment)) 
	p = p + xlim(c(0,24)) + ylim(c(0,NA)) + theme(legend.position="none") 
	# ACE_22VarLD$Temp = ACE_22VarLD$Temp/max(ACE_22VarLD$Temp) * max(light$ymax)
	# p = p + geom_line(data = ACE_22VarLD,aes(x=Hour,y=Temp),color='blue')
	class(p) = 'ggplot'
	return(p)		
}	


# I don't have this file.
gene_list = read.csv('RNAseq_analysis/Gene_list.csv',stringsAsFactors=F)

# Set up the target genes:
genes = list(Clock = c(
						LHY = 'AT1G01060.1',
						TOC1 = 'AT5G61380.1',
						PRR5 = 'AT5G24470.1',
						PRR7 = 'AT5G02810.1'
						),
			Clock2 = c(
						ELF3 = 'AT2G25930.1',
						ELF4 = 'AT2G40080.1',
						LUX = 'AT3G46640.1',
						RVE8 = 'AT3G09600.1'
						),
			Clock3 = c(
						CAB2 = 'AT1G29920.1',
						SPL9 = 'AT2G42200.1',
						TOE1 = 'AT1G01060.1',
						HY5 = 'AT5G11260.1'
						),
			Photoperiod = c(
						GI = 'AT1G22770.1',
						CO = 'AT5G15840.1',
						FKF1 = 'AT1G68050.1',
						PIF4 = 'AT2G43010.1'
						),
			Norm = c(
						TUB5 = 'AT1G20010.1',
						ACT2 = 'AT3G18780.2',
						UBC = 'AT1G14400.1',
						PP2A = 'AT1G69960.1'
						),
			MADS = c(
						FLC = 'AT5G10140.1',
						SVP = 'AT2G22540.1',
						FLM_b = 'AT1G77080.4',
						FLM_d = 'AT1G77080.2'
						),
			Temp = c(
						HSP70 = 'AT3G12580.1',
						Unk3 = 'AT5G13930.1',
						Unk1 = 'AT5G52640.1',
						Unk2 = 'AT1G05680.1'
						)
			)

pdf('Targeted_timecourse_plots.pdf')


# These are the dulplicates of the previous setting:
#load('Sleuth_exp.RData')
#sub_sample_info = droplevels(subset(sub_sample_info,Genotype == 'ColFRI' & Sampling.Day == '2014-02-26'))
#exp = exp[,sub_sample_info$sample]


# Set up the color function (but why?)
# Anyway, the colors Dan chose are good. Just use them!
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
# These are checking the colors (will be used in the for loop below.)
cols = gg_color_hue(2)
cols = c(gg_color_hue(2)[1],NA)


# Plots!
for(group in genes){
	plots = list()
	cols = gg_color_hue(2)
	for(i in 1:length(group)){
		plots[[i]] = timecourse_plot(names(group)[i],group[[i]],cols)
	}
	print(plot_grid(plotlist = plots))
	cols = c(gg_color_hue(2)[1],NA)
	for(i in 1:length(group)){
		plots[[i]] = timecourse_plot(names(group)[i],group[[i]],cols)
	}
	print(plot_grid(plotlist = plots))
}


###### The other plots:

load('Sleuth_exp.RData')

sub_sample_info = droplevels(subset(sub_sample_info,Genotype == 'Col' & Sampling.Day == '2014-02-26'))
exp = exp[,sub_sample_info$sample]


# Photoperiod only:
group = genes$Photoperiod
plots = list()
	cols = gg_color_hue(2)
	print(timecourse_plot('FT','AT1G65480.1',cols))
	cols = c(gg_color_hue(2)[1],NA)
	print(timecourse_plot('FT','AT1G65480.1',cols))
dev.off()


# Plots for clock genes:
clock_genes = c(
						LHY = 'AT1G01060.1',
						CCA1 = 'AT2G46830.2',
						PRR9 = 'AT2G46790.2',
						PRR7 = 'AT5G02810.1',
						PRR5 = 'AT5G24470.1',
						TOC1 = 'AT5G61380.1',
						LUX = 'AT3G46640.1',
						ELF4 = 'AT2G40080.1',
						ELF3 = 'AT2G25930.1',
						RVE8 = 'AT3G09600.1'
						# RVE4 = 'AT5G02840.1'
						)
plots = list()
for(i in 1:length(clock_genes)){
	cols = gg_color_hue(2)
	for(i in 1:length(clock_genes)){
		plots[[names(clock_genes)[i]]] = timecourse_plot(names(clock_genes)[i],clock_genes[[i]],cols)
	}	
}
with(plots,plot_grid(NULL,NULL,RVE8,NULL,NULL,
		  LHY,PRR5,TOC1,LUX,ELF4,
		  CCA1,PRR9,PRR7,ELF3,NULL,nrow=3))

#install.packages('patchwork')
library( ggplot2 )
library( patchwork )

visorirs <- c('IR','VIS')

# DEFINE CONSTANTS ETC
gal  <- 'SMC'
path <- '/home/aschoot/Desktop/etteren_projects'
DM   <- 18.977 # Distance modulus from graczyk+2020
AV   <- 0.35   # Extinction in mag from Schootemeijer+2021

# READ DATA
d <- read.csv( file.path(path,'/tidy_lrsg_vs_z/lrsg_meth/data/yang19_x_d18smc.csv') )

# DEFINE FUNCTIONS
get_props <- function( visorir ) {              # data and extinction are different for IR and VIS data, we do this so we can use a loop later
	if (visorir =='IR') {
		cb_poly <- c(5.00513314,-18.39883674,25.83141869,-18.2448306,8.62224833,-0.19544753)  # highest to lowest order
		mags    <- d$Kmag
		colors  <- d$Jmag - d$Kmag
		AbdivAv <- 0.243                # extinction in blue band divided by extinction in V band
		ArdivAv <- 0.078
	}

	if (visorir == 'VIS') {
		cb_poly <- c(0.12977808, -1.07136445, 3.16558804, -3.92916783, 1.26515746, 0.1510785) # highest to lowest order
		mags    <- d$imagSkyM
		colors  <- d$rmagSkyM - d$imagSkyM
		AbdivAv <- 0.843
		ArdivAv <- 0.628
	}
	return ( list(cb_poly=cb_poly, mags=mags, colors=colors, AbdivAv=AbdivAv, ArdivAv=ArdivAv ) )
}

get_logls <- function( mags,colors,DM,AV,AbdivAv,ArdivAv,cb_poly ){
	Ablue   <- AV*AbdivAv
	Ared    <- AV*ArdivAv
	mabs0   <- mags - Ared - DM
	colors0 <- colors - (Ablue-Ared)
	BCs     <- cb_poly[1]*colors0^5 + cb_poly[2]*colors0^4 + cb_poly[3]*colors0^3 + cb_poly[4]*colors0^2 + cb_poly[5]*colors0 + cb_poly[6]
	logls   <- -0.4*( mabs0 - 4.74 + BCs )
	return (logls)
}

# START CALCULATING THINGS
'Column names:'
print( names(d) )
for (i in 1:2) {
	visorir <- visorirs[i]

	results <- get_props( visorir )
	# extract mags etc from results, since R does not return multiple variables
	mags    <- results$mags
	colors  <- results$colors
	cb_poly <- results$cb_poly
	AbdivAv <- results$AbdivAv
	ArdivAv <- results$ArdivAv

	logls <- get_logls( mags,colors,DM,AV,AbdivAv,ArdivAv,cb_poly ) # this line gets the results
	#print( sort(logls,decreasing=TRUE) )
	new_col_name = paste( 'logL', visorir, sep='_' )
	print( new_col_name )
	d[[new_col_name]] <- logls                                      # add calculated logls to dataframe
}

# GET THE STDs of logL values calculated with IR and VIS data
std_IR  <- sd( d$logL-d$logL_IR )
std_VIS <- sd( d$logL-d$logL_VIS, na.rm=TRUE )
IR_text <- paste( 'SMC, IR\n std =',  round(std_IR,2 ) )
VIS_text<- paste( 'SMC, VIS\n std =', round(std_VIS,2) )

# OK NOW WE PLOT
# (in 2 panels)
#print(d[ ,6:ncol(d)])  # show new columns
p1 <- ggplot( data=d, aes( x=logL, y=logL_IR  ) ) +                     # left panel IR data
	geom_point( color="#ff7f0e" )             +
	annotate( "text", x=4.6, y=5.4, label=IR_text )
p2 <- ggplot( data=d, aes( x=logL, y=logL_VIS ) ) +                     # right panel VIS data
	geom_point( color="#ff7f0e" )             +
	annotate( "text", x=4.6, y=5.4, label=VIS_text )
p <- p1 + p2

ggsave( "compaR.png", plot=p, width=8, height=4 )

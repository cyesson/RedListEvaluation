#####RedListEvaluation - Perform red list evalation from point observation data

####required packages
library(dplyr)
library(ggplot2)
library(ggspatial)
library(raster)
library(red)
library(sf)
library(rnaturalearth)
library(stats)
library(stringr)
#' @import dplyr 
#' @import ggplot2
#' @import ggspatial
#' @import red
#' @importFrom stats median quantile lm coefficients predict
#' @importFrom sf sf_project st_crs
#' @importFrom raster rasterFromXYZ rasterToPolygons buffer disaggregate
#' @importFrom rnaturalearth ne_countries
#' @importFrom stringr str_replace

################################################################
# Helper functions
################################################################

################################################################
# get final evaluation (taking the worst criterion)
################################################################
GetFinalEvaluation<-function(dl){

    RedList.Category<-c("CR","EN","VU","NT","LC", "DD")

    # find worst of the criteria
    dl.cat<-data.frame(Criteria=c("A2b", "A2c", "A2c", "B1", "B2", "D2"), Category=NA)
    dl.cat$Category[1]<-dl$A2b$A2b
    dl.cat$Category[2]<-dl$A2c$AOO$A2c
    dl.cat$Category[3]<-dl$A2c$EOO$A2c
    dl.cat$Category[4]<-dl$B$B1
    dl.cat$Category[5]<-dl$B$B2
    dl.cat$Category[6]<-dl$D2$D2

    # factorise category based on full set of categories
    dl.cat$Category<-factor(dl.cat$Category, levels=RedList.Category)
    dl.cat.worst<-which(as.integer(dl.cat$Category) == min(as.integer(dl.cat$Category)))

    final.category<-dl.cat$Category[dl.cat.worst[1]]

    if(final.category != "DD"){

        # string final criteria together
        final.criteria <- paste(dl.cat$Criteria[dl.cat.worst], collapse="; ")

        # fix for certain combinations
        final.criteria <- stringr::str_replace(final.criteria, "A2c; A2c", "A2c")
        final.criteria <- stringr::str_replace(final.criteria, "A2b; A2c", "A2bc")
        final.criteria <- stringr::str_replace(final.criteria, "B1; B2", "B12ab(ii)")
        final.criteria <- stringr::str_replace(final.criteria, "B1", "B1ab(ii)")
        final.criteria <- stringr::str_replace(final.criteria, "B2", "B2ab(ii)")

        # if LC remove D2 from the end of the string
        if(final.category=="LC"){
            final.criteria <- stringr::str_replace(final.criteria, "; D2", "")
            # if this is LC only based on D2 then its not really an evaluation
            if(final.criteria=="D2"){
                final.category<-"DD"
                final.criteria==""
            }
        }

        final.category <- paste(final.category, final.criteria)

    }

    return(final.category)
    
}


################################################################
# internal function to count tetrads
# find how many contiguous areas of tetrads are in the distribution
################################################################
EvaluateTetradBlocks <- function(df){
    # find unique 2km grids

    # find UTM grid
    CRS=GetUTMProj(df)

    # convert lat/long to utm grid and put coordinates in dataframe
    xy<-sf::sf_project(sf::st_crs(4326), CRS, df[c("Longitude","Latitude")])
    df$X<-xy[,1]
    df$Y<-xy[,2]

    # round to the nearest 2km
    df$X2k<-2000 * floor(df$X/2000)
    df$Y2k<-2000 * floor(df$Y/2000)

    # work on unique 2km grids
    d2k <- distinct(df, X2k, Y2k)
    d2k$Z <- 1

    # create a raster grid from these points
    raster2k<-raster::rasterFromXYZ(d2k, res=c(2000,2000), crs=CRS)
    # convert to a polygon
    r2k.p<-raster::rasterToPolygons(raster2k, dissolve=T)
    # put a very small buffer around polygon to connect corners of adjacent blocks
    r2k.p.b<-raster::buffer(r2k.p, width=0.001, dissolve=T)
    # separate discontinuous areas
    r2k.d<-raster::disaggregate(r2k.p.b)

    # number of contiguous blocks of tetrads = the number of features in the polygon
    return(length(r2k.d))

}


################################################################
# internal function to guess UTM projection (based on package red)
################################################################
GetUTMProj <- function(coords){
    UTMZone = floor((min(coords$Longitude) + 180) / 6) + 1
    return(paste("+proj=utm +zone=",UTMZone," ellps=WGS84",sep=""))
}


################################################################
# main functions
################################################################

#' Evaluate criterion A2c 
#' @description Evaluate criterion A2c based on AOO & EOO trends derived from point observations
#' @param df dataframe containing columns Latitude, Longitude, Year (of observation) 
#' @param StartYear observation year to start the assessment 
#' @param EndYear observation year to end the assessment (defaults to most recent year in data)
#' @param AssessmentYears number of years over which to perform the assessment (default 10)
#' @details EvaluateA2c() employs a simple regression analysis of point distribution data
#' calculating area of occupancy and extent of occurrence values
#' to perform an assessment of distribution trend following Red List criterion A2c
#' @return A list containing the following elements
#'    Data: a data frame with columns Year, AOO, EOO
#'    AOO: a list with the following elements
#'         A2c: a2c evaluation based on AOO
#'         Area.AllTime: AOO based on all observations
#'         Area.AllAssessment: AOO over the assessment period
#'         Area.AssessmentStart: AOO at the first assessment year based on regression model
#'         Area.AssessmentEnd: AOO at the final assessment year based on regression model
#'         Slope: Slope of the regression Year~AOO (category decline per year)
#'         p: significance value of the regression Year~AOO 
#'         Change.Percent: Percentage area change from the start to the end of the assessment period
#'    EOO: a list with the following elements
#'         A2c: a2c evaluation based on EOO
#'         Area.AllTime: EOO based on all observations
#'         Area.AllAssessment: EOO over the assessment period
#'         Area.AssessmentStart: EOO at the first assessment year based on regression model
#'         Area.AssessmentEnd: EOO at the final assessment year based on regression model
#'         Slope: Slope of the regression Year~EOO (category decline per year)
#'         p: significance value of the regression Year~EOO 
#'         Change.Percent: Percentage area change from the start to the end of the assessment period
#' @examples
#' data(Alaria)
#' a2c <- EvaluateA2c(Alaria)
#' @export
EvaluateA2c<-function(df, StartYear=NA, EndYear=NA, AssessmentYears=10){

    A2c.Threshold<-c(-Inf,-0.8,-0.5,-0.3,-0.25,Inf)
    A2c.Category<-c("CR","EN","VU","NT","LC")
    A2c.Category.Threatened<-c("CR","EN","VU")

    # check if Start & End year are given
    if(is.na(StartYear) & is.na(EndYear)){
        # work out latest year and subtract Generation Time
        EndYear <- max(df$Year)
        StartYear <- EndYear-AssessmentYears+1
    } else if(is.na(StartYear)){
        # End given work out start
        StartYear <- EndYear-AssessmentYears+1
    } else if(is.na(EndYear)) {
        # start given work out end
        EndYear <- StartYear+AssessmentYears-1
    }

    # data frame of just years
    ydf=data.frame(Year=StartYear:EndYear)

    # create full output data frame
    odf<-data.frame(Year=StartYear:EndYear, AOO=NA, EOO=NA, N=NA)

    # loop through years
    for(i in 1:nrow(odf)){

        # create matrix for area calculations for coordinates in given year
        mi<-as.matrix(df[c("Latitude","Longitude")][df$Year==odf$Year[i],])

        # check we have some data
        if(!nrow(mi)==0){

            # filter to observations in given year
            odf$AOO[i]<-red::aoo(mi)
            odf$EOO[i]<-red::eoo(mi)
            odf$N[i]<-nrow(mi)
        } else {
            odf$N[i]<-0
        }

    }

    # find linear regression of area ~ year
    AOO.lm<-lm(AOO~Year, data=odf)
    EOO.lm<-lm(EOO~Year, data=odf)

    AOO.p<-coefficients(summary(AOO.lm))[2,][4]
    EOO.p<-coefficients(summary(EOO.lm))[2,][4]

    ydf=data.frame(Year=c(StartYear,EndYear))

    # predict change % over years of the evaluation
    AOO.prd<-predict(AOO.lm, newdata=ydf)
    AOO.Change.Pct<-max(AOO.prd[2],0)/max(AOO.prd[1],0)-1

    EOO.prd<-predict(EOO.lm, newdata=ydf)
    EOO.Change.Pct<-max(EOO.prd[2],0)/max(EOO.prd[1],0)-1

    # calculate area for all time, all assessment period, assessment start, assessment end
    mi<-as.matrix(df[c("Latitude","Longitude")])
    AOO.AllTime=aoo(mi)
    EOO.AllTime=eoo(mi)
    mi<-as.matrix(df[c("Latitude","Longitude")][df$Year>=StartYear & df$Year<=EndYear,])
    AOO.AllAssessment=aoo(mi)
    EOO.AllAssessment=eoo(mi)

    ## AOO
    # finally work out category based on thresholds
    A2c.AOO<-as.character(cut(AOO.Change.Pct,A2c.Threshold,A2c.Category))
    # set to NT if slope is not significant
    if(A2c.AOO %in% A2c.Category.Threatened & AOO.p>=0.05){A2c.AOO<-"NT"}
    # set to DD if no p-value
    if(is.na(AOO.p)){A2c.AOO<-"DD"}

    ## EOO
    # finally work out category based on thresholds
    A2c.EOO<-as.character(cut(EOO.Change.Pct,A2c.Threshold,A2c.Category))
    # set to NT if slope is not significant
    if(A2c.EOO %in% A2c.Category.Threatened & EOO.p>=0.05){A2c.EOO<-"NT"}
    # set to DD if no p-value
    if(is.na(EOO.p)){A2c.EOO<-"DD"}

    return(list(Data=odf,
                AOO=list(A2c=as.character(A2c.AOO),
                         Area.AllTime=AOO.AllTime,
                         Area.AllAssessment=AOO.AllAssessment,
                         Area.AssessmentStart=AOO.prd[1],
                         Area.AssessmentEnd=AOO.prd[length(AOO.prd)],
                         Slope=AOO.lm$coefficients[2],
                         p=AOO.p,
                         Change.Percent=AOO.Change.Pct*100),
                EOO=list(A2c=as.character(A2c.EOO),
                         Area.AllTime=EOO.AllTime,
                         Area.AllAssessment=EOO.AllAssessment,
                         Area.AssessmentStart=EOO.prd[1],
                         Area.AssessmentEnd=EOO.prd[length(EOO.prd)],
                         Slope=EOO.lm$coefficients[2],
                         p=EOO.p,
                         Change.Percent=EOO.Change.Pct*100)))

}


#' Evaluate criteria B1 & B2
#' @description Evaluate criteria B1 & B2 based on AOO & EOO derived from point observations
#' @param df dataframe containing columns Latitude, Longitude, Year (of observation) 
#' @param StartYear observation year to start the assessment 
#' @param EndYear observation year to end the assessment (defaults to most recent year in data)
#' @param AssessmentYears number of years over which to perform the assessment (default 10)
#' @param MinObs minimum number of observations requires to perform an assessment (default 3)
#' @details EvaluateB() calculate AOO & EOO from point distribution data
#' compare to B1 & B2 criteria thresholds and assign risk category
#' assess long term decline by comparing 
#' @return A list containing the following elements
#'    AOO: AOO area based on the assessment period
#'    AOOAllTime: AOO area over all time
#'    AOO.Decline: Boolean indicating whether there is a long term AOO decline 
#'    B2.AOO: B2 threat category derived solely from AOO calculation
#'    B2: final B2 threat category 
#'    EOO: EOO area based on the assessment period
#'    EOOAllTime: EOO area over all time
#'    B1.EOO: B1 threat category derived solely from EOO calculation
#'    B1: final B1 threat category 
#'    NTetrads: Number of tetrads (2km x 2km blocks) with observations
#'    NTetradBlocks: Number of contiguous areas based on tetrad presence
#'    B.Locations: Number of locations based on the Tetrad proxy
#'    StartYear: start year of evaluation
#'    EndYear: end year of evaluation
#' @examples
#' data(Alaria)
#' b <- EvaluateB(Alaria)
#' @export
EvaluateB<-function(df, StartYear=NA, EndYear=NA, AssessmentYears=10, MinObs=3){

    B1.Threshold<-c(0,100,5000,20000,30000,10^10) # AOO
    B2.Threshold<-c(0,10,500,2000,3000,10^10) # EOO
    B.Category<-c("CR","EN","VU","NT","LC")

    # check if Start & End year are given
    if(is.na(StartYear) & is.na(EndYear)){
        # work out latest year and subtract Generation Time
        EndYear <- max(df$Year)
        StartYear <- EndYear-AssessmentYears+1
    } else if(is.na(StartYear)){
        # End given work out start
        StartYear <- EndYear-AssessmentYears+1
    } else if(is.na(EndYear)) {
        # start given work out end
        EndYear <- StartYear+AssessmentYears-1
    }

    # convert coordinates to a matrix for area calculations
    mi<-as.matrix(df[c("Latitude","Longitude")])
    AOO.AllTime<-red::aoo(mi)
    EOO.AllTime<-red::eoo(mi)

    # convert coordinates to a matrix for area calculations
    mi<-as.matrix(df[c("Latitude","Longitude")][df$Year>=StartYear & df$Year<=EndYear,])
    AOO<-red::aoo(mi)
    EOO<-red::eoo(mi)
    
    B2.AOO<-cut(AOO, breaks=B2.Threshold, labels=B.Category)
    B1.EOO<-cut(EOO, breaks=B1.Threshold, labels=B.Category)

    # other criteria required for B1/B2 - evidence of continued decline 
    # define this as present AOO less than 1/2 of all time AOO
    AOO.Decline <- AOO < AOO.AllTime/2

    # other criteria required for B1/B2 - number of locations
    # this isn't directly measurable - but have proxies - 
    # proxy 1) number of 2km grids aka tetrads
    NTetrads <- AOO/4

    # proxy 2) number of contiguous tetrad blocks
    NTetradBlocks <- EvaluateTetradBlocks(dplyr::filter(df, Year>=StartYear & Year<=EndYear))

    # use the smallest of the two proxies
    NLocationProxy <- min(NTetrads, NTetradBlocks)

    # find category based on number of locations (proxy)
    BLocThreshold <- c(0,1,5,10,10.000001,Inf)
    BLocCat <- cut(NLocationProxy, breaks=BLocThreshold, labels=B.Category)

    # now combine area with decline and number of location proxy
    if(nrow(mi)<MinObs){
        B2="DD"
        B1="DD"
    } else if(!AOO.Decline) {
        B1="LC"
        B2="LC"
    } else {
        B2 = B.Category[max(as.integer(B2.AOO), as.integer(BLocCat))]
        B1 = B.Category[max(as.integer(B1.EOO), as.integer(BLocCat))]
    }

    return(list(AOO=AOO, AOOAllTime=AOO.AllTime, AOO.Decline=AOO.Decline,
                B2.AOO=as.character(B2.AOO), B2=B2,
                EOO=EOO, EOOAllTime=EOO.AllTime, B1.EOO=as.character(B1.EOO), B1=B1,
                NTetrads=NTetrads, NTetradBlocks=NTetradBlocks, B.Locations=BLocCat,
                StartYear=StartYear, EndYear=EndYear))
}


#' Evaluate criterion D2
#' @description Evaluate criteria D2 based on AOO
#' @param df dataframe containing columns Latitude, Longitude, Year (of observation) 
#' @param StartYear observation year to start the assessment 
#' @param EndYear observation year to end the assessment (defaults to most recent year in data)
#' @param MinObs minimum number of observations requires to perform an assessment (default 1)
#' @details EvaluateD2() calculate AOO from point distribution data and compares to area threshold D2
#' Warning! This makes no evaluation of ongoing threats. To formally apply the D2 criterion you must document a threat
#' @return A list with the following elements
#'   AOO: AOO area
#'   D2: D2 threat category based solely on AOO area (no consideration of ongoing threats)
#' @examples
#' data(Alaria)
#' b <- EvaluateD2(Alaria)
#' @export
EvaluateD2<-function(df, StartYear=NA, EndYear=NA, MinObs=1){

    D2.Threshold<-c(0,20,Inf)
    D2.Category<-c("VU", "LC")

    if(is.na(StartYear)){StartYear<-min(df$Year)}
    if(is.na(EndYear)){EndYear<-max(df$Year)}

    # convert coordinates to a matrix for area calculations
    mi<-as.matrix(df[c("Latitude","Longitude")][df$Year>=StartYear & df$Year<=EndYear,])

    if(nrow(mi)>=MinObs){
        # find AOO
        AOO<-red::aoo(mi)
        D2<-cut(AOO, breaks=D2.Threshold, labels=D2.Category)
    } else {
        AOO <- NA
        D2 <- "DD"
    }

    return(list(AOO=AOO, D2=as.character(D2)))
}


#' Evaluate criterion A2b
#' @description Evaluate criteria A2b based on SACFOR abundance indices
#' @param df dataframe containing columns Latitude, Longitude, Grid (grid reference), Year (of observation) and Abundance (SACFOR category)
#' @param StartYear observation year to start the assessment 
#' @param EndYear observation year to end the assessment (defaults to most recent year in data)
#' @param AssessmentYears number of years over which to perform the assessment (default 10)
#' @param MinGrids minimum number of repeat observations required to perform an assessment (default 5)
#' @details EvaluateA2b() Measure change in abudance based on SACFORN index
#' compare change over the assessment period with A2b criterion thresholds
#' @return A list with the following elements
#'   Data: A data frame with columns GridRef, Year & Abundance (SACFOR category)
#'   Trend: Observed trend in abundance (category change per year)
#'   Trend.Percentiles: Percentils of observed trends
#'   Percent: Percentage change in abundance over the assessment period based on the trend
#'   NGrids: Number of grid cells with repeat surveys used for assessment
#'   NObs: Number of observations used in the assessment (includes multiple observations at each grid cell)
#'   A2b: threat category based on A2b criterion thresholds
#' @examples
#' data(Alaria)
#' a2b <- EvaluateA2b(Alaria)
#' @export
EvaluateA2b<-function(df, StartYear=NA, EndYear=NA, AssessmentYears=10, MinGrids=5){

    SACFOR<-c("N", "R", "O", "F", "C", "A", "S", "SA")
    A2b.Threshold<-c(-Inf, -0.8, -0.5, -0.3, 0, +Inf)
    A2b.Category<-c("CR","EN","VU","NT","LC")

    # work out start and end points
    if(is.na(StartYear) & is.na(EndYear)){
        # work out latest year and subtract Generation Time
        EndYear <- max(df$Year)
        StartYear <- EndYear-AssessmentYears+1
    } else if(is.na(StartYear)){
        # End given work out start
        StartYear <- EndYear-AssessmentYears+1
    } else if(is.na(EndYear)) {
        # start given work out end
        EndYear <- StartYear+AssessmentYears-1
    }

    # only consider years within the start & end range
    df <- dplyr::filter(df, Year>=StartYear & Year<=EndYear & !is.na(Abundance) & Abundance %in% SACFOR)

    if(nrow(df)>0){

        # convert SACFOR abundance to integer values
        df$SACFOR<-match(as.character(df$Abundance), SACFOR)

        # start by filtering to unique combinations of grid & year
        df.grid.year <- df %>% dplyr::group_by(Grid,Year) %>% dplyr::summarise(AbMean=mean(SACFOR))

        # filter for only grids with more than one year of observations
        repeat.grids<- df.grid.year %>% dplyr::group_by(Grid) %>% dplyr::summarise(n=n()) %>% dplyr::filter(n>1)

        repeat.grids$Trend<-NA

        # loop through grids with multiple years of data
        for(i in 1:nrow(repeat.grids)){
            # subset for grid cell
            df.g <- dplyr::filter(df.grid.year, Grid==repeat.grids$Grid[i])

            df.g.trend<-lm(AbMean~Year, data=df.g)

            repeat.grids$Trend[i]<-df.g.trend$coefficients[[2]]
        }

        # find median trend
        Abundance.Trend <- median(repeat.grids$Trend)
        # convert to % change
        Abundance.Pct.Change<-sign(Abundance.Trend)*(1-(0.5 ^ abs(Abundance.Trend)))
        A2b<-A2b.Category[cut(Abundance.Pct.Change, A2b.Threshold, right = FALSE)]
        # set to DD if no trend or too few observations
        if(is.na(A2b) | nrow(repeat.grids)<MinGrids){A2b<-"DD"}


        # return data for plotting
        df.plot <- dplyr::filter(df.grid.year, Grid %in% repeat.grids$Grid)

        Trend.Percentiles <-quantile(repeat.grids$Trend,
                                     c(0.05, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95))


        return(list(Data=df.plot, Trend=Abundance.Trend,
                    Trend.Percentiles=Trend.Percentiles,
                    Percent=Abundance.Pct.Change*100,
                    NGrids=nrow(repeat.grids),
                    NObs=sum(repeat.grids$n),
                    A2b=as.character(A2b),
                    SACFOR=SACFOR))
    } else {
        return(list(Data=NA, Trend=NA,
                    Trend.Percentiles=NA,
                    Percent=NA,
                    NGrids=0,
                    NObs=0,
                    A2b="DD"))
    }
}

#' Evaluate all Red List Criteria
#' @description Evaluate all Red List criter from point observation data
#' @param df dataframe containing columns Latitude, Longitude, Year (of observation) 
#' @param StartYear observation year to start the assessment 
#' @param EndYear observation year to end the assessment (defaults to most recent year in data)
#' @param AssessmentYears number of years over which to perform the assessment (default 10)
#' @param B.MinObs minimum number of observations requires to perform criterion B evaluation (default 3)
#' @param D2.MinObs minimum number of observations requires to perform criterion D2 evaluation (default 1)
#' @param A2b.MinGrids minimum number of repeat observations required to perform criterion A2b evaluation (default 5)
#' @details EvaluateAllCriteria() Evaluate all Red List criteria by calling all criterion evaluations
#' @return A list with the following elements
#'   A2b: A list returned from function EvaluateA2b()
#'   A2c: A list returned from function EvaluateA2c()
#'   B: A list returned from function EvaluateB()
#'   D2: A list returned from function EvaluateD2()
#'   Final.Evaluation: The final evaluation based on all criteria
#' @examples
#' data(Alaria)
#' all <- EvaluateAllCriteria(Alaria)
#' @export
EvaluateAllCriteria<-function(df, StartYear=NA, EndYear=NA, AssessmentYears=10,
                         A2b.MinGrids=5, D2.MinObs=1, B.MinObs=3){

    # find all criteria and return as a list
    A2b<-EvaluateA2b(df, StartYear=StartYear, EndYear=EndYear, AssessmentYears=AssessmentYears,
                MinGrids=A2b.MinGrids)
    A2c<-EvaluateA2c(df, StartYear=StartYear, EndYear=EndYear, AssessmentYears=AssessmentYears)
    B<-EvaluateB(df, StartYear=StartYear, EndYear=EndYear,
              MinObs=B.MinObs)
    D2<-EvaluateD2(df, StartYear=StartYear, EndYear=EndYear, 
              MinObs=D2.MinObs)

    AllCriteria<-list(A2b=A2b, A2c=A2c, B=B, D2=D2)

    Final<-GetFinalEvaluation(AllCriteria)

    return(list(A2b=A2b, A2c=A2c, B=B, D2=D2, Final.Evaluation=Final))
}


#' Plot data used for criterion A2c evaluation
#' @description Plot data used for criterion A2c evaluation
#' @param A2c a list output from EvaluateA2c()
#' @param type a character string either "AOO" or "EOO" indicationg which metric to plot
#' @details Creates a scatter plot of Area~Year used for A2c evaluation
#' @return A ggplot object
#' @examples
#' data(Alaria)
#' A2c <- EvaluateA2c(Alaria)
#' PlotA2c(A2c)
#' @export
PlotA2c<-function(A2c, type="AOO"){

    if(type=="AOO"){

       # plot A2c data from AOO
        p<-ggplot(A2c$Data, aes(Year, AOO)) +
            geom_point() +
            geom_smooth(method="lm") +
            ggtitle(paste("A2c Criteria (AOO) =", A2c$AOO$A2c))

    } else if(type=="EOO"){
       # plot A2c data from AOO
        p<-ggplot(A2c$Data, aes(Year, EOO)) +
            geom_point() +
            geom_smooth(method="lm") +
            ggtitle(paste("A2c Criteria (EOO) =", A2c$EOO$A2c))
    } else {
        print(paste("Don't recognise this category: " , type, "must be AOO or EOO"))
        p<-NA
    }


    return(p)
    
}


#' Plot data used for criterion A2b evaluation
#' @description Plot data used for criterion A2b evaluation
#' @param A2b a list output from EvaluateA2b()
#' @details Creates a scatter plot of Abundance~Year used for A2b evaluation
#' @return A ggplot object
#' @examples
#' data(Alaria)
#' A2b <- EvaluateA2b(Alaria)
#' PlotA2b(A2b)
#' @export
PlotA2b<-function(A2b){

    A2bd<-as.data.frame(A2b$Data)
    A2bd$Observation<-"Site"

    yrs<-min(A2b$Data$Year):max(A2b$Data$Year)

    yrs.mm<-c(min(A2b$Data$Year),max(A2b$Data$Year))
    yr.dif<-yrs.mm[2]-yrs.mm[1]+1

    mean.line<-data.frame(Grid="", Observation="Trend", 
                          AbMean=mean(A2b$Data$AbMean)+c(-1,1)*(A2b$Trend*yr.dif/2),
                          Year=c(min(A2b$Data$Year),max(A2b$Data$Year)))

    A2bd<-rbind(A2bd, mean.line)

    p<-ggplot() +
        geom_line(data=A2bd, aes(x=Year, y=AbMean, group=Grid, colour=Observation)) +
        scale_color_manual(values = c("Site" = "grey", "Trend" = "black")) +
        scale_y_continuous(name="Abundance Category", breaks=1:8, labels=A2b$SACFOR) + 
        scale_x_continuous(breaks=yrs, labels=yrs) +
        ggtitle(paste("A2b Criteria (Abundance) =", A2b$A2b))

    return(p)

}


#' Plot map based on point observations
#' @description Plot map based on point observations used in evaluation 
#' @param df a data frame used in evaluations
#' @param B a list returned from EvaluateB()
#' @param D2 a list returned from EvaluateD2()
#' @param CRS a coordinate reference system definition recognisable by sf::st_crs - will guess a local UTM zone if left blank
#' @param GridRes a value indicating resolution 
#' @details Plot a simple map based on unique tetrads of point observations used in evaluation
#' display results of B1&2 and D2 evaluations too
#' @return A ggplot object
#' @examples
#' data(Alaria)
#' B <- EvaluateB(Alaria)
#' D2 <- EvaluateD2(Alaria)
#' PlotArea(Alaria, B, D2)
#' @export
PlotArea <- function(df, B, D2, CRS=NA, GridRes=2000){

    # load world map
    world <- rnaturalearth::ne_countries(scale="medium", returnclass="sf")

    # find utm zone
    if(is.na(CRS)){CRS=GetUTMProj(df)}

    # convert lat/long to utm grid and put coordinates in dataframe
    xy<-sf::sf_project(sf::st_crs(4326), CRS, df[c("Longitude","Latitude")])
    df$X<-xy[,1]
    df$Y<-xy[,2]

    # round to the nearest 2km
    df$X2k<-GridRes * floor(df$X/GridRes)
    df$Y2k<-GridRes * floor(df$Y/GridRes)

    # Get unique grids for before, during & after the assessment period
    df$Period <- ifelse(df$Year<B$StartYear, "Before Assessment", ifelse(df$Year>B$EndYear, "After Assessment", "During Assessment"))

    # get unique points for before, during & after the assessment period
    df.p<- distinct(df, X2k, Y2k, Period)

    p<-ggplot(data = world) +
        geom_sf(fill="darkseagreen1") +
        coord_sf(xlim=c(min(df.p$X2k), max(df.p$X2k)), ylim=c(min(df.p$Y2k), max(df.p$Y2k)),
                 crs = CRS) +
        geom_point(data=df.p %>% arrange(Period), aes(x=X2k, y=Y2k, colour=factor(Period))) +
        annotation_scale(location = "bl", width_hint = 0.5) +
        labs(colour="Observation Period", x="", y="") +
        ggtitle(paste("Area Criteria ",
                      "B1 ", B$B1, " (EOO = ", B$EOO, " = ", B$B1.EOO, ", Decline = ", substr(B$AOO.Decline,1,1), ", Locations = ", B$NTetradBlocks, " = ", B$B.Locations, ")\n                     ",
                      "B2 ", B$B2, " (AOO = ", B$AOO, " = ", B$B2.AOO, ", Decline = ", substr(B$AOO.Decline,1,1), ", Locations = ", B$NTetradBlocks, " = ", B$B.Locations, ")\n                     ",
                      "D2 (AOO = ", D2$AOO, ") = ", D2$D2, sep=""))
    
    return(p)

}

################################################################
# data layers

#' Test dataset for Alaria esculenta from uk (NBN Atlas 2021).
#'
#' Data frame containing Latitude, Longitude, Year & Abundance for Alaria esculenta from UK (NBN Atlas 2021)
#'
#' @docType data
#' @keywords datasets
#' @name Alaria
#' @usage data(Alaria)
#' @format dataframe object with columns Name, Latitude, Longitude, Year & Abundance with data for Alaria esculenta from UK 
#' @references Joint Nature Conservation Committee (2023). Marine Nature Conservation Review (MNCR) and associated benthic marine data held and managed by JNCC. Occurrence dataset on the NBN Atlas.
NULL

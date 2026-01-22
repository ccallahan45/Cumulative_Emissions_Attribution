# risk of extreme events in the multi-model large ensemble

library(dplyr)
library(ncdf4)
library(tidync)
library(extRemes)
library(nanoparquet)
library(tidyr)
library(zoo)
library(stringr)
library(MASS)

select <- dplyr::select
filter <- dplyr::filter

loc_extremes <- "/N/project/climatesociety/datasets/MMLEA/annual/"
loc_emissions <- "/N/project/climatesociety/datasets/CMIP6/Emissions/"
loc_gcp <- "/N/project/climatesociety/datasets/GlobalCarbonProject/"
loc_lsm <- "/N/project/climatesociety/datasets/LandSea/2p5degree/"
loc_out_emis <- "/N/u/ccallah/Quartz/Projects/Extremes_Emissions/Data/Emissions/"
loc_out_params <- "/N/u/ccallah/Quartz/Projects/Extremes_Emissions/Data/Fitted_Parameters/"


indices <- c("txx","rx1day")
# pool all models and scenarios
# match emissions to scenario

process_emissions <- FALSE
if (process_emissions){
    # read hist emissions
    hist_emissions <- read.csv(paste0(loc_gcp,"global_carbon_budget_extracted.csv")) 
    hist_emissions_world <- hist_emissions %>% rename(year=Year) %>% 
        mutate(land.use.change.emissions = if_else(is.na(land.use.change.emissions), 0, land.use.change.emissions)) %>%
        mutate(emissions_gtc = fossil.emissions.excluding.carbonation+land.use.change.emissions) %>%
        mutate(emissions_gtco2 = emissions_gtc*3.667) %>%
        select(year,emissions_gtco2)
    #print(hist_emissions_world)
    
    # read ssp emissions
    ssp_emissions <- read.csv(paste0(loc_emissions,"SSP_CMIP6_201811.csv"))
    ssp_emissions_world <- ssp_emissions %>% filter(REGION=="World",VARIABLE=="CMIP6 Emissions|CO2") %>%
        select(-c(REGION,MODEL,VARIABLE)) %>%
        pivot_longer(!c(SCENARIO,UNIT),names_prefix="X",names_to="year",values_to="emissions") %>%
        mutate(emissions_gtco2 = emissions/1000) %>%
        mutate(scenario=tolower(SCENARIO),year=as.numeric(year))
    #print(ssp_emissions_world)
    
    ## create full array with diff scenarios and every year, merge
    scenarios <- tolower(unique(ssp_emissions_world$SCENARIO))
    years <- c(1750:2100)
    emissions <- data.frame("year"=years,"scenario"=gsub("-","",scenarios[1])) %>%
        left_join(hist_emissions_world,by="year")
    emissions[which(emissions$year%in%c(2030,2040,2050,2060,2070,2080,2090,2100)),"emissions_gtco2"] <- 
        ssp_emissions_world[(ssp_emissions_world$year>=2030)&(tolower(ssp_emissions_world$scenario)==scenarios[1]),"emissions_gtco2"]
    
    for (s in scenarios[2:9]){
        if (s=="ssp3-lowntcf"){
            next
        } else if (s=="ssp3-70 (baseline)"){
            scen = "ssp370"
        } else if (s=="ssp5-85 (baseline)"){
            scen = "ssp585"
        } else {
            scen = gsub("-","",s)
        }
        print(scen)
        emis <- data.frame("year"=years,"scenario"=scen) %>%
            left_join(hist_emissions_world,by="year")
        emis[which(emis$year%in%c(2030,2040,2050,2060,2070,2080,2090,2100)),"emissions_gtco2"] <- 
            ssp_emissions_world[(ssp_emissions_world$year>=2030)&(tolower(ssp_emissions_world$scenario)==s),"emissions_gtco2"]
        emissions <- rbind(emissions,emis)
        
    }
    emissions_interp <- emissions %>% mutate(emissions_gtco2_interp = zoo::na.approx(emissions_gtco2)) %>%
        group_by(scenario) %>%
        mutate(cumulative_emissions_gtco2 = cumsum(emissions_gtco2_interp))
    
    write.csv(emissions_interp,paste0(loc_out_emis,"global_emissions_historical_ssp_interpolated_1750-2100.csv"))
} else {
    emissions_interp <- read.csv(paste0(loc_out_emis,"global_emissions_historical_ssp_interpolated_1750-2100.csv")) %>% 
        select(year,scenario,emissions_gtco2_interp,cumulative_emissions_gtco2) %>%
        rename(emissions=emissions_gtco2_interp,cumul_emissions=cumulative_emissions_gtco2)
}


# func
nonstationary_gev_params <- function(df,varname,covar){
    mdl <- fevd(as.numeric(unlist(df[varname])), df, 
                type = "GEV",method="MLE",
                location.fun=as.formula(paste0("~",covar)), 
                    scale.fun=as.formula(paste0("~",covar)))
    mu0 <- as.numeric(mdl$results$par["mu0"])
    mu1 <- as.numeric(mdl$results$par["mu1"])
    sigma0 <- as.numeric(mdl$results$par["sigma0"])
    sigma1 <- as.numeric(mdl$results$par["sigma1"])
    shape <- as.numeric(mdl$results$par["shape"])
    df_out <- data.frame("mu0"=mu0,"mu1"=mu1,"sigma0"=sigma0,
                            "sigma1"=sigma1,"shape"=shape)
    return(df_out)
}

nonstationary_gev_loc_params <- function(df,varname,covar){
    mdl <- fevd(as.numeric(unlist(df[varname])), df, 
                type = "GEV",method="MLE",
                location.fun=as.formula(paste0("~",covar)))
    mu0 <- as.numeric(mdl$results$par["mu0"])
    mu1 <- as.numeric(mdl$results$par["mu1"])
    scale <- as.numeric(mdl$results$par["scale"])
    shape <- as.numeric(mdl$results$par["shape"])
    df_out <- data.frame("mu0"=mu0,"mu1"=mu1,"scale"=scale,"shape"=shape)
    return(df_out)
}



## read land sea mask
lsm = tidync::hyper_tibble(paste0(loc_lsm,"IMERG_land1_sea0_mask_2p5degree.nc"),na.rm=F) %>%
    mutate(lat=as.numeric(lat),lon=as.numeric(lon)) %>% rename(land=landseamask) %>%
    select(lat,lon,land)

## first fit obs so we can compare to models
extremes <- c("txx","rx1day") 
set.seed(110224)
for (e in extremes){
    next
    if (e=="txx"){obsdata<-"era5"}else{obsdata<-"cpc"}
    loc_obs <- paste0("/N/project/climatesociety/datasets/",toupper(obsdata),"/annual/")
    y1_obs <- 1979
    y2_obs <- 2024

    # read data
    for (yy in c(y1_obs:y2_obs)){
        obs_in <- tidync::hyper_tibble(paste0(loc_obs,obsdata,".",e,".2p5.",yy,".nc"),na.rm=FALSE) %>%
            mutate(lon=as.numeric(lon),lat=as.numeric(lat),year=yy) %>% select(lat,lon,year,!!as.name(e))
        if (yy==y1_obs){obs<-obs_in}else{obs<-rbind(obs,obs_in)}
    }
    
    obs_emissions = obs %>% left_join(emissions_interp[emissions_interp$scenario=="ssp245",],by="year") %>%
            left_join(lsm,by=c("lat","lon"))
    
    # fit params
    params_obs <- obs_emissions %>% filter(land==1, lat > -65, lat < 65) %>%
            group_by(lat,lon) %>% 
            filter(!any(is.na(!!as.name(e)))) %>% # filter out grid points with any na 
            summarize(params = nonstationary_gev_params(pick(everything()),e,"cumul_emissions"),.groups="keep") %>%
            unpack(cols = params)

    final_df <- obs_emissions %>% group_by(lat,lon) %>% summarize(n=n(),.groups = "keep") %>%
            left_join(params_obs,by=c("lat","lon")) %>% select(-n)
    #print(final_df)
    
    fname <- paste0(loc_out_params,obsdata,"_",e,"_gev_cumulativeemissions_params.parquet")
    write_parquet(final_df,fname)
    print(final_df)
    print(fname)
}


## we'll do the fits to the models twice
# the main analysis over a longer period
# then a shorter period to compare to the obs
y1_fit = 1850
y2_fit = 2100
y1_fit_obs = 1979
y2_fit_obs = 2024


## now loop through models and assemble annual extremes across them, combine with emissions
extremes <- c("txx","rx1day") 
scen = "ssp585" # "ssp585"

for (e in extremes){
    print(e)
    
    # models
    loc_in <- paste0(loc_extremes,e,"/biascorrected/")
    file_list_all <- list.files(loc_in)
    file_list <- file_list_all[str_detect(file_list_all,"1850-2100")]
    file_list_s <- file_list[str_detect(file_list,scen)]

    print(length(file_list_s))
    
    for (ff in c(1:length(file_list_s))){
        mdl_fname <- file_list_s[ff]
        print(mdl_fname)
        
        strings <- str_split(mdl_fname,"_")
        model <- paste0(strings[[1]][1],"_",strings[[1]][2])
        print(model)
        
        
        model_in <- tidync(paste0(loc_in,mdl_fname))

        # the extreme event index
        model_ext <- tidync::hyper_tibble(model_in,na.rm=F) %>%
                mutate(year=as.numeric(time),lat=as.numeric(lat),lon=as.numeric(lon)) %>%
                select(year,lat,lon,!!as.name(e)) %>%
                mutate(model=model,scenario=scen)
        
        # pvals from k-s test against obs
        mdl_pvals <- model_in %>% tidync::activate("ks_pvalue") %>% hyper_tibble() %>%
                mutate(lat=as.numeric(lat),lon=as.numeric(lon))

        # join with emissions and land mask
        ptm <- as.numeric(proc.time()[3])
        ext_m <- model_ext %>% filter(year>=y1_fit,year<=y2_fit) %>%
            left_join(emissions_interp,by=c("year","scenario")) %>%
            left_join(lsm,by=c("lat","lon"))
        
        params_m <- ext_m %>% filter(land==1, lat > -65, lat < 65) %>%
            group_by(lat,lon) %>% 
            filter(!any(is.na(!!as.name(e)))) %>% # filter out grid points with any na 
            summarize(params = nonstationary_gev_params(pick(everything()),e,"cumul_emissions"),.groups="keep") %>%
            #summarize(params = nonstationary_gev_loc_params(pick(everything()),paste0(e,"_bc"),"cumul_emissions"),.groups="keep") %>%
            unpack(cols = params)
        
        print(paste0("elapsed time: ",round(as.numeric(proc.time()[3] - ptm)/60,digits=1)," minutes"))

        final_df_n <- ext_m %>% group_by(lat,lon) %>% summarize(n=n(),.groups = "keep") %>%
            left_join(params_m,by=c("lat","lon")) %>% select(-n) %>% mutate(model=model,scenario=scen) %>%
            left_join(mdl_pvals,by=c("lat","lon"))
        
        if (ff==1){final_df<-final_df_n}else{final_df<-rbind(final_df,final_df_n)}
    }

    ## write out final df as parquet
    fname <- paste0(loc_out_params,"mmlea_",e,"_gev_",scen,"_cumulativeemissions_bymodel_params_",y1_fit,"-",y2_fit,".parquet")
    write_parquet(final_df,fname)
    print(final_df)
    print(fname)
}


## lastly, do the same fits over the observational interval of 1979-2024
for (e in extremes){
    next
    print(e)

    # models
    loc_in <- paste0(loc_extremes,e,"/biascorrected/")
    file_list_all <- list.files(loc_in)
    file_list <- file_list_all[str_detect(file_list_all,"1850-2100")]
    file_list_s <- file_list[str_detect(file_list,scen)]

    print(length(file_list_s))
    
    for (ff in c(1:length(file_list_s))){
        mdl_fname <- file_list_s[ff]
        strings <- str_split(mdl_fname,"_")
        model <- paste0(strings[[1]][1],"_",strings[[1]][2])
        print(model)
        #if (model != "MIROC-ES2L_r10i1p1f2"){next}
        if (strings[[1]][1]=="MIROC-ES2L"){
            next
        }
        
        model_in <- tidync(paste0(loc_in,mdl_fname))

        # the extreme event index
        model_ext <- tidync::hyper_tibble(model_in,na.rm=F) %>%
                mutate(year=as.numeric(time),lat=as.numeric(lat),lon=as.numeric(lon)) %>%
                select(year,lat,lon,!!as.name(e)) %>%
                mutate(model=model,scenario=scen)
        
        # pvals from k-s test against obs
        mdl_pvals <- model_in %>% tidync::activate("ks_pvalue") %>% hyper_tibble() %>%
                mutate(lat=as.numeric(lat),lon=as.numeric(lon))

        # join with emissions and land mask
        ptm <- as.numeric(proc.time()[3])
        ext_m <- model_ext %>% filter(year>=y1_fit_obs,year<=y2_fit_obs) %>%
            left_join(emissions_interp,by=c("year","scenario")) %>%
            left_join(lsm,by=c("lat","lon"))
        
        params_m <- ext_m %>% filter(land==1, lat > -65, lat < 65) %>%
            group_by(lat,lon) %>% 
            filter(!any(is.na(!!as.name(e)))) %>% # filter out grid points with any na 
            summarize(params = nonstationary_gev_params(pick(everything()),e,"cumul_emissions"),.groups="keep") %>%
            #summarize(params = nonstationary_gev_loc_params(pick(everything()),paste0(e,"_bc"),"cumul_emissions"),.groups="keep") %>%
            unpack(cols = params)

        print(paste0("elapsed time: ",round(as.numeric(proc.time()[3] - ptm)/60,digits=1)," minutes"))

        final_df_n <- ext_m %>% group_by(lat,lon) %>% summarize(n=n(),.groups = "keep") %>%
            left_join(params_m,by=c("lat","lon")) %>% select(-n) %>% mutate(model=model,scenario=scen) %>%
            left_join(mdl_pvals,by=c("lat","lon"))
        
        if (ff==1){final_df<-final_df_n}else{final_df<-rbind(final_df,final_df_n)}
    }

    ## write out final df as parquet
    fname <- paste0(loc_out_params,"mmlea_",e,"_gev_",scen,"_cumulativeemissions_bymodel_params_",y1_fit_obs,"-",y2_fit_obs,".parquet")
    write_parquet(final_df,fname)
    print(final_df)
    print(fname)
}


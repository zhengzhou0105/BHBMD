## Functions for plotting BHBMD results
## Author: Zheng Zhou
## Date: Oct 24 2022

library("extraDistr")

## Arranged by dose response models to which the models were fit

## Linear-----
plot_specific_Linear <- function(stanfit,DataList,
                                 aname = "a",bname = "b",
                                 xup = NULL, yup = NULL,
                                 savetag = "Specific"){
  # stanfit posterior
  df_posteriors <- as.data.frame(stanfit)
  iterEff <- nrow(df_posteriors)
  NStudy <- length(unique(DataList$Index))          # number of studies
  G <- as.integer(table(DataList$Index)) # number of dose group in each study
  
  # curves as xy plot
  doselow <- 0
  dosehigh <- ceiling(max(DataList$dose))
  # doselow and dosehigh are theoretical range of doses to be as wide as possible
  xdose <- seq(from = doselow, to = dosehigh,
               by = round((dosehigh - doselow)/1e3,3) )
  # convert xaxis coordinates because parameters were fitted on standarized dose
  mindose <- min(DataList$dose)
  maxdose <- max(DataList$dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose
  doseinput <- dosed
  # doseinput <- xdose
  
  # extract posteriors and plot by study
  y_specific <- matrix(NA, nrow = iterEff,ncol = length(doseinput))
  
  for(s in 1:NStudy){
    
    # calculate study specific y
    a_specific <- df_posteriors[,paste0(aname,"[",s,"]")]
    b_specific <- df_posteriors[,paste0(bname,"[",s,"]")]
    
    # Two ways of plotting
    
    # Method 1: dose ~ summary of iterative RR by all parameters
    # each iter of pars fits a distribution of y then take descriptives
    # y_specific is the dataframe of y based on study-specific parameters with simulation
    for(i in 1:iterEff){
      # use standardized dose input
      y_specific[i,] <- get_Linear(doseinput, a = a_specific[i], b = b_specific[i])
    }

    # summary data
    df_plot <- data.frame(x = xdose,
                          median = apply(y_specific,MARGIN =  2,
                                         FUN= function(x) median(x,
                                                                 na.rm = T)),
                          mean = apply(y_specific,MARGIN =  2,
                                       FUN= function(x) mean(x,
                                                             na.rm = T)),
                          Q5 = apply(y_specific,MARGIN =  2,
                                     FUN= function(x) quantile(x,0.05,
                                                               na.rm = T)),
                          Q95 = apply(y_specific,MARGIN =  2,
                                      FUN= function(x) quantile(x,0.95,
                                                                na.rm = T)),
                          Q25 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.25,
                                                                 na.rm = T)),
                          Q75 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.75,
                                                                 na.rm = T))
    )
    
    # # Method 2: dose ~ RR by median parameters
    # y_median <- get_Linear(doseinput, a = median(a_specific),b = median(b_specific))
    # y_p5 <- get_Linear(doseinput, a = quantile(a_specific,0.05),b = quantile(b_specific,0.05))
    # y_p95 <- get_Linear(doseinput, a = quantile(a_specific,0.95),b = quantile(b_specific,0.95))
    # y_p25 <- get_Linear(doseinput, a = quantile(a_specific,0.25),b = quantile(b_specific,0.25))
    # y_p75 <- get_Linear(doseinput, a = quantile(a_specific,0.75),b = quantile(b_specific,0.75))
    # df_plot <- data.frame(
    #   x = xdose,
    #   median = y_median, Q5 = y_p5, Q95 = y_p95,Q25=y_p25,Q75=y_p75
    # )
    
    # extract study specific raw data
    Reference <- DataList$Study[(sum(G[1:s-1])+1)]
    Dose <-  DataList$dose[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR <- DataList$ymean[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l <- DataList$RR.l[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.u <- DataList$RR.u[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l[is.na(RR.l)] <- RR[is.na(RR.l)]
    RR.u[is.na(RR.u)] <- RR[is.na(RR.l)]
    
    # data points of each study
    df_raw_study <- data.frame(
      Dose = Dose,
      RR.l = RR.l,
      RR.u = RR.u,
      RR = RR
    )

    ## plotting parameters
    xupper <- ifelse(is.null(xup),   # if a upper bound not specified
                     ceiling (max(
                       median(Dose)*1.1,
                       median(Dose) + 1,
                       max(Dose)
                     )),
                     xup
    )
    
    yupper <- ifelse(is.null(yup),
                     ceiling(
                       min(
                         max(RR.u,na.rm = T)*1.1,
                         max(RR.u,na.rm = T) +1
                       )
                     ),
                     yup)
    
    # get a blank plot base
    graphics.off()
    myplotbase <- ggplot(data.frame())+ geom_blank()
    
    plot_specific_summary <- myplotbase+
      # title to display reference
      labs(title = paste0("Study #",s," ",Reference),
           x = "ADD(ug/kg/D)",
           y = "RR")+
      # median curve
      geom_line(data = df_plot,
                aes(x = x , y = median),
                color = "blue", size = 1)+
      # # # mean curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = mean),
      #           color = "red", size = 2)+
      # 5% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q5),
                color = "green", size = 1,linetype ="dotted")+
      # 95% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q95),
                color = "red", size = 1,linetype ="dotted")+
      # # 25% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q25),
      #           color = "green", size = 0.8,linetype ="dotted")+
      # # 75% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q75),
      #           color = "red", size = 0.8,linetype ="dotted")+
      # data points median
      geom_point(data = df_raw_study,
                 aes(x = Dose,y = RR),
                 size = 5, color = s)+
      # upper and lower bound of RR
      geom_segment(data = df_raw_study,
                   aes(x = Dose,xend = Dose,
                       y = RR.l,yend = RR.u),
                   size = 1, color= s)+
      # set axis coordinates bound
      coord_cartesian(xlim = c(0,xupper),
                      ylim = c(0,yupper))+
      scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
      scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))+
      theme_classic(base_size = 20)
    
    # create figs
    # dev.new("windows",noRStudioGD = T)
    # save as svg
    svg(filename = paste0(
      "Linear_", savetag," Study #",s," ",Reference,".svg"
      ),
        width = 15,height = 15)
    print(
      plot_specific_summary
    )
    dev.off()
  }
  
}

plot_overarching_Linear <- function(stanfit,DataList,
                                    aname = "a",bname = "b",
                                    xup = NULL, yup = NULL,
                                    savetag = "curve"){
  # posteriors
  df_posteriors <- as.data.frame(stanfit)
  iterEff <- nrow(df_posteriors)
  NStudy <- length(unique(DataList$Index))          # number of studies
  
  # Calculate overarching parameters
  a <- rlnorm(iterEff,df_posteriors$mu_a,df_posteriors$sigma_a)
  b <- extraDistr::rhcauchy(iterEff,df_posteriors$scale_b)
  
  # raw data
  dose <- DataList$dose
  doselow <- 0
  dosehigh <- ceiling(max(dose))
  xdose <- seq(from = doselow, to = dosehigh,
               by = round((dosehigh - doselow)/1e3,3))
  # if lower or upper bound of RR is missing, use RR values
  DataList$RR.l[is.na(DataList$RR.l)] <- DataList$RR[is.na(DataList$RR.l)]     
  DataList$RR.u[is.na(DataList$RR.u)] <- DataList$RR[is.na(DataList$RR.l)] 
  
  # convert xaxis coordinates
  mindose <- min(dose)
  maxdose <- max(dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose as input for y calculation
  doseinput <- dosed
  # doseinput <- xdose
  
  # plotting method 1: dose ~ median of iterative RR by all parameters
  y_overarching <- matrix(NA, nrow = iterEff,ncol = length(doseinput))

  # each iter of pars fits a distribution of y then take descriptives
  for(i in 1:iterEff){
    y_overarching[i,] <- get_Linear(doseinput,a[i],b[i])
  }

  df_plot <- data.frame(x = xdose,   # use dose as regular scale for plotting
                        median = apply(y_overarching,MARGIN =  2,
                                       FUN= function(x) median(x,
                                                               na.rm = T)),
                        mean = apply(y_overarching,MARGIN =  2,
                                     FUN= function(x) mean(x,
                                                           na.rm = T)),
                        Q5 = apply(y_overarching,MARGIN =  2,
                                   FUN= function(x) quantile(x,0.05,
                                                             na.rm = T)),
                        Q95 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.95,
                                                              na.rm = T)),
                        Q25 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.25,
                                                              na.rm = T)),
                        Q75 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.75,
                                                              na.rm = T))
  )
  
  # # # Plotting method 2: dose ~ RR by median parameters
  # y_median <- get_Linear(doseinput,a = median(a), b = median(b))
  # y_p5 <- get_Linear(x = doseinput, a = quantile(a,0.05), b = quantile(b, 0.05))
  # y_p95 <- get_Linear(x = doseinput, a = quantile(a,0.95), b = quantile(b, 0.95))
  # y_p25 <- get_Linear(x = doseinput, a = quantile(a,0.25), b = quantile(b, 0.25))
  # y_p75 <- get_Linear(x = doseinput, a = quantile(a,0.75), b = quantile(b, 0.75))
  # df_plot <- data.frame(x = xdose,    # use dose at regular scale for plotting
  #                       median = doseinput,Q5 = y_p5,Q95 = y_p95,
  #                       Q25 = y_p25,Q75 = y_p75)


  ## Plotting parameters
  # set axis bounds
  xupper <- ifelse(is.null(xup),   # if a upper bound not specified
                   ceiling (max(
                     median(dose)*1.1,
                     median(dose) + 1,
                     dosehigh
                   )),
                   xup
  )

  yupper <- ifelse(is.null(yup),
                   ceiling(
                     min(
                       max(DataList$RR.u,na.rm = T)*1.1,
                       max(DataList$RR.u,na.rm = T) +1
                     )
                   ),
                   yup)

  graphics.off()
  myplotbase <- ggplot(data.frame())+ geom_blank()
  
  plotcurves <- myplotbase+
    # title to display reference
    labs(title = "Overarching Curve",
         x = "ADD(ug/kg/D)",
         y = "Relative Risks")+
    # median curve
    geom_line(data = df_plot,
              aes(x = x , y = median),
              color = "blue", size = 2)+
    # # mean curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = mean),
    #           color = "red", size = 2)+
    # # 5% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q5),
    #           color = "grey", size = 0.8,linetype ="dotted")+
    # # 95% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q95),
    #           color = "grey", size = 0.8,linetype ="dotted")+
    # 25% curve
    geom_line(data = df_plot,
            aes(x = x , y = Q25),
            color = "brown", size = 0.8,linetype ="dotted")+
    # 75% curve
    geom_line(data = df_plot,
              aes(x = x , y = Q75),
              color = "brown", size = 0.8,linetype ="dotted")+
    # data points
    geom_point(data = DataList,
               aes(x = dose,y = ymean, group = as.factor(Study),
                   color = as.factor(Study), shape = as.factor(Study)),
               size = 5)+
    scale_shape_manual(values = 1:NStudy)+
    geom_segment(data = DataList,
                 aes(x = dose,xend = dose,
                     y = RR.l,yend = RR.u,
                     group = as.factor(Study), color = as.factor(Study)),
                 size = 1)+
    # # set axis coordinates bound
    coord_cartesian(xlim = c(0,xupper),
                    ylim = c(0,yupper))+
    scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
    scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))+
    theme_classic(base_size = 20)
  
  svg(filename = paste0("Overarching_Linear_",savetag,".svg"),
      width = 15,height = 15)
  suppressWarnings(
    print(plotcurves)
  )
  dev.off()
  
}

## Power-----
plot_specific_Power <- function(stanfit,DataList,
                                aname = "a",bname = "b",gname="g",
                                xup = NULL, yup = NULL,
                                savetag = "Specific"){
  # stanfit posterior
  df_posteriors <- as.data.frame(stanfit)
  iterEff <- nrow(df_posteriors)
  NStudy <- length(unique(DataList$Index))          # number of studies
  G <- as.integer(table(DataList$Index)) # number of dose group in each study
  
  # curves as xy plot
  doselow <- 0
  dosehigh <- ceiling(max(DataList$dose))
  # doselow and dosehigh are theoretical range of doses to be as wide as possible
  xdose <- seq(from = doselow, to = dosehigh,
               by = round((dosehigh - doselow)/1e3,3) )
  # convert xaxis coordinates because parameters were fitted on standarized dose
  mindose <- min(DataList$dose)
  maxdose <- max(DataList$dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose
  doseinput <- dosed
  # doseinput <- xdose
  
  # extract posteriors and plot by study
  y_specific <- matrix(NA, nrow = iterEff,ncol = length(doseinput))
  
  for(s in 1:NStudy){
    
    # calculate study specific y
    a_specific <- df_posteriors[,paste0(aname,"[",s,"]")]
    b_specific <- df_posteriors[,paste0(bname,"[",s,"]")]
    g_specific <- df_posteriors[,paste0(gname,"[",s,"]")]
    
    # Two ways of plotting
    
    # Method 1: dose ~ summary of iterative RR by all parameters
    # each iter of pars fits a distribution of y then take descriptives
    # y_specific is the dataframe of y based on study-specific parameters with simulation
    for(i in 1:iterEff){
      # use standardized dose input
      y_specific[i,] <- get_Power(
        doseinput, a = a_specific[i], b = b_specific[i],g=g_specific[i])
    }
    
    # summary data
    df_plot <- data.frame(x = xdose,
                          median = apply(y_specific,MARGIN =  2,
                                         FUN= function(x) median(x,
                                                                 na.rm = T)),
                          mean = apply(y_specific,MARGIN =  2,
                                       FUN= function(x) mean(x,
                                                             na.rm = T)),
                          Q5 = apply(y_specific,MARGIN =  2,
                                     FUN= function(x) quantile(x,0.05,
                                                               na.rm = T)),
                          Q95 = apply(y_specific,MARGIN =  2,
                                      FUN= function(x) quantile(x,0.95,
                                                                na.rm = T)),
                          Q25 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.25,
                                                                 na.rm = T)),
                          Q75 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.75,
                                                                 na.rm = T))
    )
    
    # ## Method 2: dose ~ RR by median parameters
    # y_median <- get_Power(
    #   doseinput, a = median(a_specific),b = median(b_specific),
    #   g = median(g_specific))
    # y_p5 <- get_Power(
    #   doseinput, a = quantile(a_specific,0.05),b = quantile(b_specific,0.05),
    #   g = quantile(g_specific,0.05))
    # y_p95 <- get_Power(
    #   doseinput, a = quantile(a_specific,0.95),b = quantile(b_specific,0.95),
    #   g = quantile(g_specific,0.95))
    # y_p25 <- get_Power(
    #   doseinput, a = quantile(a_specific,0.25),b = quantile(b_specific,0.25),
    #   g = quantile(g_specific,0.25))
    # y_p75 <- get_Power(
    #   doseinput, a = quantile(a_specific,0.75),b = quantile(b_specific,0.75),
    #   g = quantile(g_specific,0.75))
    # df_plot <- data.frame(
    #   x = xdose,
    #   median = y_median, Q5 = y_p5, Q95 = y_p95,Q25=y_p25,Q75=y_p75
    # )
    
    # extract study specific raw data
    Reference <- DataList$Study[(sum(G[1:s-1])+1)]
    Dose <-  DataList$dose[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR <- DataList$ymean[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l <- DataList$RR.l[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.u <- DataList$RR.u[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l[is.na(RR.l)] <- RR[is.na(RR.l)]
    RR.u[is.na(RR.u)] <- RR[is.na(RR.l)]
    
    # data points of each study
    df_raw_study <- data.frame(
      Dose = Dose,
      RR.l = RR.l,
      RR.u = RR.u,
      RR = RR
    )
    
    ## plotting parameters
    xupper <- ifelse(is.null(xup),   # if a upper bound not specified
                     ceiling (max(
                       median(Dose)*1.1,
                       median(Dose) + 1,
                       max(Dose)
                     )),
                     xup
    )
    
    yupper <- ifelse(is.null(yup),
                     ceiling(
                       min(
                         max(RR.u,na.rm = T)*1.1,
                         max(RR.u,na.rm = T) +1
                       )
                     ),
                     yup)
    
    # get a blank plot base
    graphics.off()
    myplotbase <- ggplot(data.frame())+ geom_blank()
    
    plot_specific_summary <- myplotbase+
      # title to display reference
      labs(title = paste0("Study #",s," ",Reference),
           x = "ADD(ug/kg/D)",
           y = "RR")+
      # median curve
      geom_line(data = df_plot,
                aes(x = x , y = median),
                color = "blue", size = 1)+
      # # # mean curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = mean),
      #           color = "red", size = 2)+
      # 5% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q5),
                color = "green", size = 1,linetype ="dotted")+
      # 95% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q95),
                color = "red", size = 1,linetype ="dotted")+
      # # 25% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q25),
      #           color = "green", size = 0.8,linetype ="dotted")+
      # # 75% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q75),
      #           color = "red", size = 0.8,linetype ="dotted")+
      # data points median
      geom_point(data = df_raw_study,
                 aes(x = Dose,y = RR),
                 size = 5, color = s)+
      # upper and lower bound of RR
      geom_segment(data = df_raw_study,
                   aes(x = Dose,xend = Dose,
                       y = RR.l,yend = RR.u),
                   size = 1, color= s)+
      # set axis coordinates bound
      coord_cartesian(xlim = c(0,xupper),
                      ylim = c(0,yupper))+
      scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
      scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))+
      theme_classic(base_size = 20)
    
    # create figs
    # dev.new("windows",noRStudioGD = T)
    # save as svg
    svg(filename = paste0(
      "Power_",savetag," Study #",s," ",Reference,".svg"
    ),
    width = 15,height = 15)
    print(
      plot_specific_summary
    )
    dev.off()
  }
  
}


plot_overarching_Power <- function(stanfit,DataList,
                                   aname = "a",bname = "b",gname = "g",
                                   xup = NULL, yup = NULL,
                                   savetag = "curve"){
  # posteriors
  df_posteriors <- as.data.frame(stanfit)
  iterEff <- nrow(df_posteriors)
  NStudy <- length(unique(DataList$Index))          # number of studies
  
  # Calculate overarching parameters
  a <- rlnorm(iterEff,df_posteriors$mu_a,df_posteriors$sigma_a)
  b <- extraDistr::rhcauchy(iterEff,df_posteriors$scale_b)
  g <- rlnorm(iterEff,df_posteriors$mu_g,df_posteriors$sigma_g)
  
  # raw data
  dose <- DataList$dose
  doselow <- 0
  dosehigh <- ceiling(max(dose))
  xdose <- seq(from = doselow, to = dosehigh,
               by = round((dosehigh - doselow)/1e3,3))
  # if lower or upper bound of RR is missing, use RR values
  DataList$RR.l[is.na(DataList$RR.l)] <- DataList$RR[is.na(DataList$RR.l)]     
  DataList$RR.u[is.na(DataList$RR.u)] <- DataList$RR[is.na(DataList$RR.l)] 
  
  # convert xaxis coordinates
  mindose <- min(dose)
  maxdose <- max(dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose as input for y calculation
  doseinput <- dosed
  # doseinput <- xdose
  
  # plotting method 1: dose ~ median of iterative RR by all parameters
  y_overarching <- matrix(NA, nrow = iterEff,ncol = length(doseinput))
  
  # each iter of pars fits a distribution of y then take descriptives
  for(i in 1:iterEff){
    y_overarching[i,] <- get_Power(doseinput,a[i],b[i],g[i])
  }
  
  df_plot <- data.frame(x = xdose,   # use dose as regular scale for plotting
                        median = apply(y_overarching,MARGIN =  2,
                                       FUN= function(x) median(x,
                                                               na.rm = T)),
                        mean = apply(y_overarching,MARGIN =  2,
                                     FUN= function(x) mean(x,
                                                           na.rm = T)),
                        Q5 = apply(y_overarching,MARGIN =  2,
                                   FUN= function(x) quantile(x,0.05,
                                                             na.rm = T)),
                        Q95 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.95,
                                                              na.rm = T)),
                        Q25 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.25,
                                                              na.rm = T)),
                        Q75 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.75,
                                                              na.rm = T))
  )
  
  # # # Plotting method 2: dose ~ RR by median parameters
  # y_median <- get_Power(doseinput,a = median(a), b = median(b),g= median(g))
  # y_p5 <- get_Power(x = doseinput, a = quantile(a,0.05),
  #                   b = quantile(b, 0.05),g= quantile(g,0.05))
  # y_p95 <- get_Power(x = doseinput, a = quantile(a,0.95),
  #                    b = quantile(b, 0.95),g=quantile(g,0.95))
  # y_p25 <- get_Power(x = doseinput, a = quantile(a,0.25),
  #                    b = quantile(b, 0.25),g=quantile(g,0.25))
  # y_p75 <- get_Power(x = doseinput, a = quantile(a,0.75),
  #                    b = quantile(b, 0.75),g=quantile(g,0.75))
  # df_plot <- data.frame(x = xdose,    # use dose at regular scale for plotting
  #                       median = doseinput,Q5 = y_p5,Q95 = y_p95,
  #                       Q25 = y_p25,Q75 = y_p75)
  
  
  ## Plotting parameters
  # set axis bounds
  xupper <- ifelse(is.null(xup),   # if a upper bound not specified
                   ceiling (max(
                     median(dose)*1.1,
                     median(dose) + 1,
                     dosehigh
                   )),
                   xup
  )
  
  yupper <- ifelse(is.null(yup),
                   ceiling(
                     min(
                       max(DataList$RR.u,na.rm = T)*1.1,
                       max(DataList$RR.u,na.rm = T) +1
                     )
                   ),
                   yup)
  
  graphics.off()
  myplotbase <- ggplot(data.frame())+ geom_blank()
  
  plotcurves <- myplotbase+
    # title to display reference
    labs(title = "Overarching Curve",
         x = "ADD(ug/kg/D)",
         y = "Relative Risks")+
    # median curve
    geom_line(data = df_plot,
              aes(x = x , y = median),
              color = "blue", size = 2)+
    # # mean curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = mean),
    #           color = "red", size = 2)+
    # # 5% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q5),
    #           color = "grey", size = 0.8,linetype ="dotted")+
    # # 95% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q95),
  #           color = "grey", size = 0.8,linetype ="dotted")+
  # 25% curve
  geom_line(data = df_plot,
            aes(x = x , y = Q25),
            color = "brown", size = 0.8,linetype ="dotted")+
    # 75% curve
    geom_line(data = df_plot,
              aes(x = x , y = Q75),
              color = "brown", size = 0.8,linetype ="dotted")+
    # data points
    geom_point(data = DataList,
               aes(x = dose,y = ymean, group = as.factor(Study),
                   color = as.factor(Study), shape = as.factor(Study)),
               size = 5)+
    scale_shape_manual(values = 1:NStudy)+
    geom_segment(data = DataList,
                 aes(x = dose,xend = dose,
                     y = RR.l,yend = RR.u,
                     group = as.factor(Study), color = as.factor(Study)),
                 size = 1)+
    # # set axis coordinates bound
    coord_cartesian(xlim = c(0,xupper),
                    ylim = c(0,yupper))+
    scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
    scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))+
    theme_classic(base_size = 20)
  
  svg(filename = paste0("Overarching_Power_",savetag,".svg"),
      width = 15,height = 15)
  suppressWarnings(
    print(plotcurves)
  )
  dev.off()
  
}

## MM-----
plot_specific_MM <- function(stanfit,DataList,
                             aname = "a",bname = "b",cname="c",
                             xup = NULL, yup = NULL,
                             savetag = "Specific"){
  # stanfit posterior
  df_posteriors <- as.data.frame(stanfit)
  iterEff <- nrow(df_posteriors)
  NStudy <- length(unique(DataList$Index))          # number of studies
  G <- as.integer(table(DataList$Index)) # number of dose group in each study
  
  # curves as xy plot
  doselow <- 0
  dosehigh <- ceiling(max(DataList$dose))
  # doselow and dosehigh are theoretical range of doses to be as wide as possible
  xdose <- seq(from = doselow, to = dosehigh,
               by = round((dosehigh - doselow)/1e3,3) )
  # convert xaxis coordinates because parameters were fitted on standarized dose
  mindose <- min(DataList$dose)
  maxdose <- max(DataList$dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose
  doseinput <- dosed
  # doseinput <- xdose
  
  # extract posteriors and plot by study
  y_specific <- matrix(NA, nrow = iterEff,ncol = length(doseinput))
  
  for(s in 1:NStudy){
    
    # calculate study specific y
    a_specific <- df_posteriors[,paste0(aname,"[",s,"]")]
    b_specific <- df_posteriors[,paste0(bname,"[",s,"]")]
    c_specific <- df_posteriors[,paste0(cname,"[",s,"]")]
    
    # Two ways of plotting
    
    # Method 1: dose ~ summary of iterative RR by all parameters
    # each iter of pars fits a distribution of y then take descriptives
    # y_specific is the dataframe of y based on study-specific parameters with simulation
    for(i in 1:iterEff){
      # use standardized dose input
      y_specific[i,] <- get_MM(
        doseinput, a = a_specific[i], b = b_specific[i],c=c_specific[i])
    }
    
    # summary data
    df_plot <- data.frame(x = xdose,
                          median = apply(y_specific,MARGIN =  2,
                                         FUN= function(x) median(x,
                                                                 na.rm = T)),
                          mean = apply(y_specific,MARGIN =  2,
                                       FUN= function(x) mean(x,
                                                             na.rm = T)),
                          Q5 = apply(y_specific,MARGIN =  2,
                                     FUN= function(x) quantile(x,0.05,
                                                               na.rm = T)),
                          Q95 = apply(y_specific,MARGIN =  2,
                                      FUN= function(x) quantile(x,0.95,
                                                                na.rm = T)),
                          Q25 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.25,
                                                                 na.rm = T)),
                          Q75 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.75,
                                                                 na.rm = T))
    )
    
    # ## Method 2: dose ~ RR by median parameters
    # y_median <- get_MM(
    #   doseinput, a = median(a_specific),b = median(b_specific),
    #   c = median(c_specific))
    # y_p5 <- get_MM(
    #   doseinput, a = quantile(a_specific,0.05),b = quantile(b_specific,0.05),
    #   c = quantile(c_specific,0.05))
    # y_p95 <- get_MM(
    #   doseinput, a = quantile(a_specific,0.95),b = quantile(b_specific,0.95),
    #   c = quantile(c_specific,0.95))
    # y_p25 <- get_MM(
    #   doseinput, a = quantile(a_specific,0.25),b = quantile(b_specific,0.25),
    #   c = quantile(c_specific,0.25))
    # y_p75 <- get_MM(
    #   doseinput, a = quantile(a_specific,0.75),b = quantile(b_specific,0.75),
    #   c = quantile(c_specific,0.75))
    # df_plot <- data.frame(
    #   x = xdose,
    #   median = y_median, Q5 = y_p5, Q95 = y_p95,Q25=y_p25,Q75=y_p75
    # )
    
    # extract study specific raw data
    Reference <- DataList$Study[(sum(G[1:s-1])+1)]
    Dose <-  DataList$dose[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR <- DataList$ymean[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l <- DataList$RR.l[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.u <- DataList$RR.u[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l[is.na(RR.l)] <- RR[is.na(RR.l)]
    RR.u[is.na(RR.u)] <- RR[is.na(RR.l)]
    
    # data points of each study
    df_raw_study <- data.frame(
      Dose = Dose,
      RR.l = RR.l,
      RR.u = RR.u,
      RR = RR
    )
    
    ## plotting parameters
    xupper <- ifelse(is.null(xup),   # if a upper bound not specified
                     ceiling (max(
                       median(Dose)*1.1,
                       median(Dose) + 1,
                       max(Dose)
                     )),
                     xup
    )
    
    yupper <- ifelse(is.null(yup),
                     ceiling(
                       min(
                         max(RR.u,na.rm = T)*1.1,
                         max(RR.u,na.rm = T) +1
                       )
                     ),
                     yup)
    
    # get a blank plot base
    graphics.off()
    myplotbase <- ggplot(data.frame())+ geom_blank()
    
    plot_specific_summary <- myplotbase+
      # title to display reference
      labs(title = paste0("Study #",s," ",Reference),
           x = "ADD(ug/kg/D)",
           y = "RR")+
      # median curve
      geom_line(data = df_plot,
                aes(x = x , y = median),
                color = "blue", size = 1)+
      # # # mean curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = mean),
      #           color = "red", size = 2)+
      # 5% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q5),
                color = "green", size = 1,linetype ="dotted")+
      # 95% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q95),
                color = "red", size = 1,linetype ="dotted")+
      # # 25% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q25),
      #           color = "green", size = 0.8,linetype ="dotted")+
      # # 75% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q75),
      #           color = "red", size = 0.8,linetype ="dotted")+
      # data points median
      geom_point(data = df_raw_study,
                 aes(x = Dose,y = RR),
                 size = 5, color = s)+
      # upper and lower bound of RR
      geom_segment(data = df_raw_study,
                   aes(x = Dose,xend = Dose,
                       y = RR.l,yend = RR.u),
                   size = 1, color= s)+
      # set axis coordinates bound
      coord_cartesian(xlim = c(0,xupper),
                      ylim = c(0,yupper))+
      scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
      scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))+
      theme_classic(base_size = 20)
    
    # create figs
    # dev.new("windows",noRStudioGD = T)
    # save as svg
    svg(filename = paste0(
      "MM_",savetag," Study #",s," ",Reference,".svg"
    ),
    width = 15,height = 15)
    print(
      plot_specific_summary
    )
    dev.off()
  }
  
}


plot_overarching_MM <- function(stanfit,DataList,
                                aname = "a",bname = "b",cname = "c",
                                xup = NULL, yup = NULL,
                                savetag = "curve"){
  # posteriors
  df_posteriors <- as.data.frame(stanfit)
  iterEff <- nrow(df_posteriors)
  NStudy <- length(unique(DataList$Index))          # number of studies
  
  # Calculate overarching parameters
  a <- rlnorm(iterEff,df_posteriors$mu_a,df_posteriors$sigma_a)
  b <- extraDistr::rhcauchy(iterEff,df_posteriors$scale_b)
  c <- rlnorm(iterEff,df_posteriors$mu_c,df_posteriors$sigma_c)
  
  # raw data
  dose <- DataList$dose
  doselow <- 0
  dosehigh <- ceiling(max(dose))
  xdose <- seq(from = doselow, to = dosehigh,
               by = round((dosehigh - doselow)/1e3,3))
  # if lower or upper bound of RR is missing, use RR values
  DataList$RR.l[is.na(DataList$RR.l)] <- DataList$RR[is.na(DataList$RR.l)]     
  DataList$RR.u[is.na(DataList$RR.u)] <- DataList$RR[is.na(DataList$RR.l)] 
  
  # convert xaxis coordinates
  mindose <- min(dose)
  maxdose <- max(dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose as input for y calculation
  doseinput <- dosed
  # doseinput <- xdose
  
  # plotting method 1: dose ~ median of iterative RR by all parameters
  y_overarching <- matrix(NA, nrow = iterEff,ncol = length(doseinput))
  
  # each iter of pars fits a distribution of y then take descriptives
  for(i in 1:iterEff){
    y_overarching[i,] <- get_MM(doseinput,a[i],b[i],c[i])
  }
  
  df_plot <- data.frame(x = xdose,   # use dose as regular scale for plotting
                        median = apply(y_overarching,MARGIN =  2,
                                       FUN= function(x) median(x,
                                                               na.rm = T)),
                        mean = apply(y_overarching,MARGIN =  2,
                                     FUN= function(x) mean(x,
                                                           na.rm = T)),
                        Q5 = apply(y_overarching,MARGIN =  2,
                                   FUN= function(x) quantile(x,0.05,
                                                             na.rm = T)),
                        Q95 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.95,
                                                              na.rm = T)),
                        Q25 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.25,
                                                              na.rm = T)),
                        Q75 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.75,
                                                              na.rm = T))
  )
  
  # # # Plotting method 2: dose ~ RR by median parameters
  # y_median <- get_MM(doseinput,a = median(a), b = median(b),c=median(c))
  # y_p5 <- get_MM(x = doseinput, a = quantile(a,0.05),
  #                b = quantile(b, 0.05),c=quantile(c,0.05))
  # y_p95 <- get_MM(x = doseinput, a = quantile(a,0.95),
  #                 b = quantile(b, 0.95),c=quantile(c,0.95))
  # y_p25 <- get_MM(x = doseinput, a = quantile(a,0.25),
  #                 b = quantile(b, 0.25),c=quantile(c,0.25))
  # y_p75 <- get_MM(x = doseinput, a = quantile(a,0.75),
  #                 b = quantile(b, 0.75),c=quantile(c,0.75))
  # df_plot <- data.frame(x = xdose,    # use dose at regular scale for plotting
  #                       median = doseinput,Q5 = y_p5,Q95 = y_p95,
  #                       Q25 = y_p25,Q75 = y_p75)
  
  
  ## Plotting parameters
  # set axis bounds
  xupper <- ifelse(is.null(xup),   # if a upper bound not specified
                   ceiling (max(
                     median(dose)*1.1,
                     median(dose) + 1,
                     dosehigh
                   )),
                   xup
  )
  
  yupper <- ifelse(is.null(yup),
                   ceiling(
                     min(
                       max(DataList$RR.u,na.rm = T)*1.1,
                       max(DataList$RR.u,na.rm = T) +1
                     )
                   ),
                   yup)
  
  graphics.off()
  myplotbase <- ggplot(data.frame())+ geom_blank()
  
  plotcurves <- myplotbase+
    # title to display reference
    labs(title = "Overarching Curve",
         x = "ADD(ug/kg/D)",
         y = "Relative Risks")+
    # median curve
    geom_line(data = df_plot,
              aes(x = x , y = median),
              color = "blue", size = 2)+
    # # mean curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = mean),
    #           color = "red", size = 2)+
    # # 5% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q5),
    #           color = "grey", size = 0.8,linetype ="dotted")+
    # # 95% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q95),
  #           color = "grey", size = 0.8,linetype ="dotted")+
  # 25% curve
  geom_line(data = df_plot,
            aes(x = x , y = Q25),
            color = "brown", size = 0.8,linetype ="dotted")+
    # 75% curve
    geom_line(data = df_plot,
              aes(x = x , y = Q75),
              color = "brown", size = 0.8,linetype ="dotted")+
    # data points
    geom_point(data = DataList,
               aes(x = dose,y = ymean, group = as.factor(Study),
                   color = as.factor(Study), shape = as.factor(Study)),
               size = 5)+
    scale_shape_manual(values = 1:NStudy)+
    geom_segment(data = DataList,
                 aes(x = dose,xend = dose,
                     y = RR.l,yend = RR.u,
                     group = as.factor(Study), color = as.factor(Study)),
                 size = 1)+
    # # set axis coordinates bound
    coord_cartesian(xlim = c(0,xupper),
                    ylim = c(0,yupper))+
    scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
    scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))+
    theme_classic(base_size = 20)
  
  svg(filename = paste0("Overarching_MM_",savetag,".svg"),
      width = 15,height = 15)
  suppressWarnings(
    print(plotcurves)
  )
  dev.off()
  
}


## Hill-----
plot_specific_Hill <- function(stanfit,DataList,
                               aname = "a",bname = "b",cname="c",gname="g",
                               xup = NULL, yup = NULL,
                               savetag = "Specific"){
  # stanfit posterior
  df_posteriors <- as.data.frame(stanfit)
  iterEff <- nrow(df_posteriors)
  NStudy <- length(unique(DataList$Index))          # number of studies
  G <- as.integer(table(DataList$Index)) # number of dose group in each study
  
  # curves as xy plot
  doselow <- 0
  dosehigh <- ceiling(max(DataList$dose))
  # doselow and dosehigh are theoretical range of doses to be as wide as possible
  xdose <- seq(from = doselow, to = dosehigh,
               by = round((dosehigh - doselow)/1e3,3) )
  # convert xaxis coordinates because parameters were fitted on standarized dose
  mindose <- min(DataList$dose)
  maxdose <- max(DataList$dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose
  doseinput <- dosed
  # doseinput <- xdose
  
  # extract posteriors and plot by study
  y_specific <- matrix(NA, nrow = iterEff,ncol = length(doseinput))
  
  for(s in 1:NStudy){
    
    # calculate study specific y
    a_specific <- df_posteriors[,paste0(aname,"[",s,"]")]
    b_specific <- df_posteriors[,paste0(bname,"[",s,"]")]
    c_specific <- df_posteriors[,paste0(cname,"[",s,"]")]
    g_specific <- df_posteriors[,paste0(gname,"[",s,"]")]
    
    # Two ways of plotting
    
    # Method 1: dose ~ summary of iterative RR by all parameters
    # each iter of pars fits a distribution of y then take descriptives
    # y_specific is the dataframe of y based on study-specific parameters with simulation
    for(i in 1:iterEff){
      # use standardized dose input
      y_specific[i,] <- get_Hill(
        doseinput, a = a_specific[i], b = b_specific[i],c=c_specific[i],
      g = g_specific[i])
    }
    
    # summary data
    df_plot <- data.frame(x = xdose,
                          median = apply(y_specific,MARGIN =  2,
                                         FUN= function(x) median(x,
                                                                 na.rm = T)),
                          mean = apply(y_specific,MARGIN =  2,
                                       FUN= function(x) mean(x,
                                                             na.rm = T)),
                          Q5 = apply(y_specific,MARGIN =  2,
                                     FUN= function(x) quantile(x,0.05,
                                                               na.rm = T)),
                          Q95 = apply(y_specific,MARGIN =  2,
                                      FUN= function(x) quantile(x,0.95,
                                                                na.rm = T)),
                          Q25 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.25,
                                                                 na.rm = T)),
                          Q75 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.75,
                                                                 na.rm = T))
    )
    
    # ## Method 2: dose ~ RR by median parameters
    # y_median <- get_Hill(
    #   doseinput, a = median(a_specific),b = median(b_specific),
    #   c = median(c_specific), g = median(g_specific))
    # y_p5 <- get_Hill(
    #   doseinput, a = quantile(a_specific,0.05),b = quantile(b_specific,0.05),
    #   c = quantile(c_specific,0.05),g = quantile(g_specific,0.05))
    # y_p95 <- get_Hill(
    #   doseinput, a = quantile(a_specific,0.95),b = quantile(b_specific,0.95),
    #   c = quantile(c_specific,0.95),g = quantile(g_specific,0.95))
    # y_p25 <- get_Hill(
    #   doseinput, a = quantile(a_specific,0.25),b = quantile(b_specific,0.25),
    #   c = quantile(c_specific,0.25), g = quantile(g_specific,0.25))
    # y_p75 <- get_Hill(
    #   doseinput, a = quantile(a_specific,0.75),b = quantile(b_specific,0.75),
    #   c = quantile(c_specific,0.75), g = quantile(g_specific,0.75))
    # df_plot <- data.frame(
    #   x = xdose,
    #   median = y_median, Q5 = y_p5, Q95 = y_p95,Q25=y_p25,Q75=y_p75
    # )
    
    # extract study specific raw data
    Reference <- DataList$Study[(sum(G[1:s-1])+1)]
    Dose <-  DataList$dose[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR <- DataList$ymean[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l <- DataList$RR.l[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.u <- DataList$RR.u[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l[is.na(RR.l)] <- RR[is.na(RR.l)]
    RR.u[is.na(RR.u)] <- RR[is.na(RR.l)]
    
    # data points of each study
    df_raw_study <- data.frame(
      Dose = Dose,
      RR.l = RR.l,
      RR.u = RR.u,
      RR = RR
    )
    
    ## plotting parameters
    xupper <- ifelse(is.null(xup),   # if a upper bound not specified
                     ceiling (max(
                       median(Dose)*1.1,
                       median(Dose) + 1,
                       max(Dose)
                     )),
                     xup
    )
    
    yupper <- ifelse(is.null(yup),
                     ceiling(
                       min(
                         max(RR.u,na.rm = T)*1.1,
                         max(RR.u,na.rm = T) +1
                       )
                     ),
                     yup)
    
    # get a blank plot base
    graphics.off()
    myplotbase <- ggplot(data.frame())+ geom_blank()
    
    plot_specific_summary <- myplotbase+
      # title to display reference
      labs(title = paste0("Study #",s," ",Reference),
           x = "ADD(ug/kg/D)",
           y = "RR")+
      # median curve
      geom_line(data = df_plot,
                aes(x = x , y = median),
                color = "blue", size = 1)+
      # # # mean curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = mean),
      #           color = "red", size = 2)+
      # 5% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q5),
                color = "green", size = 1,linetype ="dotted")+
      # 95% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q95),
                color = "red", size = 1,linetype ="dotted")+
      # # 25% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q25),
      #           color = "green", size = 0.8,linetype ="dotted")+
      # # 75% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q75),
      #           color = "red", size = 0.8,linetype ="dotted")+
      # data points median
      geom_point(data = df_raw_study,
                 aes(x = Dose,y = RR),
                 size = 5, color = s)+
      # upper and lower bound of RR
      geom_segment(data = df_raw_study,
                   aes(x = Dose,xend = Dose,
                       y = RR.l,yend = RR.u),
                   size = 1, color= s)+
      # set axis coordinates bound
      coord_cartesian(xlim = c(0,xupper),
                      ylim = c(0,yupper))+
      scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
      scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))+
      theme_classic(base_size = 20)
    
    # create figs
    # dev.new("windows",noRStudioGD = T)
    # save as svg
    svg(filename = paste0(
      "Hill_",savetag," Study #",s," ",Reference,".svg"
    ),
    width = 15,height = 15)
    print(
      plot_specific_summary
    )
    dev.off()
  }
  
}

plot_overarching_Hill <- function(stanfit,DataList,
                                  aname = "a",bname = "b",cname = "c",gname="g",
                                  xup = NULL, yup = NULL,
                                  savetag = "curve"){
  # posteriors
  df_posteriors <- as.data.frame(stanfit)
  iterEff <- nrow(df_posteriors)
  NStudy <- length(unique(DataList$Index))          # number of studies
  
  # Calculate overarching parameters
  a <- rlnorm(iterEff,df_posteriors$mu_a,df_posteriors$sigma_a)
  b <- extraDistr::rhcauchy(iterEff,df_posteriors$scale_b)
  c <- rlnorm(iterEff,df_posteriors$mu_c,df_posteriors$sigma_c)
  g <- rlnorm(iterEff,df_posteriors$mu_g,df_posteriors$sigma_g)
  
  # raw data
  dose <- DataList$dose
  doselow <- 0
  dosehigh <- ceiling(max(dose))
  xdose <- seq(from = doselow, to = dosehigh,
               by = round((dosehigh - doselow)/1e3,3))
  # if lower or upper bound of RR is missing, use RR values
  DataList$RR.l[is.na(DataList$RR.l)] <- DataList$RR[is.na(DataList$RR.l)]     
  DataList$RR.u[is.na(DataList$RR.u)] <- DataList$RR[is.na(DataList$RR.l)] 
  
  # convert xaxis coordinates
  mindose <- min(dose)
  maxdose <- max(dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose as input for y calculation
  doseinput <- dosed
  # doseinput <- xdose
  
  # plotting method 1: dose ~ median of iterative RR by all parameters
  y_overarching <- matrix(NA, nrow = iterEff,ncol = length(doseinput))
  
  # each iter of pars fits a distribution of y then take descriptives
  for(i in 1:iterEff){
    y_overarching[i,] <- get_Hill(doseinput,a=a[i],b=b[i],c=c[i],g=g[i])
  }
  
  df_plot <- data.frame(x = xdose,   # use dose as regular scale for plotting
                        median = apply(y_overarching,MARGIN =  2,
                                       FUN= function(x) median(x,
                                                               na.rm = T)),
                        mean = apply(y_overarching,MARGIN =  2,
                                     FUN= function(x) mean(x,
                                                           na.rm = T)),
                        Q5 = apply(y_overarching,MARGIN =  2,
                                   FUN= function(x) quantile(x,0.05,
                                                             na.rm = T)),
                        Q95 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.95,
                                                              na.rm = T)),
                        Q25 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.25,
                                                              na.rm = T)),
                        Q75 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.75,
                                                              na.rm = T))
  )
  
  # # # # Plotting method 2: dose ~ RR by median parameters
  # y_median <- get_Hill(doseinput,a = median(a), b = median(b),
  #                      c=median(c),g=median(g))
  # y_p5 <- get_Hill(
  #   x = doseinput, a = quantile(a,0.05),b = quantile(b, 0.05),
  #   c=quantile(c,0.05),g=quantile(g,0.05))
  # y_p95 <- get_Hill(
  #   x = doseinput, a = quantile(a,0.95),b = quantile(b, 0.95),
  #   c=quantile(c,0.95),g=quantile(g,0.95))
  # y_p25 <- get_Hill(
  #   x = doseinput, a = quantile(a,0.25),b = quantile(b, 0.25),
  #   c=quantile(c,0.25),g=quantile(g,0.25))
  # y_p75 <- get_Hill(
  #   x = doseinput, a = quantile(a,0.75),b = quantile(b, 0.75),
  #   c=quantile(c,0.75),g=quantile(g,0.75))
  # df_plot <- data.frame(x = xdose,    # use dose at regular scale for plotting
  #                       median = doseinput,Q5 = y_p5,Q95 = y_p95,
  #                       Q25 = y_p25,Q75 = y_p75)
  
  
  ## Plotting parameters
  # set axis bounds
  xupper <- ifelse(is.null(xup),   # if a upper bound not specified
                   ceiling (max(
                     median(dose)*1.1,
                     median(dose) + 1,
                     dosehigh
                   )),
                   xup
  )
  
  yupper <- ifelse(is.null(yup),
                   ceiling(
                     min(
                       max(DataList$RR.u,na.rm = T)*1.1,
                       max(DataList$RR.u,na.rm = T) +1
                     )
                   ),
                   yup)
  
  graphics.off()
  myplotbase <- ggplot(data.frame())+ geom_blank()
  
  plotcurves <- myplotbase+
    # title to display reference
    labs(title = "Overarching Curve",
         x = "ADD(ug/kg/D)",
         y = "Relative Risks")+
    # median curve
    geom_line(data = df_plot,
              aes(x = x , y = median),
              color = "blue", size = 2)+
    # # mean curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = mean),
    #           color = "red", size = 2)+
    # # 5% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q5),
    #           color = "grey", size = 0.8,linetype ="dotted")+
    # # 95% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q95),
  #           color = "grey", size = 0.8,linetype ="dotted")+
  # 25% curve
  geom_line(data = df_plot,
            aes(x = x , y = Q25),
            color = "brown", size = 0.8,linetype ="dotted")+
    # 75% curve
    geom_line(data = df_plot,
              aes(x = x , y = Q75),
              color = "brown", size = 0.8,linetype ="dotted")+
    # data points
    geom_point(data = DataList,
               aes(x = dose,y = ymean, group = as.factor(Study),
                   color = as.factor(Study), shape = as.factor(Study)),
               size = 5)+
    scale_shape_manual(values = 1:NStudy)+
    geom_segment(data = DataList,
                 aes(x = dose,xend = dose,
                     y = RR.l,yend = RR.u,
                     group = as.factor(Study), color = as.factor(Study)),
                 size = 1)+
    # # set axis coordinates bound
    coord_cartesian(xlim = c(0,xupper),
                    ylim = c(0,yupper))+
    scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
    scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))+
    theme_classic(base_size = 20)
  
  svg(filename = paste0("Overarching_Hill_",savetag,".svg"),
      width = 15,height = 15)
  suppressWarnings(
    print(plotcurves)
  )
  dev.off()
  
}

## Expo5-----
plot_specific_Expo5 <- function(stanfit,DataList,
                                aname = "a",bname = "b",cname="c",dname="d",
                                xup = NULL, yup = NULL,
                                savetag = "Specific"){
  # stanfit posterior
  df_posteriors <- as.data.frame(stanfit)
  iterEff <- nrow(df_posteriors)
  NStudy <- length(unique(DataList$Index))          # number of studies
  G <- as.integer(table(DataList$Index)) # number of dose group in each study
  
  # curves as xy plot
  doselow <- 0
  dosehigh <- ceiling(max(DataList$dose))
  # doselow and dosehigh are theoretical range of doses to be as wide as possible
  xdose <- seq(from = doselow, to = dosehigh,
               by = round((dosehigh - doselow)/1e3,3) )
  # convert xaxis coordinates because parameters were fitted on standarized dose
  mindose <- min(DataList$dose)
  maxdose <- max(DataList$dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose
  doseinput <- dosed
  # doseinput <- xdose
  
  # extract posteriors and plot by study
  y_specific <- matrix(NA, nrow = iterEff,ncol = length(doseinput))
  
  for(s in 1:NStudy){
    
    # calculate study specific y
    a_specific <- df_posteriors[,paste0(aname,"[",s,"]")]
    b_specific <- df_posteriors[,paste0(bname,"[",s,"]")]
    c_specific <- df_posteriors[,paste0(cname,"[",s,"]")]
    d_specific <- df_posteriors[,paste0(dname,"[",s,"]")]
    
    # Two ways of plotting
    
    # Method 1: dose ~ summary of iterative RR by all parameters
    # each iter of pars fits a distribution of y then take descriptives
    # y_specific is the dataframe of y based on study-specific parameters with simulation
    for(i in 1:iterEff){
      # use standardized dose input
      y_specific[i,] <- get_Expo5(
        doseinput, a = a_specific[i], b = b_specific[i],c=c_specific[i],
        d = d_specific[i])
    }
    
    # summary data
    df_plot <- data.frame(x = xdose,
                          median = apply(y_specific,MARGIN =  2,
                                         FUN= function(x) median(x,
                                                                 na.rm = T)),
                          mean = apply(y_specific,MARGIN =  2,
                                       FUN= function(x) mean(x,
                                                             na.rm = T)),
                          Q5 = apply(y_specific,MARGIN =  2,
                                     FUN= function(x) quantile(x,0.05,
                                                               na.rm = T)),
                          Q95 = apply(y_specific,MARGIN =  2,
                                      FUN= function(x) quantile(x,0.95,
                                                                na.rm = T)),
                          Q25 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.25,
                                                                 na.rm = T)),
                          Q75 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.75,
                                                                 na.rm = T))
    )
    
    # ## Method 2: dose ~ RR by median parameters
    # y_median <- get_Expo5(
    #   doseinput, a = median(a_specific),b = median(b_specific),
    #   c = median(c_specific), d = median(d_specific))
    # y_p5 <- get_Expo5(
    #   doseinput, a = quantile(a_specific,0.05),b = quantile(b_specific,0.05),
    #   c = quantile(c_specific,0.05),d = quantile(d_specific,0.05))
    # y_p95 <- get_Expo5(
    #   doseinput, a = quantile(a_specific,0.95),b = quantile(b_specific,0.95),
    #   c = quantile(c_specific,0.95),d = quantile(d_specific,0.95))
    # y_p25 <- get_Expo5(
    #   doseinput, a = quantile(a_specific,0.25),b = quantile(b_specific,0.25),
    #   c = quantile(c_specific,0.25), d = quantile(d_specific,0.25))
    # y_p75 <- get_Expo5(
    #   doseinput, a = quantile(a_specific,0.75),b = quantile(b_specific,0.75),
    #   c = quantile(c_specific,0.75), d = quantile(d_specific,0.75))
    # df_plot <- data.frame(
    #   x = xdose,
    #   median = y_median, Q5 = y_p5, Q95 = y_p95,Q25=y_p25,Q75=y_p75
    # )
    
    # extract study specific raw data
    Reference <- DataList$Study[(sum(G[1:s-1])+1)]
    Dose <-  DataList$dose[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR <- DataList$ymean[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l <- DataList$RR.l[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.u <- DataList$RR.u[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l[is.na(RR.l)] <- RR[is.na(RR.l)]
    RR.u[is.na(RR.u)] <- RR[is.na(RR.l)]
    
    # data points of each study
    df_raw_study <- data.frame(
      Dose = Dose,
      RR.l = RR.l,
      RR.u = RR.u,
      RR = RR
    )
    
    ## plotting parameters
    xupper <- ifelse(is.null(xup),   # if a upper bound not specified
                     ceiling (max(
                       median(Dose)*1.1,
                       median(Dose) + 1,
                       max(Dose)
                     )),
                     xup
    )
    
    yupper <- ifelse(is.null(yup),
                     ceiling(
                       min(
                         max(RR.u,na.rm = T)*1.1,
                         max(RR.u,na.rm = T) +1
                       )
                     ),
                     yup)
    
    # get a blank plot base
    graphics.off()
    myplotbase <- ggplot(data.frame())+ geom_blank()
    
    plot_specific_summary <- myplotbase+
      # title to display reference
      labs(title = paste0("Study #",s," ",Reference),
           x = "ADD(ug/kg/D)",
           y = "RR")+
      # median curve
      geom_line(data = df_plot,
                aes(x = x , y = median),
                color = "blue", size = 1)+
      # # # mean curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = mean),
      #           color = "red", size = 2)+
      # 5% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q5),
                color = "green", size = 1,linetype ="dotted")+
      # 95% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q95),
                color = "red", size = 1,linetype ="dotted")+
      # # 25% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q25),
      #           color = "green", size = 0.8,linetype ="dotted")+
      # # 75% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q75),
      #           color = "red", size = 0.8,linetype ="dotted")+
      # data points median
      geom_point(data = df_raw_study,
                 aes(x = Dose,y = RR),
                 size = 5, color = s)+
      # upper and lower bound of RR
      geom_segment(data = df_raw_study,
                   aes(x = Dose,xend = Dose,
                       y = RR.l,yend = RR.u),
                   size = 1, color= s)+
      # set axis coordinates bound
      coord_cartesian(xlim = c(0,xupper),
                      ylim = c(0,yupper))+
      scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
      scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))+
      theme_classic(base_size = 20)
    
    # create figs
    # dev.new("windows",noRStudioGD = T)
    # save as svg
    svg(filename = paste0(
      "Expo5_",savetag," Study #",s," ",Reference,".svg"
    ),
    width = 15,height = 15)
    print(
      plot_specific_summary
    )
    dev.off()
  }
  
}


plot_overarching_Expo5 <- function(stanfit,DataList,
                                   aname = "a",bname = "b",cname = "c",dname="d",
                                   xup = NULL, yup = NULL,
                                   savetag = "curve"){
  # posteriors
  df_posteriors <- as.data.frame(stanfit)
  iterEff <- nrow(df_posteriors)
  NStudy <- length(unique(DataList$Index))          # number of studies
  
  # Calculate overarching parameters
  a <- rlnorm(iterEff,df_posteriors$mu_a,df_posteriors$sigma_a)
  b <- extraDistr::rhcauchy(iterEff,df_posteriors$scale_b)
  c <- extraDistr::rhcauchy(iterEff,df_posteriors$scale_c)
  d <- rlnorm(iterEff,df_posteriors$mu_d,df_posteriors$sigma_d)
  
  # raw data
  dose <- DataList$dose
  doselow <- 0
  dosehigh <- ceiling(max(dose))
  xdose <- seq(from = doselow, to = dosehigh,
               by = round((dosehigh - doselow)/1e3,3))
  # if lower or upper bound of RR is missing, use RR values
  DataList$RR.l[is.na(DataList$RR.l)] <- DataList$RR[is.na(DataList$RR.l)]     
  DataList$RR.u[is.na(DataList$RR.u)] <- DataList$RR[is.na(DataList$RR.l)] 
  
  # convert xaxis coordinates
  mindose <- min(dose)
  maxdose <- max(dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose as input for y calculation
  doseinput <- dosed
  # doseinput <- xdose
  
  # plotting method 1: dose ~ median of iterative RR by all parameters
  y_overarching <- matrix(NA, nrow = iterEff,ncol = length(doseinput))
  
  # each iter of pars fits a distribution of y then take descriptives
  for(i in 1:iterEff){
    y_overarching[i,] <- get_Expo5(doseinput,a=a[i],b=b[i],c=c[i],d=d[i])
  }
  
  df_plot <- data.frame(x = xdose,   # use dose as regular scale for plotting
                        median = apply(y_overarching,MARGIN =  2,
                                       FUN= function(x) median(x,
                                                               na.rm = T)),
                        mean = apply(y_overarching,MARGIN =  2,
                                     FUN= function(x) mean(x,
                                                           na.rm = T)),
                        Q5 = apply(y_overarching,MARGIN =  2,
                                   FUN= function(x) quantile(x,0.05,
                                                             na.rm = T)),
                        Q95 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.95,
                                                              na.rm = T)),
                        Q25 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.25,
                                                              na.rm = T)),
                        Q75 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.75,
                                                              na.rm = T))
  )
  
  # # # # Plotting method 2: dose ~ RR by median parameters
  # y_median <- get_Expo5(doseinput,a = median(a), b = median(b),
  #                      c=median(c),d=median(d))
  # y_p5 <- get_Expo5(
  #   x = doseinput, a = quantile(a,0.05),b = quantile(b, 0.05),
  #   c=quantile(c,0.05),d=quantile(d,0.05))
  # y_p95 <- get_Expo5(
  #   x = doseinput, a = quantile(a,0.95),b = quantile(b, 0.95),
  #   c=quantile(c,0.95),d=quantile(d,0.95))
  # y_p25 <- get_Expo5(
  #   x = doseinput, a = quantile(a,0.25),b = quantile(b, 0.25),
  #   c=quantile(c,0.25),d=quantile(d,0.25))
  # y_p75 <- get_Expo5(
  #   x = doseinput, a = quantile(a,0.75),b = quantile(b, 0.75),
  #   c=quantile(c,0.75),d=quantile(d,0.75))
  # df_plot <- data.frame(x = xdose,    # use dose at regular scale for plotting
  #                       median = doseinput,Q5 = y_p5,Q95 = y_p95,
  #                       Q25 = y_p25,Q75 = y_p75)
  
  
  ## Plotting parameters
  # set axis bounds
  xupper <- ifelse(is.null(xup),   # if a upper bound not specified
                   ceiling (max(
                     median(dose)*1.1,
                     median(dose) + 1,
                     dosehigh
                   )),
                   xup
  )
  
  yupper <- ifelse(is.null(yup),
                   ceiling(
                     min(
                       max(DataList$RR.u,na.rm = T)*1.1,
                       max(DataList$RR.u,na.rm = T) +1
                     )
                   ),
                   yup)
  
  graphics.off()
  myplotbase <- ggplot(data.frame())+ geom_blank()
  
  plotcurves <- myplotbase+
    # title to display reference
    labs(title = "Overarching Curve",
         x = "ADD(ug/kg/D)",
         y = "Relative Risks")+
    # median curve
    geom_line(data = df_plot,
              aes(x = x , y = median),
              color = "blue", size = 2)+
    # # mean curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = mean),
    #           color = "red", size = 2)+
    # # 5% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q5),
    #           color = "grey", size = 0.8,linetype ="dotted")+
    # # 95% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q95),
  #           color = "grey", size = 0.8,linetype ="dotted")+
  # 25% curve
  geom_line(data = df_plot,
            aes(x = x , y = Q25),
            color = "brown", size = 0.8,linetype ="dotted")+
    # 75% curve
    geom_line(data = df_plot,
              aes(x = x , y = Q75),
              color = "brown", size = 0.8,linetype ="dotted")+
    # data points
    geom_point(data = DataList,
               aes(x = dose,y = ymean, group = as.factor(Study),
                   color = as.factor(Study), shape = as.factor(Study)),
               size = 5)+
    scale_shape_manual(values = 1:NStudy)+
    geom_segment(data = DataList,
                 aes(x = dose,xend = dose,
                     y = RR.l,yend = RR.u,
                     group = as.factor(Study), color = as.factor(Study)),
                 size = 1)+
    # # set axis coordinates bound
    coord_cartesian(xlim = c(0,xupper),
                    ylim = c(0,yupper))+
    scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
    scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))+
    theme_classic(base_size = 20)
  
  svg(filename = paste0("Overarching_Expo5_",savetag,".svg"),
      width = 15,height = 15)
  suppressWarnings(
    print(plotcurves)
  )
  dev.off()
  
}

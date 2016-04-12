
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyServer(function(input, output) {
  
  SimData <- function(Nyears, AgeMax, SigmaR, M, F1, S_a, 
                      SigmaF, Fdynamics, Rdynamics, R0, CVlen, mids, highs, lows,
                      W_a, L_a, Mat_a, Amat, Fequil, comp_sample=1000,
                      nburn=50){
    
    ## SB_t = spawning biomass over time
    ## F_t = fishing mortality over time
    ## Cn_at = number of individuals that die from fishing mortality
    ## N_at = abundance by number at age over time
    
    ##########################
    ## Initial calcs
    ##########################
    
    tyears <- nburn+Nyears
    
    
    ##########################
    ## Random variables
    ##########################
    
    RecDev <- rnorm(tyears, mean=-(SigmaR^2)/2, sd=SigmaR)
    FishDev <- rnorm(tyears, mean=-(SigmaF^2)/2, sd=SigmaF)
    
    ##########################
    ## Data objects
    ##########################
    SB_t <- F_t <- R_t <- rep(NA, tyears)								
    Cn_at <- N_at <- matrix(NA, nrow=AgeMax+1, ncol=tyears)
    
    #####################################
    ## Fishing and recruitment dynamics
    #####################################	
    
    if(Fdynamics=="Ramp") Framp_t <- c(rep(0.01, nburn), "rampup"=seq(F1, Fmax, length=floor(Nyears/2)), 
                                       "peak"=rep(Fmax, floor((Nyears-floor(Nyears/2))/2)), 
                                       "managed"=rep(Fmax/3, Nyears-floor(Nyears/2)-floor((Nyears-floor(Nyears/2))/2)))
    if(Fdynamics=="Constant") Fconstant_t <- rep(Fequil, tyears)
    if(Fdynamics=="Increasing") Finc_t <- c(rep(0, nburn), ((1-0.01)/(Nyears-1))*1:Nyears)
    
    if(Rdynamics=="Pulsed") Rpulse_t <- c(rep(R0, nburn), "initial"=rep(R0, floor(Nyears/3)),
                                          "pulse_up"=rep(R0*3, floor(Nyears/3)), "pulse_down"=rep(R0, Nyears-floor(Nyears/3)))
    if(Rdynamics=="Constant") Rconstant_t <- rep(R0, tyears)
    
    ##########################
    ## Initialization
    ##########################
    if(Fdynamics=="Endogenous") F_t[1] <- 0.01
    if(Fdynamics=="Ramp") F_t[1] <- Framp_t[1]
    if(Fdynamics=="Constant") F_t[1] <- Fconstant_t[1]
    if(Fdynamics=="Increasing") F_t[1] <- Finc_t[1]
    
    R_t[1] <- R0
    
    ## year 1
    for(a in 1:length(L_a)){
      if(a==1){
        N_at[a,1] <- R_t[1]
      }
      if(a>1 & a<length(L_a)){
        N_at[a,1] <- N_at[a-1,1]*exp(-M-F_t[1]*S_a[a-1])
      }
      if(a==length(L_a)){
        N_at[a,1] <- (N_at[a-1,1]*exp(-M-F_t[1]*S_a[a-1]))/(1-exp(-M-F_t[1]*S_a[a-1]))
      }
      
    }
    SB_t[1] <- sum(N_at[,1] * W_a * S_a)
    Cn_at[,1] <- N_at[,1] * (1-exp(-M - F_t[1]*S_a)) * (F_t[1]*S_a)/(M+F_t[1]*S_a)
    
    ##########################
    ## Projection
    ##########################
    Na0 <- rep(NA, length(W_a))
    if(Rdynamics=="Pulsed"){
      R0 <- median(Rpulse_t[-c(1:nburn)])
    }
    Na0[1] <- R0
    for(a in 2:length(W_a)){
      Na0[a] <- R0 * exp(-M*(a-1))
    }
    SB0 <- sum(Na0[-1]*Mat_a[-1]*W_a[-1])
    
    for(y in 2:tyears){
      ## fishing effort and recruitment, not dependent on age structure
      if(Fdynamics=="Endogenous"){
        if(y <= nburn) F_t[y] <- 0.01
        if(y > nburn) F_t[y] <- F_t[y-1]*(SB_t[y-1]/(Fequil*SB0))^Frate * exp(FishDev[y])
      }
      if(Fdynamics=="Ramp"){
        F_t[y] <- Framp_t[y] * exp(FishDev[y])
      }
      if(Fdynamics=="Constant"){
        F_t[y] <- Fconstant_t[y] * exp(FishDev[y])
      }
      if(Fdynamics=="Increasing"){
        F_t[y] <- Finc_t[y] * exp(FishDev[y])
      }
      if(Rdynamics=="Constant"){
        R_t[y] <- Rconstant_t[y] * exp(RecDev[y])
      }
      if(Rdynamics=="Pulsed"){
        R_t[y] <- Rpulse_t[y] * exp(RecDev[y])
      }
      if(Rdynamics=="BH"){
        R_t[y] <- (4 * h * R0 * SB_t[y-1] / ( SB0*(1-h) + SB_t[y-1]*(5*h-1))) * exp(RecDev[y])
      }
      
      ## age-structured dynamics
      for(a in 1:length(L_a)){
        if(a==1){
          N_at[a,y] <- R_t[y]
        }
        if(a>1 & a<length(L_a)){
          N_at[a,y] <- N_at[a-1,y-1]*exp(-M-F_t[y-1]*S_a[a-1])
        }
        if(a==length(L_a)){
          N_at[a,y] <- (N_at[a-1,y-1] + N_at[a,y-1])*exp(-M-F_t[y-1]*S_a[a-1])
        }
      }
      ## spawning biomass
      SB_t[y] <- sum((N_at[,y] * W_a * Mat_a)[-1])
      ## catch
      Cn_at[,y] <- N_at[,y] * (1-exp(-M-F_t[y]*S_a)) * (F_t[y]*S_a)/(M+F_t[y]*S_a)
    }
    Cn_t <- colSums(Cn_at)
    C_t <- Cn_t
    N_t <- colSums(N_at[-1,])
    D_t <- SB_t/SB0
    
    ################################################
    ## Probability being in a length bin given age
    ################################################
    lbprobs <- function(mnl,sdl) return(pnorm(highs,mnl,sdl)-pnorm(lows,mnl,sdl))
    vlprobs <- Vectorize(lbprobs,vectorize.args=c("mnl","sdl"))
    plba <- t(vlprobs(L_a,L_a*CVlen))
    plba <- plba/rowSums(plba)
    
    ################################################
    ## Probability being in harvested at an age
    ################################################
    page <- matrix(ncol=dim(plba)[1], nrow=tyears)
    for (y in 1:tyears) page[y,] <- N_at[,y] * S_a
    page <- page/rowSums(page)    
    
    ################################################
    ## Probability of sampling a given length bin
    ################################################
    plb <- matrix(ncol=length(mids), nrow=tyears)
    for (y in 1:tyears) plb[y,] <- page[y,] %*% plba
    plb <- plb/rowSums(plb)    
    
    #######################
    ## Length frequencies 
    #######################
    LF <- array(0, dim=dim(plb))
    rownames(LF) <- 1:tyears
    for(y in 1:tyears){
      LF[y,] <- rmultinom(n=1, size=comp_sample, prob=plb[y,])
    }
    
    ML_t <- rowMeans(LF*mids)
    
    
    ########################################################
    ## True mean length in vulnerable population each year 
    ########################################################
    L_t <- vector(length=tyears)
    VL <- matrix(NA, nrow=tyears, ncol=ncol(plb))
    for(y in 1:tyears){
      vul_pop <- sum(N_at[,y]*S_a)
      vul_lengths <- sum(vul_pop*plb[y,]*mids)
      VL[y,] <- vul_pop*plb[y,]*mids
      L_t[y] <- vul_lengths/vul_pop
    }
    
    ########################################################
    ## cut out burn-in
    ########################################################
    
    LFout <- LF[-c(1:nburn),]
    rownames(LFout) <- 1:Nyears
    
    C_tout <- C_t[-c(1:nburn)]
    
    R_tout <- R_t[-c(1:nburn)]
    N_tout <- N_t[-c(1:nburn)]
    SB_tout <- SB_t[-c(1:nburn)]
    D_tout <- D_t[-c(1:nburn)]
    F_tout <- F_t[-c(1:nburn)]
    L_tout <- L_t[-c(1:nburn)]
    N_atout <- N_at[,-c(1:nburn)]
    VL_out <- VL[-c(1:nburn),]
    
    
    DataList <- list("C_t"=C_tout, "VL"=VL_out,
                     "LF"=LFout, "SigmaR"=SigmaR, "R_t"=R_tout, "N_t"=N_tout, "SB_t"=SB_tout, "D_t"=D_tout, "F_t"=F_tout, 
                     "L_t"=L_tout, "N_at"=N_atout, "Amat"=Amat, "Mat_a"=Mat_a, "SB0"=SB0, "Nyears"=Nyears, "L_a"=L_a,
                     "W_a"=W_a, "AgeMax"=AgeMax, "M"=M, "S_a"=S_a, "plb"=plb, "plba"=plba, "page"=page, "R0"=R0)
    
    return(DataList)
  }
  
  
  get_sim <- function(error){
    
    Sa50 <- round(input$t0 - log(1-(input$S50/input$linf))/input$vbk)
    Ma50 <- round(input$t0 - log(1-(input$M50/input$linf))/input$vbk)
    
    S_a <- 1/(1+exp(Sa50 - 0:input$AgeMax))
    Mat_a <- 1/(1+exp(Ma50 - 0:input$AgeMax))
    L_a <- input$linf*(1-exp(-input$vbk*(0:input$AgeMax - input$t0)))
    W_a <- input$lwa * L_a ^ input$lwb
    binwidth <- 1
    mids <- seq((binwidth/2), input$linf*1.2, by=binwidth)
    highs <- mids + (binwidth/2)
    lows <- mids - (binwidth/2)
    ages <- c(0:input$AgeMax)
    
    if(error==TRUE){
      SigmaR <- 0.6
      SigmaF <- 0.3
    }
    if(error==FALSE){
      SigmaR <- SigmaF <- 0
    }
    
    set.seed(1234)
    simdata <- SimData(Nyears=input$nyears, AgeMax=input$AgeMax, 
                       SigmaR=SigmaR, M=input$M, F1=0.01, S_a=S_a, 
                       SigmaF=SigmaF, W_a=W_a, L_a=L_a, Mat_a=Mat_a, 
                       Amat=Ma50, Fdynamics=input$f_options, 
                       Rdynamics=input$rec_options, Fequil=0.2, R0=1000,
                       CVlen=0.2, highs=highs, mids=mids, lows=lows)

    
    return(simdata)
    
  }

  output$MLplot <- renderPlot({
    pop <- get_sim(error=FALSE)
    L_t <- get_sim(error=FALSE)$L_t
    rL_t <- L_t/L_t[1]
    if(input$error==TRUE){
      L_t2 <- get_sim(error=TRUE)$L_t
      rL_t2 <- L_t2/L_t2[1]
    }
    years <- 1:input$nyears
    plot(years, rL_t, lwd=4, pch=19, type="o", ylim=c(0, max(rL_t)*1.5),
         ylab="Observed mean length in catch", col=gray(0.2), xlab="Year", yaxt="n")
    axis(2, at=pretty(c(0, max(rL_t)*1.5)), las=2)
    if(input$error==TRUE) lines(rL_t2, lty=2, col=gray(0.2), lwd=2)
  })
  
  output$Fplot <- renderPlot({
    F_t <- get_sim(error=FALSE)$F_t
    if(input$error==TRUE) F_t2 <- get_sim(error=TRUE)$F_t
    years <- 1:input$nyears
    plot(years, F_t, lwd=4, pch=19, type="o", ylim=c(0, 1.5),
         ylab="True fishing mortality", col="forestgreen", xlab="Year", yaxt="n")
    axis(2, at=pretty(c(0,1.5)), las=2)
    if(input$error==TRUE) lines(F_t2, lty=2, col="forestgreen", lwd=2)
  })
  
  output$Rplot <- renderPlot({
    R_t <- get_sim(error=FALSE)$R_t
    if(input$error==TRUE) R_t2 <- get_sim(error=TRUE)$R_t
    years <- 1:input$nyears
    plot(years, R_t, lwd=4, pch=19, type="o", ylim=c(0, 5000),
         ylab="True recruitment", col="steelblue", xlab="Year", yaxt="n")
    axis(2, at=pretty(c(0,5000)), las=2)
    if(input$error==TRUE) lines(R_t2, lty=2, col="steelblue", lwd=2)
  })
  
  output$LFplot <- renderPlot({
    VL <- get_sim(error=FALSE)$VL
    relVL1 <- VL[1,]/max(VL[1,])
    relVL2 <- VL[round(input$nyears)/2,]/max(VL[round(input$nyears)/2,])
    
    plot(x=1, y=1, type="n", xlim=c(0, 90), ylim=c(0, 1), axes=F, xaxs="i", yaxs="i", ann=F)
    polygon(x=c(seq(0,90,length=length(relVL1)), seq(90, 0, length=length(relVL1))),
            y=c(rep(0, length(relVL1)), rev(relVL1)), col="#AAAAAA70")
    polygon(x=c(seq(0,90,length=length(relVL2)), seq(90,0,length=length(relVL2))),
            y=c(rep(0, length(relVL2)), rev(relVL2)), col="#AA000070")
    axis(1, at=pretty(c(0,90)))
    mtext(side=1, "Length", line=3)
    axis(2, at=pretty(c(0,1)), las=2)
    mtext(side=2, "Relative proportion", line=3)
    
#     barplot(VL[1,]/max(VL[1,]), col="#AAAAAA70", xlim=c(0,80), ylim=c(0,1))
#     par(new=TRUE)
#     barplot(VL[10,]/max(VL[10,]), col="#AA000070", xlim=c(0,80), ylim=c(0,1))
    
  })

  output$LFlegend <- renderPlot({
    
    plot(x=1, y=1, type="n", axes=F, ann=F, xlim=c(0,10), ylim=c(0,10))
    legend("top", legend=c("Unfished equilibrium", "Perturbed"),
           col=c("#AAAAAA70", "#AA000070"), border=NA, cex=2.2, pch=15)
    
  })
  


  

  
})

# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyServer(function(input, output) {
  
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
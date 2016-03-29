
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
source("functions.R")

shinyUI(fluidPage(

  # Application title
  titlePanel("Impacts of Fishing and Recruitment on Length"),

  column(2,
        h4("Population dynamic processes"),
        sliderInput("nyears", "Number of years:", value=20, min=0, max=100, step=1),
        sliderInput("AgeMax", "Age classes:", value=14, min=0, max=100, step=1), 
        sliderInput("linf", "Linf:", value=37, min=0, max=200, step=1),
        sliderInput("vbk", "k:", value = 0.2, min=0, max=1, step=0.05),
        sliderInput("t0", "t0:", value = -0.2, min=-1, max=-1e-3, step=0.01),
        sliderInput("M", "M:", value=0.36, min=0, max=1, step=0.01),
        sliderInput("S50", "S50:", value=8, min=0, max=200, step=1),
        sliderInput("M50", "M50:", value=4, min=0, max=200, step=1)
      ),
    column(10,
      sidebarLayout(
        sidebarPanel(
          h4("Recruitment dynamics"),
          selectInput("rec_options", "Recruitment pattern", choices=c("Constant", "Pulsed")),
          h4("Fishing dynamics"),
          selectInput("f_options", "Fishing mortality pattern", choices=c("Constant", "Increasing")),
          conditionalPanel(
            "$('li.active a').first().html()=='Time Series'",
            checkboxInput("error", "Variability?", FALSE)
          )
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Length frequency",
                     #plotOutput("LFplot"),
                     plotOutput("LFlegend")
            ),
            tabPanel("Time Series",
                     h4("Truth"),
                     column(6, plotOutput("Fplot")),
                     column(6, plotOutput("Rplot")),
                     h4("Observed"),
                     column(6, plotOutput("MLplot")),
                     column(6, plotOutput("LFplot"))
            )
            
          )
        )          
        )
      )

))

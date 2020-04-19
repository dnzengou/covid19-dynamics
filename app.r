## source: https://github.com/EconometricsBySimulation/Shiny-Demos/tree/master/Ebola-Dynamic
## inspired by Ebola Model

library(shiny)

## # Define UI for app.r 

# Define UI for slider demo application
ui <-pageWithSidebar(
  
  #  Application title
  headerPanel("COVID-19 Model"),
  
  # Sidebar with sliders that demonstrate various available options
  sidebarPanel(
    
    tags$h3("Data Generation"),
    
    sliderInput("N", "Population:", 
                min=10^4, max=10^7, value=10^6, step=10^4),
    
    sliderInput("IC", "Initial Infected Detected:", 
                min=0, max=10^3, value=100, step=1),
    
    sliderInput("np", "Number of days:", 
                min=30, max=360, value=270, step=10),
    
    textOutput("r0"), 
    
    sliderInput("P", "Transmition rate:", 
                min=.01, max=.25, value=.081, step=.001),
    
    sliderInput("Days", "Exposed Period (Days):", 
                min=0, max=40, value=18, step=1), 
    
    sliderInput("M", "Moralitity Rate:", 
                min=0, max=1, value=.6, step=.05), 
    
    sliderInput("K", "Social 'adaption' to reduce infection:", 
                min=-.0012, max=.0012, value=.0003, step=.0003), 
    
    sliderInput("DET", "Daily Detection Rate:", 
                min=0, max=1, value=.07, step=.01), 
    textOutput("LDET"), 
    
    sliderInput("bed0", "# of Quarantine Beds Available Initially:", 
                min=0, max=10^3, value=0, step=10), 
    
    sliderInput("bed1", "# of New Quarantine Beds Available at 1 months:", 
                min=0, max=10^3, value=10, step=10), 
    
    sliderInput("bed2", "# of New Quarantine Beds Available at 2 months:", 
                min=0, max=10^3, value=50, step=10), 
    
    sliderInput("bed3", "# of New Quarantine Beds Available at 3 months:", 
                min=0, max=10^3, value=50, step=10), 
    
    sliderInput("bed4", "# of New Quarantine Beds Available at 4 months:", 
                min=0, max=10^3, value=50, step=10), 
    
    sliderInput("bed5", "# of New Quarantine Beds Available at 5 months:", 
                min=0, max=10^3, value=50, step=10), 
    
    sliderInput("bed6", "# of New Quarantine Beds Available at 6 months:", 
                min=0, max=10^3, value=50, step=10), 
    
    sliderInput("bed7", "# of New Quarantine Beds Available at 7 months:", 
                min=0, max=10^3, value=50, step=10), 
    
    sliderInput("bed8", "# of New Quarantine Beds Available at 8 months:", 
                min=0, max=10^3, value=100, step=10), 
    
    sliderInput("bed9", "# of New Quarantine Beds Available at 9 months:", 
                min=0, max=10^3, value=100, step=10), 
    
    br(),
    
    h5("Created by:"),
    tags$a("@desireyavro", 
           href="https://desireyavro.com"),
    h5("Inspired by Econometrics by Simulation http://www.econometricsbysimulation.com"),
    h5(textOutput(""))
    
  ),
  
  # Show a table summarizing the values entered
  mainPanel(
    checkboxGroupInput("Indicators", "",
                       c("Susceptible", 
                         "Exposed", 
                         "Recovered",
                         "Deceased",
                         "Infected"),
                       selected=c(
                         "Exposed", 
                         "Infected",
                         "Recovered",
                         "Deceased"),
                       inline=TRUE),
    plotOutput("graph1"),
    plotOutput("graph2")
  )
)


## Define server logic reIuired
#library("shiny")
library("ggplot2")
library("scales")
#library("rmongodb"); load("mongodb-login.RData")

# Simulation and Shiny Application of Flue Season Dynamics
server <- function(input, output) {
  
  mydata <- reactive({
    # Model Parameters:
    
    IC <- input$IC  # Infected
    N  <- input$N   # Total Population
    np <- input$np  # Time periods
    
    # Infection Parameters:
    # Mortality rate Days (or epsilon)
    Mr <- input$M  
    
    # Days till resolution
    Days <- input$Days
    
    # Resolution rate per day (or gamma)
    Dr <- 1/Days 
    
    # Transmition rate per day (for those Exposed), (or beta)
    P <- input$P
    
    # Social adaption to disease rt=r0d/(1+S)^t
    
    M   <- input$M
    DET <- input$DET
    K  <- input$K
    
    # Gain in number of beds available
    bedsv <- input$bed0
    for (i in 1:9) bedsv[i+1] <- input[[paste0('bed',i)]]+bedsv[i]
    
    # Model - Dynamic Change
    
    ## Let's extrapolate to maybe 90 days, via epidemiologic modeling
    #SIR <- function(time, state, parameters) {
    #  par <- as.list(c(state, parameters))
    #  with(par, {
    #    dS <- -beta/N * I * S - epsilon*S
    #    dE <- beta/N * I * S - gamma * E - epsilon*E
    #    dI <- gamma * E - sigma * I - epsilon*I
    #    dR <- sigma * I - epsilon*R
    #    list(c(dS, dE, dI, dR))
    #  })
    #}
    # S -> (beta) E -> (gamma) I -> (sigma) R (we consider nul birth rate delta)
    
    # Change in Susceptible
    DS  <- function(t) -P*Sr*S[t]*C[t]/(S[t]+C[t]+1)
    
    # Change in Contagious, a.k.a Exposed
    DC  <- function(t) P*Sr*S[t]*C[t]/(S[t]+C[t]+1) - 
      min(C[t]*DET, max(beds-I[t]*(1-Dr),0)) - C[t]*Dr - C[t]*Dr
    
    # Change in Infected         
    DI  <- function(t) C[t]*Dr - Mr*I[t] - Mr*I - 
      min(C[t]*DET, max(beds-I[t]*(1-Dr),0))
    
    # Change in Quarantined / instead of Infected (Quarantine more reliable)
    #DI <- function(t) 
    #  min(C[t]*DET, max(beds-I[t]*(1-Dr),0)) -I[t]*Dr
    #DI <- function(t)
    #  min(C[t]*DET, max(beds-I[t]*(1-Dr),0)) -I[t]*Dr - Mr*C[t]
    
    # Change in deceased          
    DD <- function(t) (I[t]+C[t])*Mr*Dr
    #DD <- function(t) (S[t]-I[t]-C[t])*Mr*Dr
    
    # Change in recovered
    DR <- function(t) (I[t]+C[t])*(1-Mr)*Dr
    
    # Change in detection over time
    Et <- function(t) detI+(1-(1-detG)^t)*(detM-detI)
    
    S <- C <- I <- D <- R <- E <- B <- r0 <- rep(NA, np+1)
    
    # Initial populations:
    S[1]  <- N-IC # Sesceptible population
    C[1]  <- IC   # Exposed
    I[1]  <- 0    # Infected
    R[1]  <- 0    # Recovered
    D[1]  <- 0    # Deceased
    B[1]  <- input$bed0
    r0[1] <- P*Days
    
    # Loop through periods
    for (t in 1:np) {
      # Detection rate of unifected per day 
      B[t+1] <- beds <- bedsv[ceiling(t/30)]
      Sr <- (1+K)^(-t)
      r0[t+1] <- P*Sr*Days*(S[t])/(S[t]+C[t])
      
      # Calculate the change values
      dS  <- DS(t) 
      dC  <- DC(t) 
      dI  <- DI(t)
      dR  <- DR(t)
      dD  <- DD(t)
      
      # Change the total populations
      S[t+1] <- S[t] + dS
      C[t+1] <- C[t] + dC
      I[t+1] <- I[t] + dI
      R[t+1] <- R[t] + dR
      D[t+1] <- D[t] + dD
      
    }
    
    # Turn the results into a table
    long <- data.frame(
      Period=rep((0:np),6), 
      Population = c(S, C, D, R, I, B), 
      Indicator=rep(c("Susceptible", 
                      "Exposed", 
                      "Infected", #infected I <-> Deceased D
                      "Recovered",
                      "Deceased",
                      "Beds"), 
                    each=np+1))
    wide <- cbind(S, C, D, R, I, B, r0)
    
    list(long=long, wide=wide)
    
  })
  
  output$r0 <- 
    renderText(paste('Initial r0:', input$P*input$Days))
  output$LDET <- 
    renderText(paste('Likelihood of detection:', round(1-(1-input$DET)^input$Days,2)))
  
  output$datatable <- 
    renderTable({
      Tdata <- mydata()[["wide"]]
      Tdata <- cbind(day=1:nrow(Tdata), Tdata)
      Tdata[seI(1, nrow(Tdata), length.out=20),]
    })
  
  output$graph1 <- renderPlot({
    long <- mydata()[["long"]]
    p <- ggplot(long[long$Indicator %in% input$Indicators,], 
                aes(x=Period, y=Population, group=Indicator))    
    p <- p + 
      geom_line(aes(colour = Indicator), size=1, alpha=.75) + 
      ggtitle("Population Totals")+
      scale_x_continuous(name="Days")+ 
      scale_y_continuous(labels = comma, name="")
    print(p)
  })
  
  output$graph2 <- renderPlot({
    
    data2 <- mydata()[["wide"]]
    
    change <- data2[-1,]-data2[-nrow(data2),]
    
    long <- data.frame(
      Period=rep((1:nrow(change)),7), 
      Population = c(change), 
      Indicator=rep(c("Susceptible", 
                      "Exposed", 
                      "Infected",
                      "Recovered",
                      "Deceased",
                      "Beds",
                      "r0"), 
                    each=nrow(change)))
    
    p <- ggplot(long[long$Indicator %in% input$Indicators,], 
                aes(x=Period, y=Population, group=Indicator))    
    p <- p + geom_line(aes(colour = Indicator), size=1,alpha=.75) + 
      ggtitle("Change in Population")+
      scale_x_continuous(name="Days")+
      scale_y_continuous(labels = comma, name="")
    print(p)
  })
  
}


# Run App
shinyApp(ui = ui, server = server)


#Ali Haidari (alih@vt.edu)
#CMDA 4664 - Computationally Intensive Stochastic Modeling
#4/20/2023

#this is an application script for an Rshiny interactive demo
#that lets the user better understand MCMC methodology.

#the majority of code is adapted from https://shiny.rstudio.com/articles/build.html

#------------------------------------------------
#read in necessary libraries for MCMC
library(FAmle)
library(mcmc)
library(metropolis)
library(shinythemes)

# Define UI for application that draws a panel
ui <- fluidPage(theme = shinytheme("cosmo"),
    tags$style(type="text/css",
         ".shiny-output-error { visibility: hidden; }",
         ".shiny-output-error:before { visibility: hidden; }"),
    # Application title
    strong(titlePanel("Interactive MCMC Rshiny Application")),
    helpText("Ali Haidari (alih) 2023"),
    br(),
    # Sidebar with a slider input for number of bins
    sidebarLayout(
      sidebarPanel(width=2, 
                   helpText("Select from a corpus of datasets with a difficult distribution to sample from."),
                   selectInput("typeInput", "Dataset Selection",
                               choices = c("Star Data", "Poverty Data", "Yarn Data", "Flow Rate Data"),
                               selected = NULL, multiple = F),
                   sliderInput("bins",
                               "Number of bins:",
                               min = 1,
                               max = 50,
                               value = 30)),
      
    # Show a the main panel contents
    mainPanel(
      h3(textOutput("title")),
      h4(textOutput("description")),
      hr(style = "border-top: 1px solid #000000;"),
      tableOutput("table_output"),
      h4("Summary Statisitcs of Dataset"),
      verbatimTextOutput("summary"),
      tags$head(tags$style("#clickGene{color:red; font-size:12px; font-style:italic; 
            overflow-y:scroll; max-height: 50px; background: ghostwhite;}")),
      plotOutput("histogram"),
      h3("MLE estimation for Feature of Interest" ),
      hr(style = "border-top: 1px solid #000000;"),
      verbatimTextOutput("mle_sum"),
      plotOutput("mle"),
      h5("For the method of maximum likelihood estimation, we must know enough about the distribution to choose an close optimization method."),
      hr(style = "border-top: 1px solid #000000;"),
      h3("Posterior Distribution using Metropolis-Hastings"),
      h5("For the Metropolis-Hastings method with no prior distribution and a multivariate normal proposal."),
      h5("The 'metropolis' function from the 'FAmle' package was used to calculate posterior." ),
      plotOutput("metrop_mod"),
      hr(style = "border-top: 1px solid #000000;"),
      h3("Alternative Methodology using Metropolis-Hastings"),
      h5("This is altenative model methodology that calls the 'metropolis_glm()' function from the 'metropolis' package. 
         The difference is thatuses an 'adaptive' proposal distribution along with a 'guided' metropolis algorithm which 
         usually results in better distribution. This is the alternative to the standard 'random-walk' Metropolis-hastings 
         which only accepts the proposal based that minimizes the acceptance probabliity."),
      plotOutput("alt_metrop")
      
      
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  #two datasets from packages
  data("yarns")
  data("station01AJ010")
  yarn_dat <- as.data.frame(yarns)
  flow_dat <- as.data.frame(station01AJ010)
  #data preprocessing
  colnames(flow_dat) <- 'Flow_Rate'
  colnames(yarn_dat) <- 'Failure_Time'
  #data descriptions for all four datasets
  yarn_desc <- "Consists of 100 cycles-to-failure times for airplane yarns."
  flow_desc <- "Consists of 34 realized values of annual maximum daily mean flows in meters cubed per second."
  star_desc <- "Consists of 240 observations of six features of stars across the universe."
  pov_desc <- "Consists of 102 observations of eight features indicating poverty measures."
  
  #two supplementary datasets
  star_dat <- read.csv("star_data.csv")
  pov_dat <- read.csv("poverty_data.csv")
  
  #use reactive function to change dataset and description based off input
  chosen_dat <- reactive({switch(input$typeInput,
                                 "Star Data" = star_dat,
                                 "Poverty Data" = pov_dat,
                                 "Yarn Data" = yarn_dat,
                                 "Flow Rate Data" = flow_dat)})
  chosen_desc <- reactive({switch(input$typeInput,
                                  "Star Data" = star_desc,
                                  "Poverty Data" = pov_desc,
                                  "Yarn Data" = yarn_desc,
                                  "Flow Rate Data" = flow_desc)})
  
  #title chunk
  output$title <- renderText({paste("First Five Observations of the Dataset", 
                                    input$inputType)})
  #summary chunk for data
  output$summary <- renderPrint({(summary(chosen_dat()))})
  #head of dataframe chunk
  output$table_output <- renderTable(head(chosen_dat()))
  #description of dataset chunk
  output$description <- chosen_desc
  
  #mle plots chunk
  output$mle <- renderPlot({
    if (identical(chosen_dat(), pov_dat)){
      x <- chosen_dat()[, 5]
      fit.x <- mle(x,'logis',c(.1,.1))
    }
    else if (identical(chosen_dat(), star_dat)){
      x <- chosen_dat()[, 1]
      fit.x <- mle(x, 'weibull', c(.1,.1))
    }
    else{ 
      x <- chosen_dat()[, 1]
      fit.x <- mle(x,'logis',c(.1,.1))
    }
    
    plot(fit.x)
    plot(fit.x,TRUE,alpha=.01)
    p <- c(.9,.95,.99)
    distr(p,model=fit.x,type='q')
    Q.conf.int(p,fit.x,.01)
    Q.conf.int(p,fit.x,.01,T)
    
  })
  #mle summary chunk
  output$mle_sum <- renderPrint({
    if (identical(chosen_dat(), pov_dat)){
      x <- chosen_dat()[, 5]
      fit.x <- mle(x,'logis',c(.1,.1))
    }
    else if (identical(chosen_dat(), star_dat)){
      x <- chosen_dat()[, 1]
      fit.x <- mle(x, 'weibull', c(.1,.1))
    }
    else{ 
      x <- chosen_dat()[, 1]
      fit.x <- mle(x,'logis',c(.1,.1))
    }
    (fit.x)})
  
  #histogram chunk
  output$histogram <- renderPlot({
    if (identical(chosen_dat(), pov_dat)){
      x <- chosen_dat()[, 5]
      hist_main <- 'Distribution of Intensity of Deprivation Urban'
    }
    
    else if (identical(chosen_dat(), star_dat)){
      x <- chosen_dat()[, 1]
      hist_main <- 'Distribution of Temperature of Star (Kelvin)'
    }
    else{ 
      x <- chosen_dat()[, 1]
      if (identical(chosen_dat(), yarn_dat)){
        hist_main <- 'Distribution of Yarn Failure Times (seconds)'
      }
      else if (identical(chosen_dat(), flow_dat)){
        hist_main <- 'Distribution of Annual Maximum Daily Mean Flow Rate (Meters^3/s)'
      }
    }
    
    #generate bins based on input$bins from ui
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    #draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 4, border = 'white', main = hist_main)
  })
  
  #metropolis model chunk
  output$metrop_mod <- renderPlot({
    if (identical(chosen_dat(), pov_dat)){
      x <- chosen_dat()[, 5]
      fit.x <- mle(x,'logis',c(.1,.1))
      bayes.x.no.prior <- metropolis(model=fit.x, iter=500)
    }
    else if (identical(chosen_dat(), star_dat)){
      x <- chosen_dat()[, 1]
      fit.x <- mle(x, 'weibull', c(.1,.1))
      bayes.x.no.prior <- metropolis(model=fit.x, iter=500, trans.list = c(function(x) x, function(x) exp(x)))
    }
    else{
      x <- chosen_dat()[, 1]
      fit.x <- mle(x,'gamma',c(.1,.1))
      bayes.x.no.prior <- metropolis(model=fit.x, iter=500)
    }
   
    #bayes.x.no.prior$sims
    par(mfrow=c(1,2))
    plot(bayes.x.no.prior, "post.pred")
    hist(x, freq = F)
    lines(density(x), col = 4)
  })
  
  #alternative methods for metropolis-hastings sampling
  output$alt_metrop <- renderPlot({
    if (identical(chosen_dat(), pov_dat)){
      x <- chosen_dat()[, 5]
      feature_form<- as.formula(paste0(colnames(chosen_dat())[5], " ~ 1"))
    }
    else{
      x <- chosen_dat()[, 1]
      feature_form <- as.formula(paste0(sub(' ', '-', colnames(chosen_dat())[1]), " ~ 1"))
    }
    res.rw = metropolis_glm(f =  feature_form , data=chosen_dat(), family=gaussian(), iter=5000, 
                          burnin=1000, adapt=FALSE, guided=FALSE, block=TRUE, inits="glm")
    res.ga = metropolis_glm(f = feature_form, data=chosen_dat(), family=gaussian(), iter=5000, 
                          burnin=1000, adapt=TRUE, guided=TRUE, block=TRUE, inits="glm")
    par(mfrow=c(1,2))
    plot(res.ga, keepburn=TRUE)
    
  })
}
# Run the application (Last line of code)
shinyApp(ui = ui, server = server)
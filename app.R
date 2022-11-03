#
# Virulence evolution - check it out
library("shiny")
library("tidyverse")

b1 <- 1
b2 <- 0.5
c1 <- 1
c2 <- 1

chi <- function(alpha)
{
  return(c1 * alpha^(-c2))  
}

beta <- function(alpha)
{
    return(b1 * alpha^b2)    
}

R0 <- function(alpha, beta, delta, sigma, chi, x, y, h, r)
{
  alphapr <- (1.0 - r[[2]]) * (1.0 - r[[4]]) * alpha
  betapr <- (1.0 - r[[3]]) * beta((1.0 - r[[2]]) * alpha)
#  chipr <- chi((1.0 - r[[2]])*alpha)
  chipr <- chi(alpha)
  hpr <- (1.0 - r[[1]]) * betapr * y
  
  return(
       betapr * (x + sigma * y) /
      (
          delta + alphapr + chipr + sigma * hpr
      )
  )
} #   R0 







# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Virulence evolution in response to medicine"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("r1",
                        div(HTML("Anti-(super)infection immunity: <em>r<sub>1</sub></em>")),
                        min = 0,
                        max = 1,
                        value = 0),
            sliderInput("r2",
                        div(HTML("Anti-growth-rate immunity: <em>r<sub>2</sub></em>")),
                        min = 0,
                        max = 1,
                        value = 0),
            sliderInput("r3",
                        div(HTML("Transmission-blocking immunity: <em>r<sub>3</sub></em>")),
                        min = 0,
                        max = 1,
                        value = 0),
            sliderInput("r4",
                        div(HTML("Anti-toxin immunity: <em>r<sub>4</sub></em>")),
                        min = 0,
                        max = 1,
                        value = 0),
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  virulence.R0.data <- tibble(
    alpha = seq(0,10,0.01)
  )
  
  output$distPlot <- renderPlot({
    
  r <- c(input$r1
           ,input$r2
           ,input$r3
           ,input$r4)
  
      
  virulence.R0.data$R0 <- R0(alpha=virulence.R0.data$alpha, 
                            beta=1, 
                            delta=1, 
                            sigma=1, 
                            chi=1, 
                            x=10, 
                            y=10, 
                            h=10, 
                            r=r)
    
    ggplot(
      data=virulence.R0.data
      ,mapping=aes(x=alpha
                   ,y=R0)) +
      geom_line(colour="#930078",lwd=2) +
      theme_classic() +
      ylim(0,10)
      
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

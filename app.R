library("shiny")
library("tidyverse")
library("ggtext")
library("SIRImperfectVaccination")


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Virulence evolution in response to medical interventions"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          radioButtons(inputId="which_ri"
                       ,label=div(HTML("Which immunity should we vary on <em>x</em> axis?"))
                       ,c(
                         "r1" = 1
                         ,"r2" = 2
                         ,"r3" = 3
                         ,"r4" = 4
                          )
                       ),
            sliderInput("r1",
                        div(HTML("Anti-(super)infection immunity: <em>r<sub>1</sub></em>")),
                        min = 0,
                        max = 1,
                        step=0.05,
                        value = 0),
            sliderInput("r2",
                        div(HTML("Anti-growth-rate immunity: <em>r<sub>2</sub></em>")),
                        min = 0,
                        max = 1,
                        step=0.05,
                        value = 0),
            sliderInput("r3",
                        div(HTML("Transmission-blocking immunity: <em>r<sub>3</sub></em>")),
                        min = 0,
                        max = 1,
                        step=0.05,
                        value = 0),
            sliderInput("r4",
                        div(HTML("Anti-toxin immunity: <em>r<sub>4</sub></em>")),
                        min = 0,
                        max = 1,
                        step=0.05,
                        value = 0),
            sliderInput("f",
                        div(HTML("Fraction vaccinated: <em>f</em>")),
                        min = 0,
                        max = 1,
                        step=0.05,
                        value = 0.2),
            sliderInput("sigma",
                        div(HTML("Probability of superinfection: <em>&#963;</em>")),
                        min = 0,
                        max = 5,
                        step=0.25,
                        value = 1),
        ),

        # Show a plot of the generated distribution
        mainPanel=mainPanel(
          fluidRow(
            column(12,
           plotOutput("virulencePlot"),
           plotOutput("prevalencePlot")
            )
          ) # end fluidrow
        ) # end mainpanel
    ) #  end SideBarLayout
) # end fluidPage

# Define server logic required to draw the resulting plot
server <- function(input, output) {
  
    virulence_data <- reactive({
      
      # make a sequence as we need to vary at least one of the
      # ri's continuously between 0 and 1
      r_i_data <- seq(0,0.9,0.1)
      alpha_data <- rep(NA,length.out=length(r_i_data))
      prevalence_data <- rep(NA,length.out=length(r_i_data))
      
      f <- input$f
      
      sigma <- input$sigma
      
      which_r <- input$which_ri
      
      for (r_i_idx in 1:length(r_i_data))
      {
        r1 = input$r1
        
        if (which_r == 1)
        {
          r1 = r_i_data[[r_i_idx]]
        }
        
        r2 = input$r2
        
        if (which_r == 2)
        {
          r2 = r_i_data[[r_i_idx]]
        }
        
        r3 = input$r3
        
        if (which_r == 3)
        {
          r3 = r_i_data[[r_i_idx]]
        }
        
        r4 = input$r4
        
        if (which_r == 4)
        {
          r4 = r_i_data[[r_i_idx]]
        }
        
        sol <- SIRImperfectVaccination::SIRsolver(
          f = f
          ,delta = 1
          ,sigma=sigma
          ,alpha_init = 1.75
          ,r1 = r1
          ,r2 = r2
          ,r3 = r3
          ,r4 = r4
        )
        
        alpha_data[[r_i_idx]] <- sol$alpha
        prevalence_data[[r_i_idx]] <- (sol$y + sol$yprime)/(sol$y + 
                                         sol$yprime + sol$x + sol$xprime)
      } # end for()
      
      return(data.frame(r=r_i_data,alpha=alpha_data,prevalence=prevalence_data))
    })
    
  
  # make the virulence plot
  output$virulencePlot <- renderPlot({
    
    the_data <- virulence_data()
               
    ggplot(
      data=the_data
      ,mapping=aes(x=r
                   ,y=alpha)) +
      geom_line(colour="#930078",lwd=2) +
      theme_classic(base_size = 16) +
      labs(y="Virulence, *&alpha;*" 
            ,x="Vaccine efficacy, *r*") +
      theme(axis.title.x = ggtext::element_markdown()
            ,axis.title.y = ggtext::element_markdown()
            )
    })
  
  output$prevalencePlot <- renderPlot({
    
    the_data <- virulence_data()
               
    ggplot(
      data=the_data
      ,mapping=aes(x=r
                   ,y=prevalence)) +
      geom_line(colour="#568426",lwd=2) +
      theme_classic(base_size = 16) +
      labs(y="Pathogen prevalence" 
            ,x="Vaccine efficacy, *r*") +
      theme(axis.title.x = ggtext::element_markdown()) +
      ylim(0,1)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

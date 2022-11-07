library("shiny")
library("tidyverse")
library("shinydashboard")
library("ggtext")
library("shinyWidgets")
library("SIRImperfectVaccination")
beta <- function(v, intercept, slope)
{
    return(intercept + v^slope)
}

dbeta <- function(v, intercept, slope)
{
    return(slope * v^(slope - 1.0))
}

gamma <- function(v,intercept, slope)
{
    return(intercept + v^slope)
}

dgamma <- function(v, intercept, slope)
{
    return(slope * v^(slope - 1.0))
}

# expression for R0:
R0_tradeoff <- function(x, beta_intercept, gamma_intercept, beta_slope, gamma_slope, d)
{
    return(1.0 / (x + d + gamma(v=x, intercept=gamma_intercept, slope=gamma_slope))^2 * (dbeta(v=x, intercept=beta_intercept, slope=beta_slope) - beta(v=x, intercept=beta_intercept, slope=beta_slope) / (d + x + gamma(v = x, intercept=gamma_intercept, slope=gamma_slope)) * (1 + dgamma(v = x, intercept=gamma_intercept, slope=gamma_slope)) ))
} # end R0 tradeoff


# Define UI for application that draws a histogram
ui <- dashboardPage(

    # Application title
    dashboardHeader(title="Theory of virulence"),
    dashboardSidebar(
      sidebarMenu(
      menuItem("Trade-offs", tabName="tradeoffs", icon=icon("flash",lib="glyphicon"),selected=T)
      ,menuItem("Interventions", tabName="interventions", icon=icon("plus-sign",lib="glyphicon"))
      ,menuItem("Dynamics", tabName="dynamics", icon=icon("random",lib="glyphicon"))
      )
    )
    ,dashboardBody(
      # first the tab item with all the tradeoff stuff
      tabItems(
        tabItem(
          tabName = "tradeoffs"
          ,fluidRow(
            h2("Potential shapes of trade-offs that involve virulence"),
       
              box( 
               plotOutput("beta_plot"),
                sliderInput("beta_slope",
                            div(HTML("Curvature of transmission - mortality trade-off")),
                            min = -1.5,
                            max = 1.5,
                            value = 0.5),
                   ),
              box(
            plotOutput("gamma_plot"),
            sliderInput("gamma_slope",
                        div(HTML("Curvature of transmission - recovery trade-off")),
                        min = -1.5,
                        max = 1.5,
                        value = -0.5),
               ), # end box 
        
            switchInput(
                inputId="optimize"
                ,value=F
                ,label="Mess with virulence evolution just yet?"
                ,onLabel="Enjoy 🤞!"
                ,offLabel="No thanks!"),
            plotOutput("transmission_MVT")
        
            ), # end fluidrow
            
          ) # end tabItem tradeoffs 
        
      ######## INTERVENTIONS 
      ,tabItem(tabName="interventions"
               ,fluidRow(
                 h2("The evolution of virulence in response to different medical interventions")
                 ,box(
            
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
                 ), # end box
              box(
                 plotOutput("virulencePlot"),
                 plotOutput("prevalencePlot")
              )
               ) # end fluidRow
        ), # end interventions tabitem
        tabItem(tabName = "dynamics"
          ,fluidRow(
            h2("Ecological and evolutionary dynamics of a disease"),
              box(
              sliderInput(inputId = "dynamic_delta"
                          ,label = div(HTML("Mortality rate, <em>&#948;</em>"))
                          ,min = -10
                          ,max = 10
                          ,step=1
                          ,value = 1),
            sliderInput(inputId = "dynamic_lambda"
                        ,label = div(HTML("Influx of susceptibles, <em>&#955;</em>"))
                        ,min = -10
                        ,max = 30
                        ,step=1
                        ,value = 1),
            sliderInput(inputId = "variance_virulence"
                        ,label = div(HTML("Genetic variation in virulence <em>G<sub>v</sub></em>"))
                        ,min = 0
                        ,max = 0.5
                        ,step=0.01
                        ,value = 0.1)
              ),  # end box
            
              box( 
               plotOutput("plot_dynamic"),
               plotOutput("plot_dynamic_log"),
               plotOutput("alpha_dynamic")
              ) # end box
          ) # end fluidrow
        ) # end tabitem
      ) # end tabitems
    ) # end dashboard body
) # end dashboardpage

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
  
    output$beta_plot <- renderPlot({
        
        beta_slope <- input$beta_slope
        
        beta_data <- tibble(v=seq(0,10,0.01))
        
        intercept <- 0
        
        beta_data <- beta_data %>% mutate(
            beta_output = beta(v, intercept, beta_slope)
        )

        ggplot(data=beta_data
               ,mapping=aes(x=v,y=beta_output)) +
              geom_line(colour="#006373",lwd=1.5) +
              theme_classic(base_size=15) +
                xlab("Virulence, v") +
            ylab("Transmission rate") +
            labs(title = "Virulence - transmission trade-off")
    })
    
    output$gamma_plot <- renderPlot({
        
        gamma_slope <- input$gamma_slope
        
        gamma_data <- tibble(v=seq(0,10,0.01))
        
        intercept <- 0
        
        gamma_data <- gamma_data %>% mutate(
            gamma_output = gamma(v, intercept, gamma_slope)
        )

        ggplot(data=gamma_data
               ,mapping=aes(x=v,y=gamma_output)) +
              geom_line(colour="#963d00",lwd=1.5) +
              theme_classic(base_size= 15) +
                xlab("Virulence, v") +
                ylab("Host recovery rate") + 
                labs(title = "Virulence - recovery trade-off")
    })
   
    # plot of tranmission relative to host mortality rate
    output$transmission_MVT <- renderPlot({
        
        vmax <- 10
        # parameters used to optimize
        gamma_slope <- input$gamma_slope
        beta_slope <- input$beta_slope
        beta_intercept <- 0
        gamma_intercept <- 0
        
        # background mortality rate, we might set that as a value
        # but not right now
        d <- 1
        
        mort_transmission <- tibble(host_mortality=seq(-5,vmax,0.01))
        
        mort_transmission <- mutate(mort_transmission,
            transmission=ifelse(host_mortality<0
                                ,NA
                                ,NA
                                )
        )
        
        v_optimum <- NA
        
        # only run this bit if we actually want to optimize
        if (input$optimize)
        {
            
            mort_transmission <- mutate(mort_transmission,
                transmission=ifelse(host_mortality<0
                                    ,NA
                                    ,beta(host_mortality,beta_intercept,beta_slope))
            )
            
            # plot the optimal line
            v_optimum <- optimize(f=R0_tradeoff
                                  ,maximum=F
                                  ,beta_intercept=beta_intercept
                                  ,gamma_intercept=gamma_intercept
                                  ,beta_slope=beta_slope
                                  ,gamma_slope=gamma_slope
                                  ,d=d
                                  ,interval=c(0,vmax))$minimum
         
        } # end if optimize
        
        gplot <- ggplot(data=mort_transmission   
               ,mapping=aes(x=host_mortality,y=transmission)) +
              geom_line(colour="blue",lwd=1.5) +
              theme_classic(base_size= 15) +
                xlab("Host mortality rate, 𝛼") +
                ylab("Transmission rate") + 
                labs(title = "Optimal virulence")
        
        # add the point
        if (!is.na(v_optimum))
        {
            # calculate corresponding y value
            yval = y=beta(v=v_optimum
                             ,intercept=beta_intercept
                             ,slope=beta_slope)
            
            recovery = gamma(v=v_optimum
                             ,intercept=gamma_intercept
                             ,slope=gamma_slope)
            
            # calculate where line crosses x axis 
            # see Gandon et al 2001 Evolution
            x0 <- -d -recovery 
            y0 <- -dgamma(v_optimum, gamma_intercept, gamma_slope)
            
            gplot <- gplot + geom_point(mapping = aes(x=v_optimum,
                                                      y=yval)) +
                            geom_segment(mapping=aes(x=x0
                                                  ,xend=v_optimum
                                                  ,y=y0
                                                  ,yend=yval)
                                         ,linetype="dashed")
        }
        
        return(gplot)
    })
  
    # the SIR model on tab 3
    dynamicSIRdata <- reactive({
        delta = input$dynamic_delta;
        lambda = input$dynamic_lambda
        var_virulence = input$dynamic_lambda
        
        f = 0.3
        n_i_t0 <- 5
        n_s_t0 <- 100
        
        sol <- SIRImperfectVaccination::SIRsolver(
            n_susceptible_init = n_s_t0
            ,n_infected_init = n_i_t0
            ,f = f
            ,delta = delta
            ,sigma=1
            ,lambda=lambda
            ,alpha_init = 1.75
            ,r1 = 0
            ,r2 = 0
            ,r3 = 0
            ,r4 = 0
            ,maxt_eco = 1
            ,maxt_evo = 5000
            ,all_data=T)
        
        sol_l <- pivot_longer(sol
                              ,cols=c(alpha,x,xprime,y,yprime)
                              ,values_to="dynamic_value"
                              ,names_to="dynamic_variable")
        
        return(sol_l)
    }) # end dynamicSIRdata
  
  output$plot_dynamic <- renderPlot({
    
    # get the data from the SIR model
    the_data <- dynamicSIRdata()
    
    # remove alpha, we'll plot it later
    the_data_f <- filter(the_data,
                         !(as.character(dynamic_variable) %in% c("alpha")))
    
    print(sort(unique(the_data$dynamic_variable)))
    
    ggplot(data=the_data_f
           ,mapping=aes(x=t_converge
                        ,y=dynamic_value)) +
      geom_line(mapping=aes(colour=dynamic_variable)) +
      theme_classic(base_size= 15) +
        xlab("Time step, t") +
        ylab("Density") + 
        labs(title = "Dynamics over time") +
        scale_colour_manual(name=""
                            ,values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
                            ,breaks=c("x","xprime","y","yprime")
                            ,labels=c("S","S, vaccinated", "I", "I, vaccinated"))
      
  }) # end output$plot_dynamic
    
  output$plot_dynamic_log <- renderPlot({
    
    # get the data from the SIR model
    the_data <- dynamicSIRdata()
    # remove alpha, we'll plot it later
    the_data_f <- filter(the_data,
                         !(as.character(dynamic_variable) %in% c("alpha")))
    
    ggplot(data=the_data_f
           ,mapping=aes(x=t_converge
                        ,y=log(dynamic_value))) +
      geom_line(mapping=aes(colour=dynamic_variable)) +
      theme_classic(base_size= 15) +
        xlab("Time step, t") +
        ylab("Log density") + 
        scale_colour_manual(name=""
                            ,values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
                            ,breaks=c("x","xprime","y","yprime")
                            ,labels=c("S","S, vaccinated", "I", "I, vaccinated"))
  }) # end output$plot_dynamic
  
  
  output$alpha_dynamic <- renderPlot({
    
    # get the data from the SIR model
    the_data <- dynamicSIRdata()
    
    # get alpha
    the_data_f <- filter(the_data,
                         as.character(dynamic_variable) %in% c("alpha")
                         )
    
    ggplot(data=the_data_f
           ,mapping=aes(x=t_converge
                        ,y=dynamic_value)) +
      geom_line(mapping=aes(colour=dynamic_variable),show.legend=F) +
      theme_classic(base_size= 15) +
        xlab("Time step, t") +
        ylab("Virulence") 
  }) # end output$plot_dynamic
  
}

# Run the application 
shinyApp(ui = ui, server = server)

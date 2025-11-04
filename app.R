library("Rcpp")
library("shiny")
library("tidyverse")
library("shinydashboard")
library("shinyjs")
library("ggtext")
library("shinyWidgets")

if (!require("SIRImperfectVaccination")) {
    devtools::install_github("bramkuijper/SIRImperfectVaccination")
}

beta_trans <- function(v, intercept, slope)
{
    return(intercept + v^slope)
}

dbeta_trans <- function(v, intercept, slope)
{
    return(slope * v^(slope - 1.0))
}

gamma_trans <- function(v,intercept, slope)
{
    return(intercept + v^slope)
}

dgamma_trans <- function(v, intercept, slope)
{
    return(slope * v^(slope - 1.0))
}

# expression for R0:
R0_tradeoff <- function(x, beta_intercept, gamma_intercept, beta_slope, gamma_slope, d)
{
    return(beta_trans(v=x, intercept=beta_intercept, slope = beta_slope) / 
             (x + d + gamma_trans(v=x, intercept=gamma_intercept, slope=gamma_slope))) #0 * (dbeta(v=x, intercept=beta_intercept, slope=beta_slope) - beta(v=x, intercept=beta_intercept, slope=beta_slope) / (d + x + gamma(v = x, intercept=gamma_intercept, slope=gamma_slope)) * (1 + dgamma(v = x, intercept=gamma_intercept, slope=gamma_slope)) ))
    #return(1.0 / (x + d + gamma(v=x, intercept=gamma_intercept, slope=gamma_slope))^2 * (dbeta(v=x, intercept=beta_intercept, slope=beta_slope) - beta(v=x, intercept=beta_intercept, slope=beta_slope) / (d + x + gamma(v = x, intercept=gamma_intercept, slope=gamma_slope)) * (1 + dgamma(v = x, intercept=gamma_intercept, slope=gamma_slope)) ))
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
                  useShinyjs(),
                  h2("The evolution of virulence in response to different medical interventions")
                 ,box(
            
                   radioButtons(inputId="which_ri"
                         ,label=div(HTML("Which immunity should we vary on <em>x</em> axis?"))
                         ,c(
                           "r1" = 1
                           ,"r2" = 2
                           ,"r3" = 3
                           ,"r4" = 4
                            ),
                         selected=1
                         ),
              conditionalPanel(condition = "input.which_ri != '1'",
                sliderInput("r1",
                          div(HTML("Anti-(super)infection immunity: <em>r<sub>1</sub></em>")),
                          min = 0,
                          max = 1,
                          step=0.05,
                          value = 0)
                          ),
              conditionalPanel(condition = "input.which_ri != '2'",
                sliderInput("r2",
                            div(HTML("Anti-growth-rate immunity: <em>r<sub>2</sub></em>")),
                            min = 0,
                            max = 1,
                            step=0.05,
                            value = 0)
                ),
              conditionalPanel(condition = "input.which_ri != '3'",
                sliderInput("r3",
                            div(HTML("Transmission-blocking immunity: <em>r<sub>3</sub></em>")),
                            min = 0,
                            max = 1,
                            step=0.05,
                            value = 0)
              ),
              conditionalPanel(condition = "input.which_ri != '4'",
                sliderInput("r4",
                            div(HTML("Anti-toxin immunity: <em>r<sub>4</sub></em>")),
                            min = 0,
                            max = 1,
                            step=0.05,
                            value = 0)
              ),
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
            sliderInput(inputId = "dynamic_variance_virulence"
                        ,label = div(HTML("Genetic variation in virulence <em>G<sub>&#945;</sub></em>"))
                        ,min = 0
                        ,max = 0.5
                        ,step=0.01
                        ,value = 0.1),
            sliderInput(inputId = "dynamic_f"
                        ,label = div(HTML("Fraction vaccinated <em>f</em>"))
                        ,min = 0
                        ,max = 1.0
                        ,step=0.05
                        ,value = 0.2)
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
      
      which_r <- as.numeric(input$which_ri)
      
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
    }) # end virulence_data()
    
  
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
    }) #virulence vs vaccine efficacy
  
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
    }) # disease prevalence plot
  
    # plot the transmission rate versus v
    output$beta_plot <- renderPlot({
        
        beta_slope <- input$beta_slope
        
        beta_data <- tibble(v=seq(0,10,0.01))
        
        intercept <- 0
        
        beta_data <- beta_data %>% mutate(
            beta_output = beta_trans(v, intercept, beta_slope)
        )

        ggplot(data=beta_data
               ,mapping=aes(x=v,y=beta_output)) +
              geom_line(colour="#006373",lwd=1.5) +
              theme_classic(base_size=15) +
                xlab("Virulence, v") +
            ylab("Transmission rate") +
            labs(title = "Virulence - transmission trade-off")
    }) # end trade-off plot beta
    
    output$gamma_plot <- renderPlot({
        
        gamma_slope <- input$gamma_slope
        
        gamma_data <- tibble(v=seq(0,10,0.01))
        
        intercept <- 0
        
        gamma_data <- gamma_data %>% mutate(
            gamma_output = gamma_trans(v, intercept, gamma_slope)
        )

        ggplot(data=gamma_data
               ,mapping=aes(x=v,y=gamma_output)) +
              geom_line(colour="#963d00",lwd=1.5) +
              theme_classic(base_size= 15) +
                xlab("Virulence, v") +
                ylab("Host recovery rate") + 
                labs(title = "Virulence - recovery trade-off")
    }) # end trade-off plot gamma
   
    # plot of tranmission relative to host mortality rate
    # UPDATE: we will just plot per-hos transmission rat
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
                                    ,beta_trans(host_mortality,beta_intercept,beta_slope))
            )
            
            # plot the optimal line
            v_optimum <- optimize(f=R0_tradeoff
                                  ,maximum=T
                                  ,beta_intercept=beta_intercept
                                  ,gamma_intercept=gamma_intercept
                                  ,beta_slope=beta_slope
                                  ,gamma_slope=gamma_slope
                                  ,d=d
                                  ,interval=c(0,vmax))$maximum
            
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
            yval = y=beta_trans(v=v_optimum
                             ,intercept=beta_intercept
                             ,slope=beta_slope)
            
            recovery = gamma_trans(v=v_optimum
                             ,intercept=gamma_intercept
                             ,slope=gamma_slope)
            
            # calculate where line crosses x axis 
            # see Gandon et al 2001 Evolution
            x0 <- -d -recovery
            y0 <- 0 
            x1 <- v_optimum
            y1 <- yval
            
            slope = (y1 - y0) / (x1 - x0)
            
            intercept = abs((d + recovery))*slope
            
            # let max y be max_beta 
            yend = beta_trans(v = vmax, intercept = beta_intercept, slope = beta_slope)
            
            # now take inverse to calculate xend
            xend = (yend - intercept) / slope
            
            
            gplot <- gplot + geom_segment(mapping=aes(x=x0
                                                  ,xend=xend
                                                  ,y=y0
                                                  ,yend=yend)
                                         ,linetype="dashed") +
              geom_point(mapping = aes(x=v_optimum,
                                                      y=yval),
                                        pch=21,
                                        size=5,
                                        fill="white",
                                        color="black") 
        }
        
        return(gplot)
    }) # end transmission_MVT (tab 1)
  
    # the SIR model on tab 3
    dynamicSIRdata <- reactive({
        delta = input$dynamic_delta;
        lambda = input$dynamic_lambda
        var_virulence = input$dynamic_variance_virulence
        f = input$dynamic_f
        n_i_t0 <- 5
        n_s_t0 <- 100
       
        # solve over time 
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
            ,eul_evo = var_virulence
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

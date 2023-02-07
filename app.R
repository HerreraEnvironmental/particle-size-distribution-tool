#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(signal)
library(dplyr)
library(rhandsontable)
library(ggplot2)

initial_table<-
  tibble(`Screen Size (µm)`=c(500,125,63,4,1)	,
         `Mass Caught (mg)`=c(11,12,12.7,26.8,6)
  )
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

sieve_lookup<-read.csv('sieve_lookup.csv')

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Particle Size Distribution Tool"),

    # Sidebar with a slider input for number of bins 

          rHandsontableOutput('hot'),
      
    tableOutput("sieveMeshTable"),
    fluidRow(column(4,tableOutput("distrib")),
             column(4,tableOutput("coeffs"))),
    plotOutput('plot',
               width='600px')

    
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  df <- reactive({
    hot <- input$hot
    if (!is.null(hot)) hot_to_r(hot)
  })
  
  output$hot <- renderRHandsontable({
    rhandsontable(initial_table)
  })
  #2. Spreadsheet plots the seive analysis
  sieveMesh<-reactive({
    sieve_mesh<-left_join(df(),sieve_lookup,by=c(`Screen Size (µm)`='Microns')) %>% pull(Mesh)
    df() %>%
    transmute(
      `Screen Size (µm)`,
      `ASTM US Std Mesh`=sieve_mesh,
      `% Retained` = `Mass Caught (mg)`/sum(`Mass Caught (mg)`),
      `Mass Passing (mg)`=sum(`Mass Caught (mg)`)-cumsum(`Mass Caught (mg)`),
      `% Passing`=`Mass Passing (mg)`/sum(`Mass Caught (mg)`)
    )
  })
  
  output$sieveMeshTable<-renderTable({
    sieveMesh()%>%
      mutate(
        `% Retained` = percent(`% Retained`),
        `% Passing`=percent(`% Passing`)
      ) %>%
      select(-`Screen Size (µm)`)
  })
  #3. Interpolate D25, D50, D75, D84, and D90

  output$distrib<-renderTable({
    interp_out<-
      sieveMesh() %>%
      arrange(`% Passing`) %>%
      with(.,signal::interp1(x=`% Passing`,
                                y=log(`Screen Size (µm)`),
                                 xi=c(.1,.25,.3,.5,.6,.75,.87,.9),
                                )) %>%
      exp() %>%
      round(0)
    tibble(Quantile=c('D10','D25','D30','D50','D60','D75','D84','D90'),
           `Particle Size (µm)`= ifelse(is.na(interp_out),
                                        paste0('>',max(sieveMesh()$`Screen Size (µm)`)),
                                        interp_out)
                                        )
  })
  #4. Calculate Coefficient of Uniformity
  output$coeffs<-renderTable({
    interp_out<-
      sieveMesh() %>%
      arrange(`% Passing`) %>%
      with(.,signal::interp1(x=`% Passing`,
                             y=log(`Screen Size (µm)`),
                             xi=c(.1,.3,.6),
      )) %>%
      exp()
    
    tibble(Coeff=c('Coefficient of Uniformity\n(D60/D10)',
                   'Coefficient of Curve\n(D30^2/D10/D60)'),
           Value=c(interp_out[3]/interp_out[1],
                   interp_out[2]^2/interp_out[1]/interp_out[3])
      )
  })
  
  #4. Plot it up
  output$plot<-renderPlot({
    breaks_minor<-expand.grid(c(1,10,100,1000),1:9) %>% mutate(mbs=Var1*Var2) %>% arrange(mbs) %>%
      pull(mbs)
    
    reverselog_trans <- function(base = exp(1)) {
      trans <- function(x) -log(x, base)
      inv <- function(x) base^(-x)
      scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
                scales::log_breaks(base = base), 
                domain = c(1e-100, Inf))
    }
    
    ggplot(sieveMesh(),aes(`Screen Size (µm)`,`% Passing`))+
      geom_rect(ymin=-Inf,ymax=Inf,xmin=Inf,xmax=-log10(2),fill='blue',alpha=.1)+
      geom_rect(ymin=-Inf,ymax=Inf,xmin=-log10(2),xmax=-log10(62),fill='red',alpha=.1)+
      geom_rect(ymin=-Inf,ymax=Inf,xmin=-log10(62),xmax=-log10(2000),fill='green',alpha=.1)+
      geom_rect(ymin=-Inf,ymax=Inf,xmin=-log10(2000),xmax=-log10(64000),fill='grey',alpha=.1)+
      geom_rect(ymin=-Inf,ymax=Inf,xmin=-log10(64000),xmax=-Inf,fill='black',alpha=.1)+
      geom_vline(xintercept=c(1,10,100,1000,10^4))+
      geom_vline(xintercept=c(1,2,62,2000,64000),lty=2)+
      geom_line()+
      geom_point()+
      theme_bw()+
      scale_y_continuous('% Passing',labels=scales::percent,limits=c(0,1))+
      scale_x_continuous(limits=c(10^4,1),
                         breaks=c(1,10,100,1000,10^4),
                       minor_breaks=breaks_minor,
                       trans=reverselog_trans(10))+
      geom_label(data=tibble(`% Passing`=.5, `Screen Size (µm)`=#68
                               sieveMesh() %>%
                   arrange(`% Passing`) %>%
                   with(.,signal::interp1(x=`% Passing`,
                                          y=log(`Screen Size (µm)`),
                                          xi=c(.5))) %>% exp()
                   ),aes(label=paste0('D50: ',round(`Screen Size (µm)`),' µm')))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

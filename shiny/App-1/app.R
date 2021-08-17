library(shiny)
library(shinyBS)

# Define UI for app that draws a histogram ----
# ui <- fluidPage(
#   titlePanel("My Shiny App"),
#   
#   sidebarLayout(
#     sidebarPanel(h2("Installation"),
#                  p("Shiny is available on CRAN, so you can install it in the usual way from your R console:"),
#                  code("install.packages(\"shiny\")")),
#     mainPanel(
#       
#       p("Shiny is a new package from RStudio that makes it ", 
#         em("incredibly ", strong("frickin"), " easy "), 
#         "to build interactive web applications with R."),      br(),
#       p("span does the same thing as div, but it works with",
#         span("groups of words", style = "color:blue"),
#         "that appear inside a paragraph."),
#       img(src = "ArnoldBayes.png", width = 300)
#       
#     ),
#     
#   )
# )


# ui <- fluidPage(
#   titlePanel("censusVis"),
#   
#   sidebarLayout(
#     sidebarPanel(p("Create demographic maps with information from the 2010 US Census."),
#                  br(),
#                  selectInput("", label = "Choose a variable to display", 
#                              choices = c("Percent White", "Percent Black", "Percent Hispanic", "Percent Asian")),
#                  sliderInput("", "Range of Interest",
#                              min = 0, max = 100, value = c(0, 100))),
#     mainPanel(
#       
#       p("Shiny is a new package from RStudio that makes it ", 
#         em("incredibly ", strong("frickin"), " easy "), 
#         "to build interactive web applications with R."),      br(),
#       p("span does the same thing as div, but it works with",
#         span("groups of words", style = "color:blue"),
#         "that appear inside a paragraph."),
#       img(src = "ArnoldBayes.png", width = 300)
#       
#     ),
#     
#   )
# )
# 
library(shiny)

#@import url('https://fonts.googleapis.com/css2?family=Yusei+Magic&display=swap');
#h2 {
# font-family: 'Yusei Magic', sans-serif;
#}

ui <- fluidPage(
  
  tags$head(
    uiOutput("css")
  ), #cursor: url(https://emoji.discord.st/emojis/Reddit.png), auto;}
  
  titlePanel("testing out shiny functionality"),

  sidebarLayout(
    sidebarPanel(width=3,
                 verbatimTextOutput("info"),
                 actionButton("reset", "Clear Circles"),
                 br(), br(),
                 p("Hovering over the top right square -- at (1,1) -- changes the color of all squares to green. 
                   How can I change the cursor to, say, a pointing hand while that happens?")
                 
                 
    ),

    mainPanel(
      uiOutput("plot.ui"),
      tableOutput("pointlist")
    )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize=100*1024^2) # set maximum image size

  xy <- reactiveValues(x= numeric(0), y = numeric(0), col = numeric(0), line=numeric(0)) # add new points
  square_color <- reactiveValues(col = "red") # add new points

  output$css <- renderUI({
    
    if(square_color$col == "green"){
      css_string <- "
        #distplot {  
        cursor: pointer;}"
    } else {
      css_string <- "
      #distplot {  
      cursor: default;}"
    }
    tags$style(HTML(css_string))
  })
  
  output$plot.ui <- renderUI({
    plotOutput("distplot",
               click = "plot_click",
               hover = hoverOpts(
                 "plot_hover",
                 delay = 5,
                 delayType = c("debounce", "throttle"),
                 clip = TRUE,
                 nullOutside = TRUE
               ))
  })
  
  
  # observe hover values and clicks
  observe({
    
    if (is.null(input$plot_hover)){
      return()
    }
    
    if(sqrt((input$plot_hover$x - 1)^2 + (input$plot_hover$y - 1)^2) < 0.125){
      square_color$col <- "green"
    } else {
      square_color$col <- "red"
    }
    
    })
  
  observe({
    
    if (is.null(input$plot_click)){
      return()
    }

    isolate({
      dists <- sqrt((xy$x - input$plot_click$x)^2 + (xy$y - input$plot_click$y)^2)
      pts_to_remove <- which(dists < 0.1)
      if(length(pts_to_remove) > 0){
        xy$x <- xy$x[-pts_to_remove]
        xy$y <- xy$y[-pts_to_remove]
        xy$col <- xy$col[-pts_to_remove]
      } else{
        xy$x <- c(xy$x, input$plot_click$x)
        xy$y <- c(xy$y, input$plot_click$y)
        xy$col <- c(xy$col, sample(1:10, 1))
      }
      
    })
  })
  
  #clear all circles when button is pressed
  observeEvent(input$reset, {

    xy$x= numeric(0) 
    xy$y= numeric(0) 
    xy$col= numeric(0) 
    xy$line= numeric(0) 
    
    }) 

  output$distplot <- renderPlot({

    plot(xy$x, xy$y, xlim=c(-2, 2), ylim=c(-2, 2), xlab="", ylab="", pch = 19, cex = 2, col = xy$col)
    points(x = c(-1,-1,1,1), y = c(-1,1,-1,1), cex = 3, col = square_color$col, pch = 15)
    # n=5E3
    # if(square_color$col == "green"){
    #   segments(rnorm(n),rnorm(n),rnorm(n),rnorm(n))
    # }

  })

  output$pointlist <- renderTable(data.frame(x = xy$x, y = xy$y, col = xy$col))
  
  output$info <- renderText({
    xy_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("x=", round(e$x, 2), " y=", round(e$y, 2), "\n")
    }

    paste0(
      "click: ", xy_str(input$plot_click),
      "hover: ", xy_str(input$plot_hover)
    )
  })
}

# library(shiny)
# library(DT)
# 
# shinyApp(
#   ui = fluidPage(
#     DT::dataTableOutput("irisTable")
#   ),
#   server = function(input, output) {
#     
#     output$irisTable <- DT::renderDataTable({
#       DT::datatable(datasets::iris, 
#                     options = list(rowCallback = JS(
#                       "function(nRow, aData, iDisplayIndex, iDisplayIndexFull) {",
#                       "var full_text = aData[1] + ','+ aData[2] + ',' + aData[3] + ','+ aData[4];",
#                       "$('td:eq(5)', nRow).attr('title', full_text);", # Tool Tip
#                       "$('td:eq(5)', nRow).css('cursor', 'pointer');", # Cursor icon changes to hand (pointer) on Hover
#                       "}")
#                     )
#       )
#       
#     })
#   }
# )
# Run the application 
shinyApp(ui, server)
# # Define server logic required to draw a histogram ----
# server <- function(input, output) {
#   
#   
# }

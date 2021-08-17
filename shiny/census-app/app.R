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
# x = rnorm(100); y = rnorm(100)
# ui <- fluidPage(
#   titlePanel("censusVis"),
#   
#   sidebarLayout(
#     sidebarPanel(
#       helpText("Create demographic maps with 
#                information from the 2010 US Census."),
#       
#       selectInput("var", 
#                   label = "Choose a variable to display",
#                   choices = c("Percent White", 
#                               "Percent Black",
#                               "Percent Hispanic", 
#                               "Percent Asian"),
#                   selected = "Percent White"),
#       
#       selectInput("col", 
#                   label = "Choose a color to use",
#                   choices = c("red", 
#                               "blue",
#                               "green", 
#                               "black"),
#                   selected = "black"),
#       
#       sliderInput("range", 
#                   label = "Range of interest:",
#                   min = -5, max = 5, value = c(0, 1), step = 0.01)
#     ),
#     
#     mainPanel(
#       plotOutput("selected_var")
#     )
#   )
# )
# 
# server <- function(input, output) {
#   output$selected_var <- renderPlot({ 
#     plot(x, y, xlim = input$range, ylim = input$range, main = input$var, col = input$col, pch = 19)
#   })
# }

# # Define server logic required to draw a histogram ----
# server <- function(input, output) {
#   
#   
# }


library(maps)
library(mapproj)
source("helpers.R")
counties <- readRDS("data/counties.rds")

x = rnorm(100); y = rnorm(100)
ui <- fluidPage(
  titlePanel("censusVis"),

  sidebarLayout(
    sidebarPanel(
      helpText("Create demographic maps with
               information from the 2010 US Census."),

      selectInput("var",
                  label = "Choose a variable to display",
                  choices = c("Percent White",
                              "Percent Black",
                              "Percent Hispanic",
                              "Percent Asian"),
                  selected = "Percent White"),

      selectInput("col",
                  label = "Choose a color to use",
                  choices = c("red",
                              "blue",
                              "white",
                              "black"),
                  selected = "black"),

      sliderInput("range",
                  label = "Range of interest:",
                  min = 0, max = 100, value = c(0, 100), step = 0.01)
    ),

    mainPanel(
      plotOutput("map")
    )
  )
)

server <- function(input, output) {
  output$map <- renderPlot({
    data <- switch(input$var, 
                   "Percent White" = counties$white,
                   "Percent Black" = counties$black,
                   "Percent Hispanic" = counties$hispanic,
                   "Percent Asian" = counties$asian)
    percent_map(var = data, input$col, legend.title = input$var, min = input$range[1], max = input$range[2])
  })
}

shinyApp(ui = ui, server = server)

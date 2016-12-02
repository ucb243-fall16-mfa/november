#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Multiple Factor Analysis"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("checkGroup", label = h3("Assessor"),
                         choices = list("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5,
                                        "6" = 6, "7" = 7, "8" = 8, "9" = 9, "10" = 10),
                         selected = c(1:10)),
      radioButtons("radio", label = h3("Plot Type"),
                   choices = list("compromise" = 1, "partial factor scores" = 2, "variable loadings" = 3),
                   selected = 1)
        ),
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$distPlot <- renderPlot({
    if (!is.null(input$checkGroup)) {
      sets=list(2:7,8:13,14:19,20:23,24:30,31:35,36:39,40:45,46:50,51:54)
      selected <- sets[as.numeric(input$checkGroup)]
      mfa1 <- mfa(df, sets=selected)
      if (input$radio == 1) plot.comp(mfa1)
      if (input$radio == 2) plot.pfs(mfa1)
      if (input$radio == 3) plot.loadings(mfa1)
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)


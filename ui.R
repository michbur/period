library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("Periodicity in qPCR data"),
  sidebarPanel(
    fileInput("input.file", "Choose CSV File (input should contain Cq data)",
              accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
    checkboxInput("header", "Header", TRUE),
    radioButtons("csv.type", "Type of csv file",
                 c("Dec: dot (.), Sep: comma (;)" = "csv1",
                   "Dec: comma (.), Sep: semicolon (;)" = "csv2")),
      downloadButton("result.download", "Download report"),
    br(),
    p("Lost? Use button below to see an example:"),
    actionButton("run.example", "Run example")
  ),
  mainPanel(
    uiOutput("dynamic.tabset") 
  )
)
)
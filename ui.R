library(shiny)

shinyUI(pageWithSidebar(
  tags$head(h1("Periodicity in qPCR data"), includeScript("ga.js")),
  #headerPanel("Periodicity in qPCR data"),
  sidebarPanel(
    includeMarkdown("readme.md"),
    p("Lost? Use button below to see an example:"),
    actionButton("run.example", "Run example"),
    br(), br(),
    fileInput("input.file", "Choose CSV File (input should contain Cq data)",
              accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
    checkboxInput("header", "Header", TRUE),
    checkboxInput("transpose", "Cq values are in one row", FALSE),
    radioButtons("csv.type", "Type of csv file",
                 c("Dec: dot (.), Sep: comma (;)" = "csv1",
                   "Dec: comma (,), Sep: semicolon (;)" = "csv2")),
    downloadButton("result.download", "Download report")    
  ),
  mainPanel(
    uiOutput("dynamic.tabset") 
  )
)
)
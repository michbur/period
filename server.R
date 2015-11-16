library(shiny)
source("CheckPeriodApp.R")

setwd("./report_generation/")

# server for the Shiny app
shinyServer(function(input, output) {
  
  #check if no data is loaded or no example used
  null.input <- reactive({
    is.null(input[["input.file"]]) && input[["run.example"]] == 0
  })
  
  processed.data <- reactive({
    #after loading any file it would be possible to start an example
    if(is.null(input[["input.file"]])) {
      dat <- read.csv("StepOneCq.csv")
    } else {
      dat <- switch(input[["csv.type"]], 
                    csv1 = read.csv(input[["input.file"]][["datapath"]], 
                                    header = input[["header"]]),
                    csv2 = read.csv2(input[["input.file"]][["datapath"]], 
                                     header = input[["header"]]))
      if(!input[["header"]])
        colnames(dat) <- paste0("Column", 1L:ncol(dat))
    }
    
    dat
  })
  
  #dabset before and after data input
  output[["dynamic.tabset"]] <- renderUI({
    if(null.input()) {
      tabPanel("No input detected",
               HTML("No input detected. <br> Select input file or example using the left panel."))
    } else {
      tabsetPanel(
#         tabPanel("Results with graphics", 
#                  plotOutput("fit.plot"), htmlOutput("fit.text"),
#                  plotOutput("res.plot"), plotOutput("ac.plot"),
#                  plotOutput("hm.plot")),
        tabPanel("Results with graphics", htmlOutput("whole.report")),
        tabPanel("Input table", tableOutput("input.data"))
      )
    }
  })
  
  res.period <- reactive({
    CheckPeriodApp(processed.data())
  })
  
  
  create.md <- reactive({  
    knitr::knit(input = "period_report.Rmd", 
                output = "period_report.md", quiet = TRUE)
  })
  
  output[["input.data"]] <- renderTable({
    processed.data()
  })
  

  
  output[["whole.report"]] <- renderText({
    create.md()
    markdown::markdownToHTML("period_report.md", output = NULL, fragment.only = TRUE)
  })
  
  
  output[["result.download"]] <- downloadHandler(
    filename  = "period_report.html",
    content = function(file) {
      markdown::markdownToHTML("period_report.md", file)
    }
  )
  
})
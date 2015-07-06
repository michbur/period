library(shiny)
source("CheckPeriodApp.R")

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
               HTML('<p><img src="http://www.dr-spiess.de/ans.gif" width="55%" height="55%"/></p>'))
    } else {
      tabsetPanel(
        tabPanel("Input data", tableOutput("input.data")),
        tabPanel("Results with graphics", plotOutput("fit.plot"), 
                 plotOutput("res.plot"), plotOutput("ac.plot"),
                 plotOutput("hm.plot"))
      )
    }
  })
  
  res.period <- reactive({  
    dat <- processed.data()
    
    res <- CheckPeriodApp(dat)
    res
  })
  
  output[["input.data"]] <- renderTable({
    processed.data()
  })
  
  output[["fit.plot"]] <- renderPlot({
    plotFit(res.period())
  })
  
  output[["res.plot"]] <- renderPlot({
    plotRes(res.period())
  })
  
  output[["ac.plot"]] <- renderPlot({
    plotAc(res.period())
  })
  
  output[["hm.plot"]] <- renderPlot({
    plotHm(res.period())
  })
  

  
  output[["result.download"]] <- downloadHandler(
    filename  = "period_report.html",
    content <- function(file) {
      knitr:::knit(input = "period_report.Rmd", 
                   output = "period_report.md", quiet = TRUE)
      markdown:::markdownToHTML("period_report.md", file)
    }
  )
  
})
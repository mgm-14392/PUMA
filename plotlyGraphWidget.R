#but the user should not have to edit it

#Output Graph Function
graphOutput <- function(inputId, width="100%", height="550px") {
  tagList(
    singleton(tags$head(
      tags$script(src="plotlyGraphWidget.js")
    )),
    tags$iframe(id=inputId, src="https://plot.ly/~playground/7.embed",
                class="graphs", style="border:none;", seamless=TRUE, width=width, height=height)
  )
}

renderGraph <- function(expr, env=parent.frame(), quoted=FALSE) {
  ## This gets called when inputs change --
  ## Place data wrangling code in here
  ## and pass the result to the client
  ## to be graphed.
  
  
  installExprFunction(expr, "func", env, quoted)
  
  function(){
    data = func();
    ## data is the state of the widgets: see server.R
    ## this function returns a list of named lists that descripe
    ## valid postMessage commands to be sent to the embedded iframe.
    ## see binding.renderValue for the receiving JS side of this function
    ## and https://github.com/plotly/Embed-API for more about the postMessage
    ## graph messages
    return(data)
    
  }
}
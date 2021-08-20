# mushroom_test1.R - quick initial implementation of the fairy ring cellular automaton

## Convert to Rshiny
#
#

# Jayson Rook
#
# Forked by Jacob Veta
# March 9, 2020
#
# Changes:
#  - converted matrices to single time series array and applied gganimate
#  - converted if/else structure to case_when structure for speed
# 
#
#
# Using plotrix for now, needs replaced badly (but can use install.packages( "plotrix" ) to use for now)
library(plotrix)
library(tidyverse)
library(abind)
library(shiny)
library(animation)

iter=20
probSporeToHyphae = 0.3
probMushroom = 0.7
probSpore = 0.2

ui <- fluidPage(
  # Application title
  titlePanel(title = "Fairy Mushroom Explorer"),
  sidebarLayout(
    # Sidebar typically used to house input controls
    sidebarPanel(
      tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 500px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #CCFF66;
               z-index: 105;
             }
          ")),
      sliderInput(inputId="animation",
                  label="Iterations: ",
                  value=1, min=1, max=iter, step = 1,
                  animate=animationOptions(interval=500,loop=T)
      ),
      sliderInput(inputId="n",
                  label="Grid size: ",
                  value=10, min=4, max= 100, step = 1,
      ),
      sliderInput(inputId="p",
                  label="Probability of Spread",
                  value=0.5, min=0, max=2, step = 0.01,
      ),
      sliderInput(inputId="d",
                  label="Diagonal Spread Factor",
                  value=0.5, min=0, max=1, step = 0.01,
      ),
      checkboxInput("obstructions", "Generate obstructions", value = FALSE),
      actionButton("dl","Generate .gif"),
      #conditionalPanel(condition="$('html').hasClass('shiny-busy')",
       #                tags$div("Releasing Spores...",id="loadmessage"))
    ),
    
    # Main panel typically used to display outputs
    mainPanel(
      tabsetPanel(id="tabs",
                  tabPanel("Unobstructed", plotOutput(outputId="mushplot")),
                  tabPanel("Distance Plot", plotOutput(outputId="distplot"))
      )
    )
  )
)


server <- function(input, output) {
  
  sliderValues <- reactive({(input$animation)})
  
  aniGraph <- reactive({
    # output$mushplot <- renderImage({
    
    ## Tunable probability parameters #find research and convert to 2d matrix for comparison  
    n=input$n
    
    initPointRow = floor(n/2) 
    initPointCol = floor(n/2)
    initSquareSize = 1
    
    ts <- array(numeric(), c(n, n, iter))
    colormat_s <- array(character(), c(n, n, iter))
    colormat_s[,,]=0
    #Adding obstructions
    
    ts[,,]=0
    ##
    # Add parameters for units (x/y ~ in/in ~ m/m, etc)
    #
    #
    
    # Table of colors from page 718
    # Table of colors from page 718
    colorcode = c("lightgreen",
                  "black",
                  "darkgrey",
                  "lightgrey",
                  "white",
                  "lightgrey",
                  "tan",
                  "brown",
                  "darkgreen",
                  "yellow")
    
    
    
    # Initialize cell grid to have custom-sized square at initPointRow, initPointCol
    for (row in 1:n) {
      for (col in 1:n) {
        if (row >= initPointRow - floor(initSquareSize / 2 - 0.1) && row <= initPointRow + floor(initSquareSize / 2 + 0.1)
            && col >= initPointCol - floor(initSquareSize / 2 - 0.1) && col <= initPointCol + floor(initSquareSize / 2 + 0.1)) {
          ts[row,col,1] = 1  
          colormat_s[row,col,1] = colorcode[2]
        } else {
          colormat_s[row,col,1] = colorcode[1]
        }
      }
    }
    
    if(input$obstructions){
      if (row < n * 0.80 && row > n * 0.60 && col < n * 0.10 && col > n * 0.00) { 
        ts[row,col,1] = 9
        colormat_s[row,col,1] = colorcode[10]}
      if (row < n * 0.80 && row > n * 0.60 && col < n * 0.40 && col > n * 0.20) { 
        ts[row,col,1] = 9
        colormat_s[row,col,1] = colorcode[10]}
      if (row < n * 0.80 && row > n * 0.60 && col < n * 0.70 && col > n * 0.50) { 
        ts[row,col,1] = 9
        colormat_s[row,col,1] = colorcode[10]}
      if (row < n * 0.80 && row > n * 0.60 && col < n * 1.00 && col > n * 0.80) { 
        ts[row,col,1] = 9
        colormat_s[row,col,1] = colorcode[10]}
      else {
        colormat_s[row,col,1] = colorcode[1]
      }
    }
    
    # Check if the cell (row, col) has a Young cell in its Moore neighborhood at iteration j
    hasYoungMooreNeighbor <- function(row, col,j) {
      if (row > 1) {
        if (col > 1) {
          if (ts[row-1,col-1,j] == 2) {
            return(TRUE)
          }
        }
        if (col < n) {
          if (ts[row-1,col+1,j] == 2) {
            return(TRUE)
          }
        }
        if (ts[row-1,col,j] == 2) {
          return(TRUE)
        }
      }
      if (row < n) {
        if (col > 1) {
          if (ts[row+1,col-1,j] == 2) {
            return(TRUE)
          }
        }
        if (col < n) {
          if (ts[row+1,col+1,j] == 2) {
            return(TRUE)
          }
        }
        if (ts[row+1,col,j] == 2) {
          return(TRUE)
        }
      }
      if (col > 1) {
        if (ts[row,col-1,j] == 2) {
          return(TRUE)
        }
      }
      if (col < n) {
        if (ts[row,col+1,j] == 2) {
          return(TRUE)
        }
      }
      return(FALSE)
    }
    
    # Check if the cell (row, col) has a Young cell in its Von Neumann neighborhood at iteration j
    hasYoungVonNeumannNeighbor <- function(row, col, j) {
      if (col > 1 && ts[row,col-1,j] == 2) {
        return(TRUE)
      }
      if (row > 1 && ts[row-1,col,j] == 2) {
        return(TRUE)
      }
      if (col < n && ts[row,col+1,j] == 2) {
        return(TRUE)
      }
      if (row < n && ts[row+1,col,j] == 2) {
        return(TRUE)
      }
      return(FALSE)
    }
    
    # Determined whether the fungus should spread to cell (row, col) on iteration j
    spreadToCell <- function(row, col, j) {
      randomValue = runif(1)
      if (hasYoungVonNeumannNeighbor(row, col, j)) {
        # There's a young cell directly adjacent, so use higher probability
        if (randomValue < probSpread) {
          return(TRUE)
        }
      } else if (hasYoungMooreNeighbor(row, col, j) && randomValue < probSpread * input$d) {
        return(TRUE)
      }
      return(FALSE)
    }
    
    
    # ---------------------------------------------------------------------------------------------------------------
    
    
    # Run iter iterations
    for (t in 2:iter) {
      
      # Iterate over every cell
      for (row in 1:n) {
        for (col in 1:n) {
          # Follow different rules based on old state, following state diagram on page 717
          oldstate = ts[row,col,t-1]
          ts[row,col,t]= case_when(oldstate==0 ~ 
                                     if (spreadToCell(row,col,t-1)) {
                                       2 # young
                                     } else {
                                       0 # empty
                                     },
                                   oldstate==1 ~ # cell is previously SPORE
                                     if (runif(1) < probSporeToHyphae) {
                                       2 # young
                                     } else {
                                       1 # spore
                                     },
                                   oldstate==2 ~ 3, # maturing# cell is previously YOUNG
                                   oldstate==3 ~ if (runif(1) < probMushroom) {
                                     4 # mushrooms
                                   } else {
                                     5 # older
                                   },
                                   oldstate==4 ~ 6, # decaying
                                   oldstate==5 ~ 6, # decaying
                                   oldstate==6 ~ 7, # dead1
                                   oldstate==7 ~ 8, # dead2
                                   oldstate==8 ~ 0, # empty
                                   oldstate==9 ~ 9 # inert
          )
          # Update the color matrix to reflect the actual cell state value
          colormat_s[row,col,t] = colorcode[ts[row,col,t]+1]
        }
      
      
    }
      
    }
    newList <- list(colormat_s=colormat_s, ts=ts)
    return(newList)
  })
  
  output$mushplot <- renderPlot({
    n=input$n
    newList = aniGraph()
    colormat_s=newList$colormat_s
    color2D.matplot(matrix(1:(n*n),nrow=n),cellcolors=colormat_s[,,sliderValues()],yrev=FALSE,
                    border="lightgray",xlab="",ylab="",main="Mushroom Growth Over Time",axes=FALSE)
    
  })
  
  output$distplot <- renderPlot({
    n=input$n
    newList = aniGraph()
    ts=newList$ts
    nPoints = 0
    tVals = c()
    rVals = c()
    
    for (t in 1:sliderValues()) {
      for (x in 1:n) {
        for (y in 1:n) {
          if (ts[x,y,t] == 4) {
            nPoints = nPoints + 1
            tVals[nPoints] = t
            rVals[nPoints] = as.numeric(sqrt((x - (n/2))^2 + (y - (n/2))^2))
          }
        }
      }
    }
    
    plot(tVals, rVals,xlab="time (iters)",ylab="distance (cells)",
         main="Mushroom Growth Distances", xlim=c(0,sliderValues()), ylim=c(0,n / sqrt(2)))
  })
  
  observeEvent(input$dl,{
    n=input$n
    newList = aniGraph()
    colormat_s=newList$colormat_s
    saveGIF({
      for(i in 1:iter){
        color2D.matplot(matrix(1:(n*n),nrow=n),cellcolors=colormat_s[,,i],yrev=FALSE,
                        border="lightgray",xlab="",ylab="",main="Mushroom Growth Over Time",axes=FALSE)
      }
    }, movie.name = "mushroom.gif", interval = 0.5, ani.width = 550, ani.height = 350)  #relate this to real time scale with multiplier for simulation
    
  })
}


# Run the application 
runApp(list(ui = ui, server = server))

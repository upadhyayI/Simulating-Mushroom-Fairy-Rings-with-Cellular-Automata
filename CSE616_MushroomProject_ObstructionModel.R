# Michael Larson
# CSE616 Simulation of Physical Systems
# Mushroom Cellular Automota Code -- Obstructions to Growth

# Jayson Rook
#
# Forked by Jacob Veta
#
# Changes:
#  - converted matrices to single time series array and applied gganimate
#  - converted if/else structure to case_when structure for speed

# Using plotrix for now, needs replaced badly (but can use install.packages( "plotrix" ) to use for now)
library(plotrix)
library(tidyverse)
library(abind)

## Tunable probability parameters #find research and convert to 2d matrix for comparison  
probSporeToHyphae = 0.30     #affects spore to young Hyphae
probMushroom = 0.70          #affects maturing Hyphae to mushroom
probSpread = 0.825           #affects empty to young
probSpore = 0.2              #for random initialization
iter = 65

## Graphics Parameters

speed = 0.50  #GIF output transition interval
n = 65  # Side length for square grid
nsq = n*n
color_choice = "complex" # determines whetehr all or few colors are used: choose "simple" or "complex"

ts <- array(numeric(), c(n, n, iter))
colormat_s <- array(character(), c(n, n, iter))
colormat_s[,,]=0
ts[,,]=0

##
# Add parameters for units (x/y ~ in/in ~ m/m, etc)

# Table of colors from page 718
if (color_choice=="complex") { 
  colorcode = c("lightgreen",   #EMPTY
                "black",        #SPORE
                "darkgrey",     #YOUNG
                "lightgrey",    #MATURING
                "white",        #MUSHROOM(S)
                "lightgrey",    #OLDER
                "tan",          #DECAYING
                "brown",        #DEAD1
                "darkgreen",    #DEAD2
                "yellow")       #INERT
}
if (color_choice=="simple") {
  colorcode = c("lightgreen",   #EMPTY
                "black",        #SPORE
                "lightgreen",   #YOUNG
                "lightgreen",   #MATURING
                "white",        #MUSHROOM(S)
                "lightgreen",   #OLDER
                "lightgreen",   #DECAYING
                "lightgreen",   #DEAD1
                "lightgreen",   #DEAD2
                "yellow")       #INERT
}

# Initialize to square of spore cells in the center
for (row in 1:n) {
    for (col in 1:n) {
      # initializing spores
      
      # selection of percentage of rows
      # if (row < n * 0.20 && row > n * 0.15 && col < n * 0.50 && col > n * 0.45) {
      #   ts[row,col,1] = 1  
      #   colormat_s[row,col,1] = colorcode[2]
      # } 
      
      # selecting a point for a spore
      if (row==round(0.15*n) && col==round(0.75*n)) {
        ts[row,col,1] = 1
      }
      
      # initializing obstructions with inert space
      if (row < n * 0.65 && row > n * 0.45 && col < n * 0.10 && col > n * 0.00) { 
        ts[row,col,1] = 9
        colormat_s[row,col,1] = colorcode[10]}
      if (row < n * 0.70 && row > n * 0.50 && col < n * 0.40 && col > n * 0.20) { 
        ts[row,col,1] = 9
        colormat_s[row,col,1] = colorcode[10]}
      if (row < n * 0.75 && row > n * 0.55 && col < n * 0.70 && col > n * 0.50) { 
        ts[row,col,1] = 9
        colormat_s[row,col,1] = colorcode[10]}
      if (row < n * 0.80 && row > n * 0.60 && col < n * 1.00 && col > n * 0.80) { 
        ts[row,col,1] = 9
        colormat_s[row,col,1] = colorcode[10]}
      else {
            colormat_s[row,col,1] = colorcode[1]
        }
    }
}

# Random initialization
# for (row in 1:n) {
#     for (col in 1:n) {
#         if (runif(1) < probSpore) {
#             mat[row,col] = 1
#             colormat[row,col] = colorcode[2]
#         } else {
#             mat[row,col] = 0
#             colormat[row,col] = colorcode[1]
#         }
#     }
# }


# Messy function for determining if any of a specified cell's neighbors are "young" (state 2)
# Several problem variants possible just by modifying which cells are checked here
# Thought it could be a single boolean expression (commented at the end), but didn't seem to be boolean short-circuiting properly
hasYoungNeighbor <- function(row, col,j) {
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
    
    #return ((row > 1 && ((col > 1 && mat[row-1,col-1] == 2) || (mat[row-1,col] == 2) || (col < n && mat[row-1,col+1] == 2)))
    #|| ((col > 1 && mat[row,col-1] == 2) || (col < n && mat[row,col+1] == 2))
    #|| (row < n && ((col > 1 && mat[row+1,col+1] == 2) || (mat[row+1,col] == 2) || (col < n && mat[row+1,col+1] == 2))))
}	

# ---------------------------------------------------------------------------------------------------------------


# Run 1 iterations
for (t in 2:iter) {
    # Delay some amount of time each iteration (probably not the best way to do this)
    #Sys.sleep(0.5)
    
    # Iterate over every cell
    for (row in 1:n) {
        for (col in 1:n) {
            # Follow different rules based on old state, following state diagram on page 717
            oldstate = ts[row,col,t-1]
            ts[row,col,t]= case_when(oldstate==0 ~ # cell is previously EMPTY
                                         if (runif(1) < probSpread && hasYoungNeighbor(row,col,t-1)) {
                                             2 # becomes young
                                         } else {
                                             0 # stays empty
                                         },
                                     oldstate==1 ~ # cell is previously SPORE
                                         if (runif(1) < probSporeToHyphae) {
                                             2 # becomes young
                                         } else {
                                             1 # stays spore
                                         },
                                     oldstate==2 ~ 3, # maturing cell is previously YOUNG
                                     oldstate==3 ~ if (runif(1) < probMushroom) {
                                         4 # mushrooms
                                     } else {
                                         5 # older
                                     },
                                     oldstate==4 ~ 6, # decaying1
                                     oldstate==5 ~ 6, # decaying2
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
library(animation)

saveGIF({
  for(i in 1:iter){
    color2D.matplot(matrix(1:nsq,nrow=n),cellcolors=colormat_s[,,i],yrev=FALSE,
                    border="lightgray",xlab="",ylab="",
                    main="Shroomage over time",
                    sub=paste("Time = Approx. Year #",1+(i%/%4)),
                    axes=FALSE)
  }
}, interval = speed, ani.width = 600, ani.height = 600)  #relate this to real time scale with multiplier for simulation


cat('measuring distances')
nPoints = 0
tVals = c()
rVals = c()
for (t in 1:iter) {
  for (x in 1:n) {
    for (y in 1:n) {
      if (ts[x,y,t] == 4) {
        nPoints = nPoints + 1
        tVals[nPoints] = t
        rVals[nPoints] = sqrt((x - (n/2))^2 + (y - (n/2))^2)
      }
    }
  }
}

tVals_real = tVals/4
rVals_real = rVals*2
plot(tVals_real, rVals_real,xlab="Time (Years)",ylab="Distance from Center of Plane (in)",
     main=sprintf("Mushroom distances - probSpread=%.2f", probSpread), xlim=c(0,max(tVals_real)), ylim=c(0,max(rVals_real)))


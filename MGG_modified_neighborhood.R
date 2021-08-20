
library(plotrix)
library(tidyverse)
library(abind)

## Tunable probability parameters 
probSporeToHyphae = 0.3
probMushroom = 0.7
probSpread = 0.9
probSpore = 0.2
# new constant describing proportionally lower spread to diagonal cells than adjacent ones
diagonalSpreadFactor = 0.3

## Graphics Parameters
speed = 0.5  # GIF output transition interval
n = 70       # Side length for square grid
iter = 50    # Duration of simulation
nsq = n*n

## Initialization parameters
initPointRow = floor(n/2) 
initPointCol = floor(n/2)
initSquareSize = 3

# Initialize matrices for cell states and colors
ts <- array(numeric(), c(n, n, iter))
colormat_s <- array(character(), c(n, n, iter))
colormat_s[,,]=0
ts[,,]=0


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
  } else if (hasYoungMooreNeighbor(row, col, j) && randomValue < probSpread * diagonalSpreadFactor) {
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

# Create GIF
library(animation)

saveGIF({
  for(i in 1:iter){
    color2D.matplot(matrix(1:nsq,nrow=n),cellcolors=colormat_s[,,i],yrev=FALSE,
                    border="lightgray",xlab="",ylab="",main="shroomage over time",axes=FALSE)
  }
}, interval = speed, ani.width = 550, ani.height = 350)  #relate this to real time scale with multiplier for simulation

# Plot mushroom distances over time
# cat('Measuring distances...')
# nPoints = 0
# tVals = c()
# rVals = c()
# for (t in 1:iter) {
#   for (x in 1:n) {
#     for (y in 1:n) {
#       if (ts[x,y,t] == 4) {
#         nPoints = nPoints + 1
#         tVals[nPoints] = t
#         rVals[nPoints] = sqrt((x - initPointRow)^2 + (y - initPointCol)^2)
#       }
#     }
#   }
# }
# 
# plot(tVals, rVals,xlab="time (iters)",ylab="distance (cells)",
#      main=sprintf("Mushroom distances - probSpread=%.2f", probSpread), xlim=c(0,iter), ylim=c(0,n / sqrt(2)))


# Plot mushroom distances vs angle, for a small interval of time
cat('Measuring distances vs angle...')
nPoints = 0
angVals = c()
rVals = c()
timeIntervalStart = floor(n / 2.7) + 2 # try to start when the ring is at a good radius
for (t in timeIntervalStart:(timeIntervalStart+2)) {
  for (x in 1:n) {
    for (y in 1:n) {
      if (ts[x,y,t] == 4) {
        # Mushroom here, so record its polar coords for plotting
        nPoints = nPoints + 1
        angVals[nPoints] = atan2(y - initPointCol, x - initPointRow)
        rVals[nPoints] = sqrt((x - initPointRow)^2 + (y - initPointCol)^2)
      }
    }
  }
}

plot(angVals, rVals,xlab="angle (rad)",ylab="distance (cells)",
     main=sprintf("Mushroom distance vs angle - probSpread=%.2f, alpha=%.2f", probSpread, diagonalSpreadFactor),
     xlim=c(-pi,pi), ylim=c(0,n / sqrt(2)))

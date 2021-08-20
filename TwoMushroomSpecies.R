# mushroom_test1.R - quick initial implementation of the fairy ring cellular automaton

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


n = 100  # Side length for square grid
nsq = n*n
iter = 100

ts <- array(numeric(), c(n, n, iter))
colormat_s <- array(character(), c(n, n, iter))
colormat_s[,,]=0
ts[,,]=0


# Tunable probability parameters
probSporeToHyphae_S1 = 0.3
probMushroom_S1 = 0.7
probSpread_S1 = 0.9
probSpore_S1 = 0.05

probSporeToHyphae_S2 = 0.3
probMushroom_S2 = 0.7
probSpread_S2 = 0.7
probSpore_S2 = 0.1

diagonalSpreadFactor = 0.3

# Table of colors from page 718
# colorcode = c("lightgreen",
#               "black",
#               "darkgrey",
#               "lightgrey",
#               "white",
#               "lightgrey",
#               "tan",
#               "brown",
#               "darkgreen",
#               "yellow")

# Table of colors from page 718
colorcode = c("lightgreen", 
              "lightgreen", 
              "lightgreen", 
              "lightgreen", 
              "white", 
              "lightgreen", 
              "lightgreen", 
              "lightgreen", 
              "lightgreen", 
              "lightgreen", 
              "lightgreen", 
              "lightgreen", 
              "lightgreen", 
              "red"
)

# Initialize to square of spore cells in the center
# for (row in 1:n) {
#   for (col in 1:n) {
#     if (row < n * 0.85 && row > n * 0.8 && col < n * 0.85 && col > n * 0.8) {
#       ts[row,col,1] = 1  
#       colormat_s[row,col,1] = colorcode[2]
#     } 
#     else if (row < n * 0.25 && row > n * 0.2 && col < n * 0.25 && col > n * 0.2) {
#       ts[row,col,1] = 10  
#       colormat_s[row,col,1] = colorcode[11]
#     }  
#     else {
#       colormat_s[row,col,1] = colorcode[1]
#     }
#   }
# }

#Random initialization
for (row in 1:n) {
    for (col in 1:n) {
        if (runif(1) < probSpore_S1) {
            ts[row,col,1] = 1
            colormat_s[row,col,1] = colorcode[2]}
        else if (runif(1) < probSpore_S2) {
              ts[row,col,1] = 10
              colormat_s[row,col,1] = colorcode[11]
        } else {
            ts[row,col,1] = 0
            colormat_s[row,col,1] = colorcode[1]
        }
    }
}


hasYoungMooreNeighbor_S1 <- function(row, col,j) {
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

hasYoungMooreNeighbor_S2 <- function(row, col,j) {
  if (row > 1) {
    if (col > 1) {
      if (ts[row-1,col-1,j] == 11) {
        return(TRUE)
      }
    }
    if (col < n) {
      if (ts[row-1,col+1,j] == 11) {
        return(TRUE)
      }
    }
    if (ts[row-1,col,j] == 11) {
      return(TRUE)
    }
  }
  if (row < n) {
    if (col > 1) {
      if (ts[row+1,col-1,j] == 11) {
        return(TRUE)
      }
    }
    if (col < n) {
      if (ts[row+1,col+1,j] == 11) {
        return(TRUE)
      }
    }
    if (ts[row+1,col,j] == 11) {
      return(TRUE)
    }
  }
  if (col > 1) {
    if (ts[row,col-1,j] == 11) {
      return(TRUE)
    }
  }
  if (col < n) {
    if (ts[row,col+1,j] == 11) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# Check if the cell (row, col) has a Young cell in its Von Neumann neighborhood at iteration j
hasYoungVonNeumannNeighbor_S1 <- function(row, col, j) {
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

hasYoungVonNeumannNeighbor_S2 <- function(row, col, j) {
  if (col > 1 && ts[row,col-1,j] == 11) {
    return(TRUE)
  }
  if (row > 1 && ts[row-1,col,j] == 11) {
    return(TRUE)
  }
  if (col < n && ts[row,col+1,j] == 11) {
    return(TRUE)
  }
  if (row < n && ts[row+1,col,j] == 11) {
    return(TRUE)
  }
  return(FALSE)
}

# Determined whether the fungus should spread to cell (row, col) on iteration j
spreadToCell_S1 <- function(row, col, j) {
  randomValue = runif(1)
  if (hasYoungVonNeumannNeighbor_S1(row, col, j)) {
    # There's a young cell directly adjacent, so use higher probability
    if (randomValue < probSpread_S1) {
      return(TRUE)
    }
  } else if (hasYoungMooreNeighbor_S1(row, col, j) && randomValue < probSpread_S1 * diagonalSpreadFactor) {
    return(TRUE)
  }
  return(FALSE)
}	

spreadToCell_S2 <- function(row, col, j) {
  randomValue = runif(1)
  if (hasYoungVonNeumannNeighbor_S2(row, col, j)) {
    # There's a young cell directly adjacent, so use higher probability
    if (randomValue < probSpread_S2) {
      return(TRUE)
    }
  } else if (hasYoungMooreNeighbor_S2(row, col, j) && randomValue < probSpread_S2 * diagonalSpreadFactor) {
    return(TRUE)
  }
  return(FALSE)
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
      
      ts[row,col,t]= case_when(oldstate==0 ~ 
                                 if (spreadToCell_S1(row,col,t-1)) {
                                   2 } # young species 2
                                 else if (spreadToCell_S2(row,col,t-1)) {
                                   11} # young species 1
                                  else {
                                   0 # empty
                                 },
                               oldstate==1 ~ # cell is previously SPORE S1
                                 if (runif(1) < probSporeToHyphae_S1) {
                                   2} # young
                                  else {
                                   1 # spore
                                 },
                               oldstate==10 ~ # cell is previously SPORE S2
                                 if (runif(1) < probSporeToHyphae_S2) {
                                   11} # young
                                  else {
                                   10 # spore
                               },
                               oldstate==2 ~ 3, # maturing# cell is previously YOUNG S1
                               oldstate==11 ~ 12, # maturing# cell is previously YOUNG S2
                               oldstate==3 ~ if (runif(1) < probMushroom_S1) {
                                 4 # mushrooms S1
                               } else {
                                 5 # older S1/S2
                               },
                               oldstate==12 ~ if (runif(1) < probMushroom_S2) {
                                 13 # mushrooms S2
                               } else {
                                 5 # older S1/S2
                               },
                               oldstate==4 ~ 6, # decaying S1
                               oldstate==13 ~ 6, # decaying S2
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

library(animation)

saveGIF({
  for(i in 1:iter){
    color2D.matplot(matrix(1:nsq,nrow=n),cellcolors=colormat_s[,,i],yrev=FALSE,
                    border="lightgray",xlab="",ylab="",main="Mushroom Growth",axes=FALSE)
  }
}, interval = 0.5, ani.width = 550, ani.height = 550)



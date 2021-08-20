# CSE616 Simulation of Physical Systems
# Mushroom Cellular Automota Code -- Probabilistic Death
# Indrima Upadhyay
# Jayson Rook
# Forked by Jacob Veta
#
#
library(plotrix)
library(tidyverse)
library(abind)

## Tunable probability parameters #find research and convert to 2d matrix for comparison  
probSporeToHyphae = 0.3
probMushroom = 0.7
probSpread = 0.9
probSpore = 0.2
probDeath = 0.5
iter = 65


## Graphics Parameters

speed = 0.40  #GIF output transition interval
n = 65  # Side length for square grid
nsq = n*n

ts <- array(numeric(), c(n, n, iter))
colormat_s <- array(character(), c(n, n, iter))
colormat_s[,,]=0
ts[,,]=0



# Table of color
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

# Initialize to square of spore cells in the center
for (row in 1:n) {
    for (col in 1:n) {
        if (row==round(0.5*n) && col==round(0.5*n)) {
            ts[row,col,1] = 1  
            colormat_s[row,col,1] = colorcode[2]
        } else {
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
            ts[row,col,t]= case_when(oldstate==0 ~ 
                                         if (runif(1) < probSpread && hasYoungNeighbor(row,col,t-1)) {
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
                                     oldstate==6 ~ if (runif(1) < probDeath) {
                                       7 # dead1
                                     } else {
                                       6 # Decaying
                                     },
                                     oldstate==7 ~ if (runif(1) < probDeath ){
                                       8 # dead2
                                     } else {
                                       7 # dead1
                                     },
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
                    border="lightgray",xlab="",ylab="",main="mushroom growth over time",axes=FALSE)
  }
}, interval = speed, ani.width = 550, ani.height = 350)  #relate this to real time scale with multiplier for simulation


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



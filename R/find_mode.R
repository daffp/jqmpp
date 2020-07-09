#' Find the nearest maximum in a matrix
#' 
#' Find the nearest maximum in a matrix to a given point by hillclimbing.
#'
#' @param A A matrix.
#' @param px,py A point within the dimensions of the matrix.
#' @return Returns the x and y matrix positions of the nearest maximum that was found.
#' @details A hillclimbing search from the input point across the array to find the nearest local maximum. 
#' The search is through each square surrounding the next input point in the search (accounting for the array boundaries).
#' Ties are broken by taking the nearest euclidean distance (i.e. point which are vertical or horizontal from the input rather than diagonal)
#' and then random sampling. To make this reproducible the seed should be set.
#' @examples
#' set.seed(444)
#' msize <- 20
#' A <- matrix(sample(seq_len(msize), msize^2, replace=TRUE), ncol=msize)
#' find_nearest_mode2D(A, 10, 1)
#' find_nearest_mode2D(A, 2, 4)
#' @importFrom stats dist
#' @export
find_nearest_mode2D = function(A, px, py){
  
  # Dimensions of A -- needed for boundary check
  nc = ncol(A)
  nr = nrow(A)
  
  # value at current point 
  current_max = A[px, py]
  
  # NEIGHBOURS: 
  pos = create_indices_for_submat(px, py, nr, nc)
  submat = A[pos$x, pos$y]
  max_submat = max(submat)   
  
  # REPEAT
  repeat{
    
    # CHECK CONDITION
    if(current_max == max_submat) break # RETURN CURRENT POSITION
    
    # RESET MAX
    current_max = max_submat
    
    # need to get position of max in main matrix & break ties
    max_subpos = which(submat == max_submat, arr.ind=TRUE)
    new_pos = cbind(pos$x[max_subpos[,1]], pos$y[max_subpos[,2]])
    dist_to_max = unname(as.matrix(dist(rbind(c(px, py), new_pos))))[,1][-1]
    closest_max = which(dist_to_max == min(dist_to_max))
    new_pos = new_pos[sample(closest_max, 1),]
    px = new_pos[1]; py = new_pos[2]
    
    # NEIGHBOURS
    pos = create_indices_for_submat(px, py, nr, nc)
    submat = A[pos$x, pos$y]
    max_submat = max(submat)   
  }
  
  return(list(x=px, y=py))
}


#' Find the nearest maximum in an array
#' 
#' Find the nearest maximum in an array to a given point by hillclimbing.
#'
#' @param A A 3D array.
#' @param px,py,pz A point within the dimensions of the array.
#' @return Returns the x, y and z array positions of the nearest maximum that was found.
#' @details A hillclimbing search from the input point across the array to find the nearest local maximum. 
#' The search is through each cube surrounding the next input point in the search (accounting for the array boundaries).
#' Ties are broken by taking the nearest euclidean distance (i.e. point which are vertical or horizontal from the input rather than diagonal)
#' and then random sampling. To make this reproducible the seed should be set.
#' @examples
#' set.seed(444)
#' msize <- 20
#' A <- array(sample(seq_len(msize), msize^3, replace=TRUE), dim = rep(msize, 3))
#' find_nearest_mode3D(A, 3,3,3)
#' find_nearest_mode3D(A, 4,1,5)
#' @export
find_nearest_mode3D = function(A, px, py, pz){
  
  # Dimensions of A -- needed for boundary check
  d = dim(A)
  nr = d[1]; nc = d[2]; nh = d[3]
  
  # value at current point 
  current_max = A[px, py, pz]
  
  # NEIGHBOURS: 
  pos = create_indices_for_subarray(px, py, pz, nr, nc, nh)
  subarray = A[pos$x, pos$y, pos$z]
  max_subarray = max(subarray)   
  
  # REPEAT
  repeat{
    
    # CHECK CONDITION
    if(current_max == max_subarray) break # RETURN CURRENT POSITION
    
    # RESET MAX
    current_max = max_subarray
    
    # need to get position of max in main matrix & break ties
    max_subpos = which(subarray == max_subarray, arr.ind=TRUE)
    new_pos = cbind(pos$x[max_subpos[,1]], pos$y[max_subpos[,2]], pos$z[max_subpos[,3]])
    dist_to_max = unname(as.matrix(dist(rbind(c(px, py, pz), new_pos))))[,1][-1]
    closest_max = which(dist_to_max == min(dist_to_max))
    new_pos = new_pos[sample(closest_max, 1),]
    px = new_pos[1]; py = new_pos[2]; pz=new_pos[3]
    
    # NEIGHBOURS
    pos = create_indices_for_subarray(px, py, pz, nr, nc, nh)
    subarray = A[pos$x, pos$y, pos$z]
    max_subarray = max(subarray)      
  }
  
  return(list(x=px, y=py, z=pz))
}



# Create indices range used to search neighbours in a matrix
# accounts for boundary effects 
# x,y: point position
# NRA: number of rows of array A (first dim)
# NCA: number of columns of array A (2nd)
create_indices_for_submat <- function(x, y, NRA, NCA){
  minx = pmax(1, x-1); maxx = pmin(NRA, x+1)
  miny = pmax(1, y-1); maxy = pmin(NCA, y+1)
  xrange = minx:maxx
  yrange = miny:maxy
  return(list(x=xrange, y=yrange))
}


# Create indices range used to search neighbours in a 3D array
# accounts for boundary effects 
# x,y,z: point position
# NRA: number of rows of array A (first dim)
# NCA: number of columns of array A (2nd)
# NHA: height of array A(3rd)  
create_indices_for_subarray <- function(x, y, z, NRA, NCA, NHA){
  minx = max(1, x-1); maxx = min(NRA, x+1)
  miny = max(1, y-1); maxy = min(NCA, y+1)
  minz = max(1, z-1); maxz = min(NHA, z+1)
  xrange = minx:maxx
  yrange = miny:maxy
  zrange = minz:maxz
  return(list(x=xrange, y=yrange, z=zrange))
}


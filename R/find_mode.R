#' Find the nearest maximum in an array
#' 
#' Find the nearest maximum in an array to a given point by hillclimbing.
#'
#' @param A An array.
#' @param ... The coordinates of a point within the dimensions of the array.
#' @return Returns a list of the array positions of the nearest maximum that was found.
#' @details A hillclimbing search from the input point across the array to find the nearest local maximum. 
#' The search is through each unit n-cube surrounding the next input point in the search (accounting for the array boundaries).
#' Ties are broken by taking the nearest euclidean distance (i.e. point which are vertical or horizontal from the input rather than diagonal)
#' and then random sampling. To make this reproducible the seed should be set.
#' @examples
#' set.seed(444)
#' msize <- 20
#' A <- array(sample(seq_len(msize), msize^3, replace=TRUE), dim = rep(msize, 3))
#' find_nearest_mode(A, 3,3,3)
#' find_nearest_mode(A, 4,1,5)
#' @importFrom stats dist
#' @export
find_nearest_mode = function(A, ...){
  
  # The coordinates
  mouse = list(...)
  
  # Dimensions of A -- needed for boundary check
  dm = dim(A)
  
  # value at current point 
  current_max = do.call(`[`, c(list(A), mouse))
  
  # NEIGHBOURS: 
  pos = create_indices_for_subarray(dm, mouse)
  subarray = do.call(`[`, c(list(A), pos))
  max_subarray = max(subarray)   
  new_pos = mouse # return value if already at max
  
  # REPEAT
  repeat{
    
    # CHECK CONDITION
    if(current_max == max_subarray) break # RETURN CURRENT POSITION
    
    # RESET MAX
    current_max = max_subarray
    
    # need to get position of max in main matrix & break ties
    max_subpos = which(subarray == max_subarray, arr.ind=TRUE)
    new_pos = matrix(mapply(function(X,Y) X[Y], X=pos, Y=asplit(max_subpos, MARGIN=2)), ncol=length(mouse))
    dist_to_max = unname(as.matrix(dist(rbind(unlist(mouse), new_pos))))[,1][-1]
    closest_max = which(dist_to_max == min(dist_to_max))
    new_pos = new_pos[sample(closest_max, 1),]
    
    # NEIGHBOURS
    pos = create_indices_for_subarray(dm, new_pos)
    subarray = do.call(`[`, c(list(A), pos))
    max_subarray = max(subarray)   
  }
  
  return(as.list(new_pos))
}


# Create indices range used to search neighbours in an array
# accounts for boundary effects 
# dims: the dimensions of the array
# pts:  An ordered list or vector of positions to subset an array i.e. pts=list(x=1, y=2, z=3). It isn't required to be named.
create_indices_for_subarray <- function(dims, pts){
  pos = Map(function(SIDE, point) max(1, point-1): min(SIDE, point+1), SIDE=dims, point=pts)
  return(pos)
}



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


#' Find the nearest maximum in a scaled array
#' 
#' Find the nearest maximum in an array to a given point by hillclimbing accounting for different distances in each dimension.
#'
#' @param A An array.
#' @param scale A vector giving the value for one unit of distance. The length of \code{state} should be the same as the number of dimensions in \code{A}.
#' @param ... The coordinates of a point within the dimensions of the array.
#' @return Returns a list of the array positions of the nearest maximum that was found.
#' @details A hillclimbing search from the input point across the array to find the nearest local maximum. 
#' The search is through each unit n-cube surrounding the next input point in the search (accounting for the array boundaries).
#' Differences in the scaled distance along each dimension are accounted for by finding the next maximum which is scaled by the distance e.g. (current_max - neighbour_max)/distance,
#' and select the value with the largest gradient. Ties are broken by random sampling. To make this reproducible the seed should be set.
#' @examples
#' set.seed(444)
#' msize <- 20
#' A <- array(sample(seq_len(msize), msize^3, replace=TRUE), dim = rep(msize, 3))
#' find_nearest_scaled_mode(A, c(1,1,1), 11, 5, 5)
#' @importFrom stats dist
#' @export
find_nearest_scaled_mode = function(A, scale, ...){ 

  # The coordinates
  mouse = list(...)
  
  # Dimensions of A -- needed for boundary check
  dm = dim(A)
  
  # value at current point 
  current_max = do.call(`[`, c(list(A), mouse))

  # NEIGHBOURS: 
  pos = create_indices_for_subarray(dm, mouse)
  subarray = do.call(`[`, c(list(A), pos))
  
  # max scaled subarray: just need to grab one max if ties
  steepest_inc = get_steepest_increase(dm, mouse, scale, subarray)
  max_subarray = do.call(`[`, c(list(subarray), steepest_inc$max_ind[1,]))  
  new_pos = mouse # the return value if already at max  
  
   # REPEAT
  repeat{

    # CHECK CONDITION
    if(current_max == max_subarray) break # RETURN CURRENT POSITION

    # RESET MAX
    current_max = max_subarray 
    
    # need to get position of max in main matrix & break ties
    # just break tie randomly if equal gradient (or sh/could select nearest)
    max_subpos = which(subarray == max_subarray, arr.ind=TRUE)
    new_pos = matrix(mapply(function(X,Y) X[Y], X=pos, Y=asplit(max_subpos, MARGIN=2)), ncol=length(mouse))
    new_pos = new_pos[sample(nrow(new_pos), 1),]
    
    # NEIGHBOURS
    pos = create_indices_for_subarray(dm, new_pos)
    subarray = do.call(`[`, c(list(A), pos))
    steepest_inc = get_steepest_increase(dm, new_pos, scale, subarray)
    max_subarray = do.call(`[`, c(list(subarray), steepest_inc$max_ind[1,]))  
  }
  
  return(as.list(new_pos))
}  


# Helper functions

# Create indices range used to search neighbours in an array
# accounts for boundary effects 
# dims: the dimensions of the array
# pts:  An ordered list or vector of positions to subset an array i.e. pts=list(x=1, y=2, z=3). It isn't required to be named.
create_indices_for_subarray <- function(dims, pts){
  pos = Map(function(SIDE, point) max(1, point-1): min(SIDE, point+1), SIDE=dims, point=pts)
  return(pos)
}


# We now need to create delta(intensity)/distance and select the maximum
# dims: dim of A
# mouse:  current position
# rescale: the interval size in each array dimension
# subarray: the array of values surrounding the current position
get_steepest_increase <- function(dims, mouse, rescale, subarray){
    
          # Distances of surrounding points to current position
          indices_of_sub_array = t(which(subarray >= -Inf, arr.ind=TRUE))
          indices_for_current_position = indices_for_current_position_in_subarray(dims, mouse)
          
          # distance in matrix units
          dist_mat_unit = abs(indices_of_sub_array - indices_for_current_position)
        
          # multiply this by the scaling for each axis
          dist_real_unit = rescale* dist_mat_unit
          
          # calculate distance
          distances = sqrt(colSums(dist_real_unit^2))
          attributes(distances) = attributes(subarray)
          
          # Find steepest increase
          current_max = do.call(`[`, c(list(subarray), indices_for_current_position))
          scaled_distance = (subarray - current_max) / distances  
          scaled_distance[is.nan(scaled_distance)] = 0
          
          # Get index where the scaled maximum is
          max_subarray = max(scaled_distance) 
          max_subpos = which(scaled_distance == max_subarray, arr.ind=TRUE)     
          
          # output
          lst = list(max_ind = max_subpos, scales_dist = scaled_distance)
          
          return(lst)
        }

  
# Find the index of the current position in the subsetted array
# dims: dimensions of A
# mouse: the current position in intensity matrix dimension
indices_for_current_position_in_subarray <- function(dims, mouse){
  subarrID = create_indices_for_subarray(dims = dims, pts = mouse)
  mapply(function(X,Y) match(X,Y), X=mouse, Y=subarrID)
}

# This module consists of all utility functions.

get_circular_index = function(pos, step, size){
  new = pos + step
  if (new<1){
    new = new + size
  } else if (new>size){
    new = new - size
  }
  
  return(new)
}

get_beta_groups = function(beta_hat, tol=1e-8){
  group = rep(NA, N)
  group[1] = 1     # initialize first beta into first group
  
  for (i in 2:N){
    is_grouped = FALSE      # flag whether beta_i has been assigned a group
    collected_groups = unique(group)     # groups gathered so far
    collected_groups = collected_groups[!is.na(collected_groups)]   # remove NA
    
    for (l in collected_groups){
      is_group_l = TRUE     # flag whether beta_i belongs to group_l
      members = beta_hat[,which(group==l), drop=FALSE]     # all members of group l
      
      # if beta_i is not equal to ANY member, then not group_l
      for (j in 1:ncol(members)){       
        if (max(abs(beta_hat[,i] - members[,j])) > tol){
          is_group_l = FALSE
          break
        }
      }
      
      # if beta_i is equal to ALL members, then assign to group_l
      if (is_group_l){
        group[i] = l
        is_grouped = TRUE
        break
      }
    }
    
    # if beta_i does not belong to any collected groups, then assign to a new group
    if (!is_grouped) group[i] = max(collected_groups) + 1
  }
  
  return(group)
}

assert_dimensions = function(A, B=NULL, dim_=NULL){
  dim_A = if (is.array(A) || inherits(A, "CVXR::Expression")) dim(A) else length(A)
  
  if (!is.null(B) && is.null(dim_)){
    dim_B = if (is.array(B)) dim(B) else length(B)
    stopifnot(identical(dim_A, dim_B))
    
  } else if (is.null(B) && (!is.null(dim_))){
    dim_ = as.integer(dim_)
    stopifnot(identical(dim_A, dim_))
    
  } else{
    stop("One and only one of B or dim_ should be supplied. Or some other error.")
  }
}

# # Some test cases
# # Good
# assert_dimensions(matrix(NA,2,3), matrix(NA,2,3))
# assert_dimensions(matrix(NA,2,3), dim_=c(2,3))
# assert_dimensions(array(NA, c(1,2,3)), array(NA, c(1,2,3)))
# assert_dimensions(array(NA, c(1,2,3)), dim_=c(1,2,3))
# assert_dimensions(c(NA,NA), c(3,5))
# assert_dimensions(c(NA,NA), dim_=2)
# 
# # Bad
# assert_dimensions(matrix(NA,2,3), matrix(NA,3,3))
# assert_dimensions(matrix(NA,2,3), dim_=c(3,3))
# assert_dimensions(matrix(NA,2,3), matrix(NA,3,3), c(1,2))
# assert_dimensions(matrix(NA,2,3))
# assert_dimensions(array(NA, c(1,2,3)), array(NA, c(1,2,4)))
# assert_dimensions(array(NA, c(1,2,3)), dim_=c(1,2,4))
# assert_dimensions(c(1,2), c(3,4,5))
# assert_dimensions(c(NA,NA), dim_=3)









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

list2str = function(x, digits=NULL){
  out = ""
  for (item in names(x)){
    val = ifelse(is.null(digits), x[[item]], round(x[[item]], digits))
    out = paste0(out, item, "=", val, ", ")
  }
  out = substring(out, 1, nchar(out)-2) # remove trailing comma
  return(out)
}

clog = function(msg, file=LOG_FILE){
  logline = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%OS6 "), msg, "\n")
  cat(logline)   # print to console
  cat(logline, file=file, append=TRUE)    # print to log
}

make_csv_beta = function(N, p){
  header = ""
  for (i in 1:N){
    for (k in 1:p){
      header = paste0(header, "beta_", i, "_", k, ",")
    }    
  }
  header = paste0(header, "\n")
  cat(header, file=CSV_BETA, append=TRUE)
}

timer = function(f){
  # Can only log scalar and named parameters
  fn_name = deparse(substitute(f))
  force(f)
  function(...){
    start = Sys.time()
    fn_params = list(...)
    clog(paste0(fn_name, ": status=start, ", list2str(fn_params)))
    output = f(...)
    
    elapsed = difftime(Sys.time(), start, units="secs")
    clog(paste0(fn_name, ": status=end, elapsed=", round(elapsed, 2), ", ", list2str(fn_params)))
    return(output)
  }
}

read_args = function(args){
  args = sub("^--", "", args)
  params = strsplit(args, "=", fixed=TRUE)
  
  params = setNames(
    lapply(params, `[`, 2),
    sapply(params, `[`, 1)
  )
  
  out = list(
    n_rep=as.integer(params$n_rep),
    N=as.integer(params$N),
    T_=as.integer(params$T_),
    r=as.integer(params$r),
    p=as.integer(params$p),
    phi=as.numeric(params$phi),
    K0=as.integer(params$K0),
    pi_=as.numeric(params$pi_),
    L=as.integer(params$L),
    kappa=as.integer(params$kappa),
    is_test=as.logical(params$is_test)
  )
  
  return(out)
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


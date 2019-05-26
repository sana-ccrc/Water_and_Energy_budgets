
library(raster)
library(lubridate)
library(numDeriv)
library(gmp)
library(plyr)

## Calculate Rn and ds given the value of lh, h, g, p and q
r_obs_function <- function(x){
  
  
  ### Matrix representing the right side of the water and energy budget equations where,
  ### Rn = 1.lh + 1.h + 1.g + 0.p + 0.q
  ### ds - 1.lh + 0.h + 0.g + 1.p + 0.q
  A <- rbind(c(1, 1, 1, 0, 0),c(-1, 0, 0, 1, -1)) 
  
  return(c(A%*% cbind(x)))
}   

### Enforce the closure of the water and energy budget simulatneously at a given grid cell
### Calculate the optimized fluxes and their uncertainites as well as the goodness of the fit of the original fluxes 
### Each flux is adjusted based on its relative uncertainty

Solve_Budgets_Grid <- function(a){
  

  lh <- a[1] ## latent heat flux
  h <- a[2] ## sensible heat flux
  g <- a[3] ## ground heat flux
  p <- a[4] ## precipitation
  q <- a[5] ## runoff
  rn <- a[6] ## net radiation flux
  ds <- a[7] ## change in water storage
  lh_u <- a[8] ## uncertainty of latent heat flux
  h_u <- a[9]  ## uncertainty of sensible heat flux
  g_u <- a[10] ## uncertainty of ground heat flux 
  p_u <- a[11] ## uncertainty of precipitation
  q_u <- a[12] ## uncertainty of runoff
  rn_u <- a[13] ## uncertainty of net radiation
  ds_u <- a[14] ## uncertainty of change in water storage
  
  
  ### vector of fluxes that need to be optimized
  f_obs_v <- c(lh, h, g, p, q)
  
  ### vector of storage terms
  r_obs_v <- c(rn, ds)
  
  ### vectors uncertainties
  f_u_v <- c(lh_u, h_u, g_u, p_u, q_u)
  r_u_v <- c(rn_u, ds_u)
  
  ### this transposes vectors fluxes and storage terms
  f_obs <- cbind(f_obs_v)
  r_obs <- cbind(r_obs_v)
  
  
  
  ### Matrix representing the right side of the water and energy budget equations where,
  ### Rn = 1.lh + 1.h + 1.g + 0.p + 0.q
  ### ds - 1.lh + 0.h + 0.g + 1.p + 0.q
  A <- rbind(c(1, 1, 1, 0, 0),c(-1, 0, 0, 1, -1)) ### matrix 
  
  
  ### calculate the gradient of r_obs_function at the location of vector of fluxes
  k <- jacobian(r_obs_function, f_obs_v, method="complex" )
  
  
  ## covariance of robs
  s_robs <- matrix(0,nrow=2, ncol=2)
  # s_robs[1,1] <- (rn_u)^2
  # s_robs[2,2] <- (ds_u)^2
  
  s_robs[1,1] <- (r_u_v[1])^2
  s_robs[2,2] <- (r_u_v[2])^2
  
  ## covariance of fobs 
  s_fobs <-  matrix(0,nrow=5, ncol=5)
  ## diagonals equal to error variance of the fluxes, i.e. uncertainty^2, as uncertainties here are standard deviation
  s_fobs[1,1] <- (f_u_v[1])^2
  s_fobs[2,2] <- (f_u_v[2])^2
  s_fobs[3,3] <- (f_u_v[3])^2
  s_fobs[4,4] <- (f_u_v[4])^2
  s_fobs[5,5] <- (f_u_v[5])^2
  
  
  ### This is necessary when one of the fluxes is NA. Otherwise an error will be generated
  ### In case one flux is NA, return NA for all the optimized terms
  out <- tryCatch(solve(s_robs) %*% s_robs, error = function(e) e)
  result <- any(class(out) == "error")
  if(result == TRUE){ ## there are errors
    return(rep(NA,13))
  }
  out <- tryCatch(solve(s_fobs) %*% s_fobs, error = function(e) e)
  result <- any(class(out) == "error")
  if(result == TRUE){ ## there are errors
    return(rep(NA,13))
  }
  
  
  ### Error covariance for the component fluxes after optimization
  #S_f <- (k(transpose) x s_robs(inverse) x k + S_fobs(inverse))(inverse)
  s_f <- t(k) %*% solve(s_robs) %*% k + solve(s_fobs)
  s_f <- solve(s_f)
  
  
  ### calculate the adjusted fluxes (optimal)
  f <- f_obs + s_f %*% t(k)%*% solve(s_robs) %*% (r_obs - k %*% f_obs)
  ### now calculate the adjusted storage terms (optimal)
  r <- A %*% f
  
  
  ### optimal uncertainties
  lh_u_o <- sqrt(s_f[1,1])
  h_u_o <- sqrt(s_f[2,2])
  g_u_o <- sqrt(s_f[3,3])
  p_u_o <- sqrt(s_f[4,4])
  q_u_o <- sqrt(s_f[5,5])
  
  ### Chi-squared test value: a measure of the goodness of the fit
  x2 <- t(f-f_obs) %*% solve(s_fobs)%*%(f-f_obs) + t(r-r_obs)%*% solve(s_robs) %*% (r-r_obs)
  
  ### write the optimal fluxes, optimal uncertainties and x2 in a  vector
  output_vector <- c(f,r, lh_u_o, h_u_o, g_u_o, p_u_o, q_u_o, x2)
  
  return(output_vector)
  
}


# define the number of cores you want to use
 ncores <- 2 
 beginCluster(ncores)

for(year in start_year:end_year){

  for(month in 1:12){
    
    
    ## read rasters containting the component fluxes of the water and energy budgets and their uncertainties at a year and a month
    ## All rasters have same grid, and a common mask has been applied to all of them.
    ## All rasters have the same unit W/m2
    
    raster_lh <- raster(paste0("lh_",year,"_",month,".nc")) # latent heat flux
    raster_h <- raster(paste0("h_",year,"_",month,".nc")) # sensible heat flux
    raster_g <- raster(paste0("g_",year,"_",month,".nc")) # ground heat flux
    raster_p <- raster(paste0("p_",year,"_",month,".nc")) # precipitation
    raster_q <- raster(paste0("q_",year,"_",month,".nc")) # runoff
    raster_rn <- raster(paste0("rn_",year,"_",month,".nc")) # net radiation flux
    raster_ds <- raster(paste0("ds_",year,"_",month,".nc")) # change in water storage
    raster_lh_u <- raster(paste0("lh_u_",year,"_",month,".nc")) # uncertainty of latent heat flux
    raster_h_u <- raster(paste0("h_u_",year,"_",month,".nc")) # uncertainty of sensible heat flux
    raster_g_u <- raster(paste0("g_u_",year,"_",month,".nc")) # uncertainty of ground heat flux 
    raster_p_u <- raster(paste0("p_u_",year,"_",month,".nc")) # uncertainty of precipitation
    raster_q_u <- raster(paste0("q_u_",year,"_",month,".nc")) # uncertainty of runoff
    raster_rn_u <- raster(paste0("rn_u_",year,"_",month,".nc")) # uncertainty of net radiation
    raster_ds_u <- raster(paste0("ds_u_",year,"_",month,".nc")) # uncertainty of change in water storage
    
    ## stack rasters -- necessary to apply overlay function
    s <- stack(raster_lh, raster_h, raster_g, raster_p, raster_q, raster_rn, raster_ds, raster_lh_u, raster_h_u, raster_g_u, raster_p_u, raster_q_u, raster_rn_u, raster_ds_u )
    
    
    ### Enforce the closure of the water and energy budgets simultaneously. Apply the function at each grid 
    raster_all_opt <- overlay(s, fun=Solve_Budgets_Grid , unstack=TRUE)
    
   
    ### optimized fluxes
    raster_lh_opt <- raster_all_opt[[1]]
    
    raster_h_opt <- raster_all_opt[[2]]
    
    
    raster_g_opt <- raster_all_opt[[3]]
    
    raster_p_opt <- raster_all_opt[[4]]
    
    raster_q_opt <- raster_all_opt[[5]]
    
    raster_rn_opt <- raster_all_opt[[6]]
    
    raster_ds_opt  <- raster_all_opt[[7]]
    
    ### optimized uncertainties
    
    raster_lh_u_opt <- raster_all_opt[[8]]
    
    raster_h_u_opt <- raster_all_opt[[9]]
   
    raster_g_u_opt <- raster_all_opt[[10]]
    
    raster_p_u_opt <- raster_all_opt[[11]]
    
    raster_q_u_opt <- raster_all_opt[[12]]
    
    
    ### chi_squared test -- goodness of the fit (the fit is good when the values of chi_squared <= 7)
    
    raster_x2_opt <- raster_all_opt[[13]]
    
    
  }
}


endCluster()
setwd("C:/Users/user/Desktop/SpatialR")
library("Matrix")

# calculate matrix of transmission likelihood (consider generation interval)
calc.lik_GI <- function(t, lpdf_GI , cond=1){
  
  t <- as.numeric(t)
  N <- length(t)
  
  
  if(cond==1){
    dts <- pmax(rep(t, times=N) - rep(t, each=N), 0)
    dT <- Matrix(dts, nrow = N, ncol = N, byrow = T)
    dT@x <-  lpdf_GI(dT@x)
    return(dT)
    
  } else if( cond==2){
    
    ind <- c()
    for(i in 2:N) ind <- c(ind, 1:(i-1))
    dts <- rep(t, times=c(1:N)-1) - t[ind] 
    
    ind <- c()
    for(i in 1:(N-1)) ind <- c(ind, (1:i)-1L)
    pind <- c(0L)
    for(i in 1:(N)) pind <- c(pind, pind[i]+i-1L)
    dT <- new("dtCMatrix", i = ind , p = pind, x= lpdf_GI(dts),
              uplo = "U", diag = "N",
              Dim = c(N, N), Dimnames = list(NULL, NULL))
    return(dT)
  }
}


# calculate matrix of transmission likelihood (consider spatial proximity)
calc.lik_SP <- function(x1, x2, lpdf_SP, lon_lat = T ){
  
  N <- length(x1)
  if(N != length(x2)) stop("length differs in x- and y-coordinte")
  
  x1 <- as.numeric(x1)*pi/180
  x2 <- as.numeric(x2)*pi/180
  
  dx1s <- rep(x1, times=N) - rep(x1, each=N)
  x2s_1 <- rep(x2, times=N)
  x2s_2 <-  rep(x2, each=N)
  dx2s <- x2s_1 - x2s_2
  dss <- sin(dx2s/2)^2+ cos(x2s_1)*cos(x2s_2)*sin(dx1s/2)^2
  dss <- 2*6371E3*asin(sqrt(dss)) 
  
  dss[dss==0] <- 1
  ind <- c(1)
  for(i in 2:N) ind <- c(ind, (i-1)*N + c(1:i))
  dss[ind] <- 0
  
  dS <- Matrix(dss, nrow = N, ncol = N, byrow = T)
  dS@x <-  lpdf_SP(dS@x)
  
  return(dS)
}


# calculate individual R and matrix of transmission probability
calc.Rj <- function(t, x1=NULL, x2=NULL, lpdf_GI, lpdf_SP=NULL, adj.sp=T){
  
  dd <- data.frame(t = as.numeric(t))
  dd$count <- 1
  if(adj.sp){ 
    dd$x1 <- as.numeric(x1)*pi/180 
    dd$x2 <- as.numeric(x2)*pi/180
  }
  dd<- (dd[order(dd$t),])
  
  lik_gi <- calc.lik_GI(t = dd$t, lpdf_GI)
  if(adj.sp){ 
    lik_sp <- calc.lik_SP(x1 = dd$x1, x2 = dd$x2, lpdf_SP)
    lik_sp[lik_gi==0] <- 0
    lik_sp <- drop0(lik_sp)
    lik_t <- lik_gi + lik_sp
  } else { 
    lik_t <- lik_gi }
  
  lik_t@x <- exp(lik_t@x)
  
  Pij <- lik_t %*% Diagonal(x = 1 / colSums(lik_t,na.rm = T) )
  
  Rj <- rowSums(Pij,na.rm = T)
  
  return(list(Rj=Rj, Pij=Pij)) 
} 

##########


dir()
load("DATA_Trans_Dengue_target.RData")
df <- dfList[['TN2015']]
df <- df[,c("Date_Onset","Enumeration_unit_long", "Enumeration_unit_lat")]

lpdf_GI <- function(dt) dgamma(dt, 
                               shape=(20^2)/9, scale=9/20, 
                               log = T)
lpdf_SP <- function(dd) dexp(dd,
                             rate = 1/125, log = T)


Res1 <- calc.Rj(t = df$Date_Onset, x1 = df$Enumeration_unit_long, x2 = df$Enumeration_unit_lat,
               lpdf_GI = lpdf_GI, lpdf_SP = lpdf_GI, adj.sp = T)
saveRDS(Res1, "RjEST_TN2015_adj,RDS")

Res0 <- calc.Rj(t = df$Date_Onset, x1 = df$Enumeration_unit_long, x2 = df$Enumeration_unit_lat,
                lpdf_GI = lpdf_GI, lpdf_SP = lpdf_GI, adj.sp = F)
saveRDS(Res0, "RjEST_TN2015_unadj,RDS")


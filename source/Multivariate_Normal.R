library("MASS")

random_variables = 5
n = 100
rep = 1000
params = list()

data = array(matrix(0L, nrow = random_variables, ncol = n), dim = c(random_variables, n, rep))

for (k in 1 : rep){
  
  set.seed(k)
  mu = c(0, 0, 0, 0, 0) 
  sigma2 = (matrix(c(1, 0.5, 0.25, 0.125, 0.0625, 
                     0.5, 1, 0.5, 0.25, 0.125,
                     0.25, 0.5, 1, 0.5, 0.25,
                     0.125, 0.25, 0.5, 1, 0.5,
                     0.0625, 0.125, 0.25, 0.5, 1), ncol = 5))
  
  data[ , , k] = matrix(mvrnorm(n = n, mu = mu, Sigma = sigma2), nrow = random_variables, ncol = n, byrow = T)
  
  f_optim = function(parms, data){
    
    mu = list()
    mu[[1]] = parms[1]
    mu[[2]] = parms[2]
    mu[[3]] = parms[3]
    mu[[4]] = parms[4]
    mu[[5]] = parms[5]
    
    sig = list()
    sig[[1]] = parms[6]
    sig[[2]] = parms[7]
    sig[[3]] = parms[8]
    sig[[4]] = parms[9]
    sig[[5]] = parms[10]
    
    rho  = parms[11]
    
    x = list()
    x[[1]] = data[1,,k]
    x[[2]] = data[2,,k]
    x[[3]] = data[3,,k]
    x[[4]] = data[4,,k]
    x[[5]] = data[5,,k]
    
    soma = 0
    
    for (i in 1:4){
      for (j in (i + 1):5){
        
        func = -n/2 * log(2*pi) - n * log(sig[[j]] * sig[[i]]) - n/2 * log(1 - rho^(2 * abs(j-i))) - 
        1/(2 * (1 - rho^(2 * abs(j-i)))) * (sum(((x[[i]] - mu[[i]]) / sig[[i]])^2) + sum(((x[[j]] - mu[[j]])/sig[[j]])^2) -
        2 * rho^abs(j-i) * sum((x[[i]] - mu[[i]])/sig[[i]] *(x[[j]] - mu[[j]])/sig[[j]]))
        
        soma = func + soma 
        
      }
    }
    
    return(soma)
  }
  
  params[[k]] = optim(c(rep(0.5, 5), rep(0.5, 5), .1), f_optim, 
                      method = "L-BFGS-B", 
                      data = data, control = list(fnscale=-1), 
                      lower = c(rep(-10, 5), rep(.01, 5), -0.95), 
                      upper = c(rep( 10, 5), rep(10, 5), 0.95))$par
}  

test = lapply(params, mean)

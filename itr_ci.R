## Create confidence bounds for ITR


value = function(data, d, tau, h0, h1, gR0, gR1, gA1){
  
  data = transformData(data, 12)
  
  id = data$id
  n = length(unique(id))
  m = as.numeric(data$m)
  K = max(m)
  A = data$A
  
  gA0 = 1 - gA1
  h  = A * h1 + (1 - A) * h0
  gR = A * gR1 + (1 - A) * gR0
  
  gAd = d * gA1 + (1 - d) * gA0
  gRd = d * gR1 + (1 - d) * gR0
  hd = d * h1 + (1 - d) * h0
  Ad = d * A + (1 - d) * (1 - A)
  
  crit = TRUE
  iter = 1
  
  ind = outer(m, 1:K, '<=')
  
  while(crit && iter <= 20){
    
    S1 = tapply(1 - h1, id, cumprod, simplify = FALSE)
    S0 = tapply(1 - h0, id, cumprod, simplify = FALSE)
    G1 = tapply(1 - gR1, id, cumprod, simplify = FALSE)
    G0 = tapply(1 - gR0, id, cumprod, simplify = FALSE)
    
    St1 = do.call('rbind', S1[id])
    St0 = do.call('rbind', S0[id])
    Sm1 = unlist(S1)
    Sm0 = unlist(S0)
    Gm1 = unlist(G1)
    Gm0 = unlist(G0)
    
    Z1 = -rowSums((ind * St1)[, 1:(tau-1)]) / bound(Sm1 * gA1 * Gm1)
    Z0 = -rowSums((ind * St0)[, 1:(tau-1)]) / bound(Sm0 * gA0 * Gm0)
    Zd = d * Z1 + (1 - d) * Z0
    
    eps   = coef(glm(Lm ~ 0 + offset(qlogis(hd)) + Zd,
                      family = binomial(), subset = Im == 1 & Ad == 1,
                      data = data))
    
    h1old = h1
    h0old = h0
    
    h1 = bound01(plogis(qlogis(h1) + eps * Z1))
    h0 = bound01(plogis(qlogis(h0) + eps * Z0))
    hd = d * h1 + (1 - d) * h0
    
    iter =  iter + 1
    
    crit = abs(eps) > 1e-4/n
    
  }
  
  
  S1 = tapply(1 - h1, id, cumprod, simplify = FALSE)
  S0 = tapply(1 - h0, id, cumprod, simplify = FALSE)
  G1 = tapply(1 - gR1, id, cumprod, simplify = FALSE)
  G0 = tapply(1 - gR0, id, cumprod, simplify = FALSE)
  
  St1 = do.call('rbind', S1[id])
  St0 = do.call('rbind', S0[id])
  Std = d * St1 + (1 - d) * St0
  Sm1 = unlist(S1)
  Sm0 = unlist(S0)
  Gm1 = unlist(G1)
  Gm0 = unlist(G0)
  
  Z1 = -rowSums((ind * St1)[, 1:(tau-1)]) / bound(Sm1 * gA1 * Gm1)
  Z0 = -rowSums((ind * St0)[, 1:(tau-1)]) / bound(Sm0 * gA0 * Gm0)
  Zd = d * Z1 + (1 - d) * Z0
  
  DT = with(data, tapply(Im * Ad * Zd * (Lm - hd), id, sum))
  print(mean(DT))
  DWd = with(data, rowSums(Std[m == 1, 1:(tau-1)]))
  theta = 1 + mean(DWd)
  D = DT + DWd - theta
  
  sdn = sqrt(var(D)/n)
  
  out = c(theta, sdn)
  names(out) = c('value', 'sd')
  return(out)
  
}
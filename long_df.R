

### create function that tranforms initial data frame to a long data frame

long_data = function(df){
  n = dim(df)[1]
  K = max(df$time)
  N = n*K
  ## create columns of the long data frame
  id = rep(df$id, each = K)
  # arm = rep(df$arm, each  = K)
  m  = rep(1:K, n)
  Lm = Rm = Im = Jm = rep(NA, N)

  for(t in 1:K){
    Rm[m == t] = (1 - df$censored) * (df$time == t)
    Lm[m == t] = df$censored * (df$time == t)
    Im[m == t] = 1*(df$time >= t)
    Jm[m == t] = (df$time > t) * df$censored + (df$time >= t) * (1 - df$censored)
  }

  long_df = data.frame(id, m, Rm, Lm, Im, Jm)
  return(long_df)
}

#' Make a variogram map with ggplot
#' @return plot
#' @export
ggvariogram_map = function(formula, sp_df, cutoff=NA, width=NA, threshold=NA){
  if(is.na(cutoff)){
    # assume cutoff is data range
    cutoff = mean(sp_df@bbox[3]-sp_df@bbox[1], sp_df@bbox[4]-sp_df@bbox[2])
  }

  if(is.na(width)){
    # assume width is 5% of cutoff (1/20)
    width = cutoff/20
  }

  vmap = variogram(formula, sp_df, cutoff=cutoff, width=width, map=TRUE)

  if(is.na(threshold)){
    # assume threshold is first quartile
    threshold = quantile(vmap$map$np.var1,0.25)
  }

  vmap_df = data.frame(
    semivariance = vmap$map$var1,
    no_pairs = vmap$map$np.var1,
    x = vmap$map@coords[,1],
    y = vmap$map@coords[,2]
  ) %>%
    drop_na() %>%
    filter(no_pairs > threshold)

  ggplot(vmap_df) +
    geom_tile(aes(x=x/1000,y=y/1000,fill=semivariance)) +
    coord_equal() +
    theme_minimal() +
    ylab('Northing (km)') +
    xlab('Easting (km)') +
    theme(legend.position = 'top') +
    scale_fill_viridis(option='cividis')
}

#' Compute covariance for reduced distance h
#' @param type 0:constant, 1:nugget, 2:spherical, 3:exponential, 4:gaussian, else: linear
#' @param h lag distance
#' @param epsilon absolute comparison accuracy
#' @return covariance based on reduced lag distance h
#' @export
#'
cova <- function(type, h, epsilon = 1e-12){
  if (type == 0){
    # Constant value
    return(0*h+1)
  } else if (type == 1){
    # Nugget
    return(as.integer(h<epsilon))
  } else if (type == 2){
    # Spherical model
    return(
      1
      - 1.5*ifelse(h>1,1,h)
      + 0.5*ifelse(h^3>1,1,h^3)
    )
  } else if (type == 3){
    # Exponential
    return(exp(-h))
  } else if (type == 4){
    # Gaussian model
    return(exp(-h^2))
  } else if (type == 5){
    # Linear Model
    return(-h)
  } else {
    warning('Unavailable covariance model, using linear model')
    # Linear Model
    return(-h)
  }
}

#' Setup a rotated matrix of a geostatistical model
#' @param model the geostatistical model matrix
#' @param it the specific rows (the nested model) to rotate
#' @return rotated matrix of the model
#' @export
#'
setrot <- function(model, it){
  deg2rad = pi/180
  ranges = model[it,2:4]
  angles = model[it,5:7]

  # matrix of coordinate reduction
  redmat = matrix(0,nrow=3,ncol=3)
  diag(redmat) = (1./(eps_thresh+ranges))

  a = (90-angles[1])*deg2rad
  b = -angles[2]*deg2rad
  c = angles[3]*deg2rad

  cosa = cos(a)
  sina = sin(a)
  cosb = cos(b)
  sinb = sin(b)
  cosc = cos(c)
  sinc = sin(c)

  rotmat = matrix(0,nrow=3,ncol=3)
  rotmat[1,1] = cosb * cosa
  rotmat[1,2] = cosb * sina
  rotmat[1,3] = -sinb
  rotmat[2,1] = -cosc*sina + sinc*sinb*cosa
  rotmat[2,2] = cosc*cosa + sinc*sinb*sina
  rotmat[2,3] =  sinc * cosb
  rotmat[3,1] = sinc*sina + cosc*sinb*cosa
  rotmat[3,2] = -sinc*cosa + cosc*sinb*sina
  rotmat[3,3] = cosc * cosb

  solve(rotmat)%*%redmat
}

#' Calculate weights for variogram for simulated annealing
#' @param weighting_method 0:constant weight, 1: ~npairs, 2:~/distance 3: ~npairs/distance
#' @param gam the variogram/correlogram
#' @return vector of weights, equal to nrows of variogram
#' @export
#'
calculate_weights <- function(weighting_method, gam){
  if (weighting_method == 0){
    weights = gam[,2]*0+1
  } else if (weighting_method == 1){
    weights = gam[,2]
  } else if (weighting_method == 2){
    weights = 1/(1e-10+gam[,1])
  } else {
    weights = gam[,2]/(1e-10+gam[,1])
  }
  weights/sum(weights)
}

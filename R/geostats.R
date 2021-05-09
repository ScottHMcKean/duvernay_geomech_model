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
#' Ported from Emery at al. 2009
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
#' Ported from Emery at al. 2009
#' @param model the geostatistical model matrix
#' @param it the specific rows (the nested model) to rotate
#' @return rotated matrix of the model
#' @export
#'
setrot <- function(model, it, eps_thresh){
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
#' Ported from Emery at al. 2009
#' Added log methods to soften later portion of curve
#' @param weighting_method 0:constant weight, 1: ~npairs, 2:~/distance 3: ~npairs/distance, 4: log(~npairs), 5:~/log(distance) 6: ~log(npairs)/log(distance)
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
  } else if (weighting_method == 3) {
    weights = gam[,2]/(1e-10+gam[,1])
  } else if (weighting_method == 4){
    weights = log(gam[,2])
  } else if (weighting_method == 5){
    weights = 1/(1e-10+log(gam[,1]))
  } else if (weighting_method == 6){
    weights = log(gam[,2])/(1e-10+log(gam[,1]))
  } else {
    weights = gam[,2]*0+1
  }
  weights/sum(weights)
}

#' Fit Linear Model of Coregionalization
#' Ported from Emergy et al. 2009
#' @param azimuth azimuth for sample variogram/covariance (nvariogram * 1 vector), degrees
#' @param dip dips for sample variogram/covariance (nvariog * 1 vector), degrees
#' @param nlag number of lags for each variogram/correlogram (nvariog * 1 vector)
#' @param tail tail variables (nvariog * 1 vector), an integer vector
#' @param head head variables (nvariog * 1 vector), an integer vector
#' @param gam dataframe of simple and cross variograms/covariances. The shape of the dataframe should be sum(nlag) rows x 3 columns
#' @param model the variogram/covariance model matrix with dimensions of the number of structures (n_structures) * 8
#' @param vartype the script handles traditional variograms (1) or centered covariances (2)
#' @param maxiterations maximum number of consecutive iterations with no accepted transitions
#' @param max_time_s maximum processing time in seconds
#' @param p_acc_0 probability of accepting a non-favorable transition at iteration 0
#' @param p_acc_2k probability of accepting a non-favorable transition at iteration 2000
#' @param gen_vec_corr correlation between generated Gaussian vectors
#' @param eps_thresh absolute comparison threshold for machine accuracy)
#' @param weighting_method 0:constant weight, 1: ~npairs, 2:~/distance 3: ~npairs/distance
#' @return sills of nested structures (n_structures * n_fields^2 matrix) and weighted sum of squares for optimal fit
#' @export
#'
fit_linear_model_of_coregionalization <- function(
  azimuth, dip, nlag, tail, head, gam, model, vartype,
  maxiterations = 5e3, max_time_s = 100, p_acc_0 = 0.9,
  p_acc_2k = 0.1, gen_vec_corr = 0.9, eps_thresh = 1E-12,
  weighting_method = 1
){
  azm = pi/180*azimuth
  dip = pi/180*dip
  directing_vector = matrix(
    c(sin(azm)*cos(dip), cos(azm)*cos(dip), sin(dip)),
    ncol=3
    )

  nvariog = length(azm)

  nlag = c(0,cumsum(nlag))
  n = nlag[nvariog+1]
  n_fields = max(tail,head)
  n_structures = dim(model)[1]
  if (dim(model)[1]<8){model[,8] = 1}
  if (max(model[,1])>11){stop('Error in model type')}

  sills = matrix(0, nrow=n_structures, ncol=n_fields^2)

  # remove zero distances from variogram
  is_zero_dist = (abs(gam[,1])<eps_thresh) & (abs(gam[,3])<eps_thresh)
  gam = data.matrix(gam[!is_zero_dist,])

  # calculate variogram/correlogram weights
  weights = calculate_weights(weighting_method, gam)

  for (l in 1:nvariog){
    lag = gam[(nlag[l]+1):nlag[l+1]]
    lag_mat = matrix(c(lag,lag,lag), ncol=3)
    dim(lag_mat)
    u_mat = matrix(
      rep(directing_vector[l,], dim(lag_mat)[1])
      ,ncol=3, byrow=TRUE
    )
    h_mat = lag_mat * u_mat
    if (l == 1){
      h = rbind(c(0,0,0),h_mat)
    } else {
      h = rbind(h, h_mat)
    }
  }

  g=matrix(0, nrow=dim(h)[1], ncol=n_structures)

  for (i in 1:n_structures){
    rot_model = setrot(model,i, eps_thresh)
    h_adj = h %*% t(rot_model)
    h_adj = sqrt(colSums(t(h_adj)^2))
    C = cova(model[i,1],h_adj)
    g[,i] = C
  }

  if (vartype == 1){
    # traditional variograms
    g = matrix(1, nrow=dim(h)[1]-1, ncol=n_structures)-g[2:(n+1),]
  } else {
    # centered covariance
    g = g[2:n+1,]
  }

  A = cbind(
    diag(n_fields),
    matrix(0, nrow=n_fields,ncol=(n_structures-1)*n_fields)
  )

  for (i in 1:n_structures){
    A_p = A[,((i-1)*n_fields+1):(i*n_fields)]
    B = A_p * t(A_p)
    if (i==1){
      sills = as.vector(B)
    } else {
      sills = rbind(sills, as.vector(B))
    }
  }

  gen_vec_init = A
  gamma = g %*% sills
  WSS = 0
  for (l in 1:nvariog){
    this_gam = gamma[(nlag[l]+1):(nlag[l+1]),]
    gammail = array(this_gam,
      dim = c(nlag[l+1]-nlag[l],n_fields,n_fields)
      )
    gammail = gammail[,tail[l],head[l]]
    WSS = (WSS
      + t(weights[(nlag[l]+1):nlag[l+1]])
      %*% (gam[(nlag[l]+1):nlag[l+1],3]-gammail)^2
    )
  }
  t0 = -WSS/log(p_acc_0)
  alpha = exp(log(-WSS/t0/log(p_acc_2k))/2000)

  t_start = Sys.time()
  rejections = 0
  iterations = 0
  while ((iterations < maxiterations) & ((Sys.time()-t_start) < max_time_s)){
    if (iterations %% 50 == 0){
      print(str_c(iterations,": ", WSS))
    }
    iterations = iterations +1
    temperature = t0*alpha^(iterations-1)
    this_field = ceiling(n_fields*runif(1))
    old_sills = sills
    Aprime = A
    gen_vec = (
      gen_vec_corr*gen_vec_init[this_field,]
      + sqrt(1-gen_vec_corr*gen_vec_corr)*rnorm(n_structures*n_fields)
    )
    Aprime[this_field,] = gen_vec/sqrt(sum(gen_vec^2))
    for (i in 1:n_structures){
      B = (Aprime[,((i-1)*n_fields+1):(i*n_fields)]
           %*% t(Aprime[,((i-1)*n_fields+1):(i*n_fields)])
      )
      old_sills[i,] = t(B)
    }
    gammaprime = g %*% old_sills
    WSSprime = 0
    ## THERE IS A PROBLEM HERE
    for (l in 1:nvariog){
      this_gam = gammaprime[(nlag[l]+1):(nlag[l+1]),]
      gammail = array(this_gam,
        dim = c(nlag[l+1]-nlag[l],n_fields,n_fields)
      )
      gammail = gammail[,tail[l],head[l]]
      WSSprime = (WSSprime
        + t(weights[(nlag[l]+1):nlag[l+1]])
        %*% (gammaprime[(nlag[l]+1):nlag[l+1],3]-gammail)^2
      )
    }
    accept = (runif(1) < exp((WSS-WSSprime)/temperature))
    if (accept == TRUE){
      A = Aprime
      sills = old_sills
      WSS = WSSprime
      rejections = 0
      gen_vec_init[this_field,] = gen_vec
    } else {
      rejections = rejections+1
    }
  }

  results = list(sills, WSS)
  names(results) = c('sills','WSS')
  results
}

#' Helper to map tail id from gstat to LMC code, rowwise implemented
#' @param id id value from gstat
#' @return tail
#' @export
get_tail_id = function(id, data_ids){
  split = strsplit(toString(id),"\\.")[[1]]
  if (length(split) == 1){
    tail = data_ids[split[1]]
  } else {
    tail = data_ids[split[1]]
  }
  tail
}

#' Helper to map head id from gstat to LMC code, rowwise implemented
#' @param id id value from gstat
#' @return head
#' @export
get_head_id = function(id, data_ids){
  split = strsplit(toString(id),"\\.")[[1]]
  if (length(split) == 1){
    head = data_ids[split[1]]
  } else {
    head = data_ids[split[2]]
  }
  head
}

#' Make an input list for the LMC algorithm from gstat
#' @param g gstat object and data
#' @param vmodel common variogram model used to build gstat object
#' @return lmc_input list (named for use in lmc function)
#' @export
make_lmc_input <- function(g, vmodel){
  # make variogram
  vgm = variogram(g, cressie=TRUE)

  # extract ids
  data_ids = seq(1, length(g$data))
  names(data_ids) = names(g$data)

  # extract information for Emery algorithm
  vgm_stats = v %>% group_by(id) %>%
    summarize(
      azimuth = mean(dir.hor),
      dip = mean(dir.ver),
      nlag = n()
    ) %>%
    rowwise() %>%
    mutate(tail = get_tail_id(id, data_ids)) %>%
    mutate(head = get_head_id(id, data_ids))

  # generate variogram with proper column names
  vgm_sa = v %>%
    dplyr::select(lag_dist = dist, npairs = np, variance = gamma)

  lmc_input = list(
    vgm_stats %>% pull(azimuth),
    vgm_stats %>% pull(dip),
    vgm_stats %>% pull(nlag),
    vgm_stats %>% pull(tail),
    vgm_stats %>% pull(head),
    vgm_sa %>% as_tibble(),
    make_lmc_model(vmodel),
    1)

  names(lmc_input) = c(
    'azimuth', 'dip', 'nlag', 'tail', 'head', 'gam', 'model', 'vartype'
  )

  lmc_input
}

#' Make an LMC compliant model from the gstat common variogram model
#' @param vmodel common variogram model used to build gstat object
#' @return lmc compliant model
#' @export
make_lmc_model = function(vmodel){
  mod_cova = sapply(vmodel$model, model_type_to_cova)

  lmc_model = cbind(
    mod_cova,
    vmodel$range,
    vmodel$range*vmodel$anis1,
    vmodel$range*vmodel$anis2,
    vmodel$ang1,
    vmodel$ang2,
    vmodel$ang3,
    c(0,0,0)
  )
  rownames(lmc_model) = vmodel$model
  colnames(lmc_model) = NULL

  lmc_model
}

#' Convert gstat factors to LMC COVA model types, list applied
#' @param mod_factor the factor from gstat for model type
#' @return integer of COVA type
#' @export
model_type_to_cova <- function(mod_factor){
  mod_str = toString(mod_factor)
  case_when(
    mod_str == "Nug" ~ 1,
    mod_str == "Sph" ~ 2,
    mod_str == "Exp" ~ 3,
    mod_str == "Gau" ~ 4,
    mod_str == "Lin" ~ 5,
    TRUE ~ 0
  )
}

#' Assign the LMC sills back to the gstat object
#' @param g gstat object
#' @return sills produced from LMC algorithm
#' @export
assign_lmc_sills = function(g, sills){
  data_ids = seq(1, length(g$data))
  names(data_ids) = names(g$data)
  for (id in names(g$model)){
    head_n = get_head_id(id, data_ids)
    tail_n = get_tail_id(id, data_ids)
    coord = (head_n-1)*length(data_ids)+tail_n
    g$model[[id]]$psill = sills[,coord]
  }
  g
}

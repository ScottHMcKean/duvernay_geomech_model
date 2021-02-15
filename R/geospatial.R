#' Make a 3D plot using Plot3D with a viridis colour scale and sf dataframe
#' @return plot
#' @export
make_3d_plot <- function(sf_df, z_col, col_col, col_func=viridis(50), size=0.1){
  coords = sf_df %>% sf::st_coordinates(.)
  scatter3D(
    x=coords[,1],
    y=coords[,2],
    z=sf_df %>% pull(z_col),
    colvar=sf_df %>% pull(col_col),
    pch=19, cex=size, bty='b2',
    theta=0, phi=30,
    ylab='N', xlab='E', zlab='',
    ticktype = 'detailed',
    col = col_func
  )
}

# make a variogram map with ggplot
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


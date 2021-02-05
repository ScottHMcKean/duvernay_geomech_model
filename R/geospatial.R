#' Make a 3D plot using Plot3D with a viridis colour scale and sf dataframe
#' @param filepath
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

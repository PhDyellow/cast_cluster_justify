make_sim_mat <- function(equal_points, point_sigma) {
#Take points and sd info, then return the similarity matrix, where similarity is the Bhattacharyya coefficient
  row_pairs <- expand.grid(i = seq_along(equal_points$x), j = seq_along(equal_points$x))
  dist_long <- purrr::map2_dbl(row_pairs$i, row_pairs$j,
              ~ {
                if(.x > .y) {
                    b_dist <- my_bhattacharyya_dist(
                                   t(equal_points[.x,]),
                                   point_sigma$sigma[[.x]],
                                   point_sigma$sigma_det[[.x]],
                                   t(equal_points[.y,]),
                                   point_sigma$sigma[[.y]],
                                   point_sigma$sigma_det[[.y]]
                                 )
                  return(b_dist)
                } else {
                  return(NA)
                }
              }, predicted_stats)
  ## Bhattacharyya coefficient = exp(-Bhattacharyya distance)
  sim_mat <- matrix(exp(-dist_long), nrow(equal_points), nrow(equal_points))
  sim_mat[upper.tri(sim_mat)] <- t(sim_mat)[upper.tri(sim_mat)]
  diag(sim_mat) <- 1

return(sim_mat)
}

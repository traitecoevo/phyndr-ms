venn <- function() {
  cols <- venn_cols()

  spp <- list(diversity = 350000, # Will's estimate
              accepted  = 312069, # From TaxonLookup
              genbank   = 92704,  # Scrubbed genbank
              traits    = 40159,  # Zanne et al
              overlap   = 28868)

  border <- 0.025
  grid.newpage()
  vp_d <- sviewport(1 - 2 * border, x=unit(border, "snpc"), just="left")
  pushViewport(vp_d)
  grid.circle(gp=gpar(col=cols$diversity, fill=cols$diversity))

  vp_a <- sviewport(sqrt(spp$accepted) / sqrt(spp$diversity), .49, .51)
  pushViewport(vp_a)
  grid.circle(gp=gpar(col=cols$accepted, fill=cols$accepted))

  ## Radii of the two circles:
  r_g <- (sqrt(spp$genbank) / sqrt(spp$accepted)) / 2
  r_t <- (sqrt(spp$traits)  / sqrt(spp$accepted)) / 2

  ## This is the area that we *want*:
  area <- spp$overlap / spp$genbank
  target <- spp$overlap / spp$genbank * pi * r_g^2
  d <- uniroot(function(x) lens_area(x, r_g, r_t) - target, c(0.1, 0.3))$root

  x_g <- 0.4
  y <- 0.6
  theta <- 40

  vp_venn <- viewport(angle=theta)
  pushViewport(vp_venn)
  grid.circle(x_g, y, r_g, gp=gpar(col=cols$genbank, fill=cols$genbank))
  grid.circle(x_g + d, y, r_t, gp=gpar(col=cols$traits, fill=cols$traits))

  lens_xy <- lens(x_g, x_g + d, y, r_g, r_t)
  grid.polygon(lens_xy[, "x"], lens_xy[, "y"],
               gp=gpar(col=cols$overlap, fill=cols$overlap))

  popViewport() # venn
  popViewport() # accepted
  popViewport() # border
  pushViewport(viewport(x=unit(border, "snpc"), y=border,
                        height=1 - 2 * border,
                        width=1 - 2 * border,
                        just=c("left", "bottom")))

  lab <- c("Total diversity",
           "Accepted names",
           "Trait data",
           "Overlap",
           "Genetic data")

  lab_y <- seq(0, 1, length.out=length(lab) + 2)
  lab_y <- lab_y[-c(1, length(lab_y))]
  lab_x <- unit(1, "npc") - unit(border, "snpc")

  grid.text(lab, lab_x, lab_y, just="right")
}

venn_cols <- function() {
  cols_flatui <- flatui_colour_scheme()
  list(diversity = cols_flatui$grey_2,
       accepted  = cols_flatui$yellow,
       genbank   = cols_flatui$orange,
       traits    = cols_flatui$blue,
       overlap   = cols_flatui$purple_dk)
}

circular_segment_area <- function(r, d) {
  r^2 * acos(d / r) - d * sqrt(r^2 - d^2)
}

lens_area <- function(d, r_g, r_t) {
  x <- (d^2 - r_t^2 + r_g^2) / (2 * d)
  chord_len <- sqrt(r_g^2 - x^2)
  circular_segment_area(r_g, r_g - abs(x - r_g)) +
    circular_segment_area(r_t, r_t - abs(x - (d - r_t)))
}

lens <- function(x_g, x_t, y, r_g, r_t) {
  d <- abs(x_g - x_t)
  d_x <- (d^2 - r_t^2 + r_g^2) / (2 * d)
  chord_len <- sqrt(r_g^2 - d_x^2)

  p0 <- c(x_g + d_x, y + chord_len)
  p1 <- c(x_g + d_x, y - chord_len)

  t_g0 <- atan2(p0[2] - y, p0[1] - x_g)
  t_g1 <- atan2(p1[2] - y, p1[1] - x_g)

  t_t0 <- pi - atan2(y - p0[2], x_t - p0[1])
  t_t1 <- pi - atan2(y - p1[2], x_t - p1[1])

  t_g <- seq(t_g0, t_g1, length.out=100)
  t_t <- seq(t_t0, t_t1, length.out=100)
  lens_g <- cbind(x=x_g + r_g * cos(t_g), y=y + r_g * sin(t_g))
  lens_t <- cbind(x=x_t + r_t * cos(t_t), y=y + r_t * sin(t_t))
  rbind(lens_g, lens_t)
}

sviewport <- function(p, ...) {
  viewport(height=unit(p, "snpc"), width=unit(p, "snpc"), ...)
}

#' @title PrepPalette
#'
#' @description Given a palette name and data frame with clustering or grouping
#' information in an 'ident' column, return a palette containing a color for each
#' unique cluster identity.
#'
#' @param palette_use The name of a palette to use.  Must be a palette available in
#' the Paletteer package.  If there are more unique identities than colors in the
#' palette, additional values will be created by interpolation.
#' @param bins Number of discrete colors to generate from palette.
#'
#' @return A list containing color values.
#' @export
#'
#' @import paletteer
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#'
PrepPalette <- function(palette_use,
                        bins) {

  if (palette_use %in% palettes_d_names$palette) {
    color.package <- palettes_d_names$package[which(palette_use == palettes_d_names$palette)]
    pal <- paletteer_d(
      package = !!color.package,
      palette = !!palette_use
    )
  } else if (palette_use %in% palettes_c_names$palette) {
    color.package <- palettes_c_names$package[which(palette_use == palettes_c_names$palette)]
    pal <- paletteer_c(
      package = !!color.package,
      palette = !!palette_use
    )
  } else {
    pal <- palette_use
  }
  pal <- colorRampPalette(pal)(bins)
  return(pal)
}

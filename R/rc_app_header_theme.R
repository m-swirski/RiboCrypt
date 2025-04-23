rc_header_image <- function() {
  tags$head(
    tags$link(rel = "icon",
              href = file.path("images", "favicon.png"),
              type = "image/x-icon"))
}

rc_theme <- function() {
  bslib::bs_theme(
    version = 5,
    primary = "#6dbaff", secondary = "#ff7e7e",
    success = "#c0ffa4", font_scale = 1.2, bootswatch = "zephyr")
}

rc_title <- function() {
  withTags(
    a(img(src = file.path("images", "logo_traces_update.png"),
          alt = "RiboCrypt",
          height = 50)))
}

rc_header_styling <- function() {
  tags$style(HTML("
  /* Reduce height of navbar */
  .navbar {
    min-height: 50px !important;
    padding-top: 5px !important;
    padding-bottom: 5px !important;
  }

  .navbar-brand {
    padding-top: 5px !important;
    padding-bottom: 5px !important;
    font-size: 18px;
  }

  .navbar-nav > li > a {
    padding-top: 10px !important;
    padding-bottom: 10px !important;
  }
"))
}




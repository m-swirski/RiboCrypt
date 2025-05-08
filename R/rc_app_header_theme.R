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

browser_ui_settings_style <- function() {
  tags$style(HTML("
      [id$='floating_settings'] {
        position: absolute;
        top: 80px;
        left: 20px;
        width: 350px;
        z-index: 1000;
        background-color: rgba(255, 255, 255, 0.95);
        border: 1px solid #ddd;
        padding: 15px;
        border-radius: 10px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.2);
      }

      .floating_settings_panel.hidden {
        display: none !important;
      }

      .floating_settings_panel {
        max-width: 700px;
        max-height: 85vh;
        overflow-y: auto;
      }

      .tab-content > .tab-pane {
        padding: 10px !important;
      }

      .form-group {
        margin-bottom: 10px !important;
      }

      #clip {
        background-color: orange;
      }

      #settingsCollapse .panel-title {
        font-size: 0.75em;  /* Adjust this value to reduce the font size */
      }

    .selectize-input {
      font-size: 15px; /* Shrinks the font inside the selectize input box */
      height: 30px; /* Shrinks the height of the input box */
    }
    .selectize-dropdown {
      font-size: 15px; /* Shrinks the font inside the dropdown list */
    }
    .selectize-input .selectize-input-inner {
      padding: 2px; /* Adjust padding to fit the box */
    }
    .selectize-input .selectize-input-default {
      padding: 2px; /* Adjust padding inside input field */
    }
    .shiny-input-container:has(.selectize-control) > label {
      font-size: 17px;  /* Only affects labels for selectizeInput */
    }
    "))
}




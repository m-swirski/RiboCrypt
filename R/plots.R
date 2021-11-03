createSinglePlot <- function(display_range, reads, withFrames, colors, kmers = 1,
                             kmers_type = "mean", ylabels, lines, type = "lines"){
  count <- NULL # Avoid data.table warning
  if (withFrames) {

    if (type %in% c("stacks", "area" )) {
      profile <- getStackProfile(display_range, reads, kmers, kmers_type = kmers_type)
    } else profile <- getRiboProfile(display_range, reads, kmers, kmers_type = kmers_type)
  } else {
    profile <- getCoverageProfile(display_range, reads, kmers, kmers_type = kmers_type)
  }
  profile_plot <- ggplot(profile)
  if (length(lines) > 0) profile_plot <- profile_plot + geom_vline(xintercept = lines, col = names(lines), linetype = 4)
  profile_plot <- profile_plot + geom_hline(yintercept = 0)
  if (!withFrames) {
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position), fill = colors, col = "black", size = 0.1, alpha = 0.8, position = "identity")
  } else if (type == "lines") {
  profile_plot <- profile_plot +
    geom_line(aes(y = count, x = position,color = frame), size = 0.5)
  } else if (type == "stacks") {
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position,fill = frame), size = 0.1, alpha = 0.8, col = "black")
    } else if (type == "columns") {
    profile_plot <- profile_plot +
      geom_col(aes(y = count, x = position, fill = frame))
    } else if (type == "area") {
      profile_plot <- profile_plot +
        geom_area(aes(y = count, x = position,fill = frame), size = 0.1, alpha = 0.8, col = "black", position = 'identity')
  }
  profile_plot <- profile_plot +
    theme(legend.position = "none") +
    ylab(ylabels) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 8, face = "bold"),
          axis.text.y = element_text(size = 6)) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    scale_x_continuous(expand = c(0,0))
  profile_plot <- automateTicksRNA(profile_plot)

  profile_plot
}


#profiles <- mapply(function(x,y,z) getProfileAnimate(display_range, x, y, z), reads, withFrames, kmers,  SIMPLIFY = FALSE)

getPlotAnimate <- function(profile, withFrames, colors, ylabels, lines){
  count <- NULL # Avoid data.table warning
  if (withFrames) {
    profile_plot <- ggplot(profile)
    if (length(lines) > 0) profile_plot <- profile_plot + geom_vline(xintercept = lines, col = names(lines), linetype = 4)
    profile_plot <- profile_plot +
      geom_line(aes(y = count, x = position,color = frame, frame = file), size = 0.75) +
      theme(legend.position = "none") +
      ylab(ylabels) +
      xlab("position [nt]") +
      theme(plot.margin = unit(c(0,0,0,0), "pt")) +
      scale_x_continuous(expand = c(0,0)) +
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
    profile_plot <- automateTicks(profile_plot)

  } else {
    profile_plot <- ggplot(profile)
    if (length(lines) > 0) profile_plot <- profile_plot + geom_vline(xintercept = lines, col = names(lines), linetype = 4)
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position, frame = file), fill = colors, position = "identity") +
      ylab(ylabels) +
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()) +
      theme(plot.margin = unit(c(0,0,0,0), "pt"))+
      scale_x_continuous(expand = c(0,0))
    profile_plot <- automateTicksRNA(profile_plot)
  }
  profile_plot
}

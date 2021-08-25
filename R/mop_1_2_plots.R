createSinglePlot <- function(target_range,reads,withFrames, colors, kmer = 1, ylabels, lines, type = "lines"){
  if (withFrames) {

    if (type %in% c("stacks", "area" )) {
      profile <- getStackProfile(target_range, reads, kmer)

    } else profile <- getRiboProfile(target_range, reads, kmer)
  } else {
    profile <- coveragePerTiling(target_range, subsetByOverlaps(reads, target_range), as.data.table = TRUE)

  }
    profile_plot <- ggplot(profile)
    if (length(lines) > 0) profile_plot <- profile_plot + geom_vline(xintercept = lines, col = "black", linetype = 4)
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
          geom_area(aes(y = count, x = position,fill = frame),position = 'identity', size = 0.1, alpha = 0.8, col = "black")
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


#profiles <- mapply(function(x,y,z) getProfileAnimate(target_range, x, y, z), reads, withFrames, kmers,  SIMPLIFY = FALSE)

getPlotAnimate <- function(profile,withFrames, colors, ylabels, lines){
  if (withFrames) {
    profile_plot <- ggplot(profile)
    if (length(lines) > 0) profile_plot <- profile_plot + geom_vline(xintercept = lines, col = "black", linetype = 4)
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
    if (length(lines) > 0) profile_plot <- profile_plot + geom_vline(xintercept = lines, col = "black", linetype = 4)
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

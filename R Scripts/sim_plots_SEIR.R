library(raster)
library(rgdal)
library(tidyverse)
library(RANN)
library(openxlsx)
library(ggmap)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)

sim_all_seir <- readRDS("sim_all_seir.rds")
sim_no_public_seir <- readRDS("sim_no_public_seir.rds")
sim_only_essential_seir <- readRDS("sim_only_essential_seir.rds")
sim_lockdown_40_seir <- readRDS("sim_lockdown_40_seir.rds")
sim_lockdown_60_seir <- readRDS("sim_lockdown_60_seir.rds")
sim_lockdown_40_stop_seir <- readRDS("sim_lockdown_40_stop_seir.rds")

tokyo <- raster("japan_map/greater_tokyo_3km.tif")

# Plot Simulation Curves --------------------------------------------------

# plotting the infection curves
library(plotly)
library(plyr)
Sys.setenv("plotly_username"="")
Sys.setenv("plotly_api_key"="")

# plot(susceptible_pop_norm)
# plot(infected_pop_norm)

plot_curves <- function(sim_output, plot_name, file_name){
  
  dfPlot <- sim_output$curves
  # "#bc5090", "#003f5c", "#ffa600"
  
  margin <- list(autoexpand = T,
                 l = 10,
                 r = 130,
                 t = 50)
  
  Susceptible <- list(
    xref = 'paper',
    x = 1,
    y = tail(dfPlot$s_p, 1) * 1,
    xanchor = 'left',
    yanchor = 'middle',
    text = paste('Susceptible ', round(tail(dfPlot$s_p, 1) * 100, 1), '%'),
    font = list(size = 12,
                color = '#ffa600'),
    showarrow = FALSE)
  
  Infected <- list(
    xref = 'paper',
    x = 1,
    y = tail(dfPlot$i_p, 1) * 1,
    xanchor = 'left',
    yanchor = 'middle',
    text = paste('Infected ', round(tail(dfPlot$i_p, 1) * 100, 1), '%'),
    font = list(size = 12,
                color = '#bc5090'),
    showarrow = FALSE)
  
  Recovered <- list(
    xref = 'paper',
    x = 1,
    y = tail(dfPlot$r_p, 1) * 1,
    xanchor = 'left',
    yanchor = 'middle',
    text = paste('Recovered ', round(tail(dfPlot$r_p, 1) * 100, 1), '%'),
    font = list(size = 12,
                color = '#003f5c'),
    showarrow = FALSE)
  
  peak_day <- dfPlot$day[dfPlot$i_p == max(dfPlot$i_p)]
  
  Peak_Line <- list(
    type = "line",
    line = list(color = "rgba(58, 71, 80, 0.6)"),
    xref = 'x',
    yref = 'y',
    x0 = peak_day,
    y0 = 0,
    x1 = peak_day,
    y1 = 1)
  
  Peak <- list(
    xref = 'x',
    x = peak_day + 1,
    y = 0.95,
    xanchor = 'left',
    yanchor = 'middle',
    text = paste('Day: ', peak_day, "\nInfected: ", round(dfPlot$i_p[peak_day] * 100, 1), "%"),
    font = list(size = 12,
                color = "rgba(58, 71, 80, 0.6)"),
    showarrow = FALSE)
  
  p <- plot_ly(data = dfPlot, x= ~day) %>%
    add_lines(y = ~s_p*1, name = "Susceptible", line = list(shape = "spline", width = 4, color = '#ffa600')) %>%
    add_lines(y = ~i_p*1, name = "Infected", line = list(shape = "spline", width = 4, color = '#bc5090')) %>%
    add_lines(y = ~r_p*1, name = "Recovered", line = list(shape = "spline", width = 4, color = '#003f5c')) %>%
    
    layout(yaxis = list(title = '% Population', tickformat = "%", range = c(0, 1.05)),
           xaxis = list(title = 'Day'),
           title = plot_name,
           legend = list(x = 0.03, y = -0.5, orientation = 'h', bgcolor = 'rgba(0,0,0,0)'),
           autosize = FALSE,
           margin = margin,
           showlegend = FALSE,
           paper_bgcolor = "transparent",
           plot_bgcolor = "transparent") %>%
    
    layout(annotations = Susceptible) %>%
    layout(annotations = Infected) %>%
    layout(annotations = Recovered) %>%
    layout(shapes  = Peak_Line) %>%
    layout(annotations = Peak)
  p
  
  api_create(p, filename = file_name)
  
}

plot_curves(sim_output = sim_all, 
            plot_name = "Simulation of coronavirus outbreak in central Tokyo area \n(all traffic as usual)", 
            file_name = "tokyo_sim_all_traffic_curves")

plot_infected_curves <- function(sim_output1, sim_output2, sim_output3, sim_output4,
                                 plot_name, file_name, label1, label2, label3, label4,
                                 hospital_curve = TRUE, 
                                 hospital_max_capacity = 0.02, hospital_ratio = .2){
  
  na.pad <- function(x,len){
    x[1:len]
  }
  
  makePaddedDataFrame <- function(l,...){
    maxlen <- max(sapply(l,length))
    data.frame(lapply(l,na.pad,len=maxlen),...)
  }
  
  percent <- function(x, digits = 4, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  
  # set up plot data
  dfPlot1 <- sim_output1$curves
  dfPlot2 <- sim_output2$curves
  dfPlot3 <- sim_output3$curves
  dfPlot4 <- sim_output4$curves
  
  max_day <- max(dfPlot1$day, dfPlot2$day, dfPlot3$day, dfPlot4$day)
  
  if (hospital_curve == FALSE) {
    dfPlot <- makePaddedDataFrame(list(day = 1:max_day,
                                       i_p1 = dfPlot1$i_p,
                                       i_p2 = dfPlot2$i_p,
                                       i_p3 = dfPlot3$i_p,
                                       i_p4 = dfPlot4$i_p))
  } else {
    dfPlot <- makePaddedDataFrame(list(day = 1:max_day,
                                       i_p1 = dfPlot1$i_p * hospital_ratio,
                                       i_p2 = dfPlot2$i_p * hospital_ratio,
                                       i_p3 = dfPlot3$i_p * hospital_ratio,
                                       i_p4 = dfPlot4$i_p * hospital_ratio))
  }
  
  # "#bc5090", "#003f5c", "#ffa600"
  
  # set up plot parameters
  margin <- list(autoexpand = T,
                 l = 10,
                 r = 10,
                 t = 50)
  
  peak_day1 <- dfPlot1$day[dfPlot1$i_p == max(dfPlot1$i_p, na.rm = T)]
  peak_day2 <- dfPlot2$day[dfPlot2$i_p == max(dfPlot2$i_p, na.rm = T)]
  peak_day3 <- dfPlot3$day[dfPlot3$i_p == max(dfPlot3$i_p, na.rm = T)]
  peak_day4 <- dfPlot4$day[dfPlot4$i_p == max(dfPlot4$i_p, na.rm = T)]
  
  peak_value1 <- dfPlot$i_p1[peak_day1]
  peak_value2 <- dfPlot$i_p2[peak_day2]
  peak_value3 <- dfPlot$i_p3[peak_day3]
  peak_value4 <- dfPlot$i_p4[peak_day4]
  
  max_value <- max(peak_value1, peak_value2, peak_value3, peak_value4)
  # print(peak_value1)
  
  x_min <- 0
  x_max <- plyr::round_any(max_day, 10, ceiling)
  
  y_min <- 0
  y_max <- plyr::round_any(max_value, 0.02, ceiling)
  
  # set up annotations (e.g. lines and text)
  Peak_Line1 <- list(
    type = "line",
    line = list(color = "rgba(58, 71, 80, 0.3)"),
    xref = 'x',
    yref = 'y',
    x0 = peak_day1,
    y0 = 0,
    x1 = peak_day1,
    y1 = peak_value1)
  
  Peak1 <- list(
    xref = 'x',
    x = peak_day1 + 1,
    y = peak_value1,
    xanchor = 'left',
    yanchor = 'middle',
    text = paste0(label1, "\n", "Day: ", peak_day1, " ", 
                  percent(peak_value1, 1)),
    font = list(size = 12,
                color = "rgba(58, 71, 80, 0.6)"),
    showarrow = FALSE)
  
  Peak_Line2 <- list(
    type = "line",
    line = list(color = "rgba(58, 71, 80, 0.3)"),
    xref = 'x',
    yref = 'y',
    x0 = peak_day2,
    y0 = 0,
    x1 = peak_day2,
    y1 = peak_value2)
  
  Peak2 <- list(
    xref = 'x',
    x = peak_day2 + 1,
    y = peak_value2,
    xanchor = 'left',
    yanchor = 'middle',
    text = paste0(label2, "\n", "Day: ", peak_day2, " ", 
                  percent(peak_value2, 1)),
    font = list(size = 12,
                color = "rgba(58, 71, 80, 0.6)"),
    showarrow = FALSE)
  
  Peak_Line3 <- list(
    type = "line",
    line = list(color = "rgba(58, 71, 80, 0.3)"),
    xref = 'x',
    yref = 'y',
    x0 = peak_day3,
    y0 = 0,
    x1 = peak_day3,
    y1 = peak_value3)
  
  Peak3 <- list(
    xref = 'x',
    x = peak_day3 + 1,
    y = peak_value3,
    xanchor = 'left',
    yanchor = 'middle',
    text = paste0(label3, "\n", "Day: ", peak_day3, " ", 
                  percent(peak_value3, 1)),
    font = list(size = 12,
                color = "rgba(58, 71, 80, 0.6)"),
    showarrow = FALSE)
  
  Peak_Line4 <- list(
    type = "line",
    line = list(color = "rgba(58, 71, 80, 0.3)"),
    xref = 'x',
    yref = 'y',
    x0 = peak_day4,
    y0 = 0,
    x1 = peak_day4,
    y1 = peak_value4)
  
  Peak4 <- list(
    xref = 'x',
    x = peak_day4 + 1,
    y = peak_value4,
    xanchor = 'left',
    yanchor = 'middle',
    text = paste0(label4, "\n", "Day: ", peak_day4, " ", 
                  percent(peak_value4, 1)),
    font = list(size = 12,
                color = "rgba(58, 71, 80, 0.6)"),
    showarrow = FALSE)
  
  Capacity_Line <- list(
    type = "line",
    line = list(color = "rgba(58, 71, 80, 1)", dash = 'dash'),
    xref = 'x',
    yref = 'y',
    x0 = 0,
    y0 = hospital_max_capacity,
    x1 = x_max,
    y1 = hospital_max_capacity)
  
  Capacity <- list(
    xref = 'x',
    yref = 'y',
    x = x_max / 2,
    y = hospital_max_capacity + 0.004,
    xanchor = 'left',
    yanchor = 'middle',
    text = "Supposed maximum capacity of\nthe city's health care system",
    font = list(size = 13,
                color = "rgba(58, 71, 80, 1)"),
    showarrow = FALSE)
  
  # plot graph
  p <- plot_ly(data = dfPlot, x= ~day) %>%
    add_lines(y = ~ i_p1, name = "All traffic as usual", 
              line = list(shape = "spline", width = 1, color = "rgba(152, 0, 67, 1)"), 
              fill = "tozeroy", fillcolor = "rgba(152, 0, 67, 0.5)") %>%
    add_lines(y = ~ i_p2, name = "No public transport", 
              line = list(shape = "spline", width = 1, color = "rgba(221, 28, 119, 1)"), 
              fill = "tozeroy", fillcolor = "rgba(221, 28, 119, 0.5)") %>%
    add_lines(y = ~ i_p3, name = "Only essential travel", 
              line = list(shape = "spline", width = 1, color = "rgba(223, 101, 176, 1)"), 
              fill = "tozeroy", fillcolor = "rgba(223, 101, 176, 0.5)") %>%
    add_lines(y = ~ i_p4, name = "Lockdown", 
              line = list(shape = "spline", width = 1, color = "rgba(215, 181, 216, 1)"), 
              fill = "tozeroy", fillcolor = "rgba(215, 181, 216, 0.5)") %>%
    
    layout(yaxis = list(title = '% Population', tickformat = ".1%", range = c(y_min, y_max)),
           xaxis = list(title = 'Day', range = c(x_min, x_max)),
           title = plot_name,
           # legend = list(x = 0.03, y = -0.15, orientation = 'h', bgcolor = 'rgba(0,0,0,0)'),
           # autosize = FALSE,
           margin = margin,
           showlegend = FALSE,
           paper_bgcolor = "transparent",
           plot_bgcolor = "transparent") %>%
    
    layout(
      # shapes = list(Peak_Line1, Peak_Line2, Peak_Line3, Peak_Line4),
      annotations = list(Peak1, Peak2, Peak3, Peak4))
  
  if (hospital_curve == TRUE) {
    p <- p %>%
      layout(shapes = list(Capacity_Line),
             annotations = list(Capacity))
  }
  
  # upload to plotly
  api_create(p, filename = file_name)
  
}

plot_infected_curves(sim_all_seir, sim_no_public_seir, sim_only_essential_seir, sim_lockdown_40_seir,
                     plot_name = "Infection curves under different scenarios", 
                     file_name = "tokyo_sim_comparison_curves_seir", 
                     label1 = "All as usual", label2 = "No public transport", 
                     label3 = "Only essential travel", label4 = "Lockdown", 
                     hospital_curve = F, hospital_max_capacity = 0.02, hospital_ratio = 0.2)

# Plot Simulations --------------------------------------------------------

generate_sim_visual <- function(plot_path, raster_map, sim_output, every_n_day = 2,
                                dot_scale = 1000, dot_size = .7, dot_alpha = .7,
                                output_dpi = 200, 
                                plot_title = NULL){
  
  generate_dots <- function(ras, dot_scale = 200) {
    # given a raster, generate dots in each raster cell that represent the value
    # of that raster cell
    dx <- diff(c(xmin(ras), xmax(ras))) / ncol(ras) / 1.5 # 2/3 of horizontal width
    dy <- diff(c(ymin(ras), ymax(ras))) / nrow(ras) / 1.5 # 2/3 of vertical width
    r <- sqrt(dx^2 + dy^2)
    xy_mat <- coordinates(ras)                            # 2-column matrix of coordinates of cells centroid
    n_cell <- nrow(xy_mat)                                # number of cells
    n_dot_vec <- (values(ras) / dot_scale) %>% 
      ceiling %>% 
      ifelse(is.na(.), 0, .)                              # a vector of number of dots for each cell
    if (sum(n_dot_vec) >= 1) {
      dfDots <- data.frame(dot_id = 1:sum(n_dot_vec), x = NA, y = NA)
      dot_id <- 0
      for (i in 1:n_cell){
        x <- xy_mat[i, "x"]
        y <- xy_mat[i, "y"]
        n_dot <- n_dot_vec[i]
        if (n_dot > 0) {
          for (d in 1:n_dot){
            dot_id <- dot_id + 1
            dfDots[dot_id, "x"] <- x + runif(1, -r, r)
            dfDots[dot_id, "y"] <- y + runif(1, -r, r)
          }
        }
      }
    } else {
      dfDots <- NULL
    }
    return(dfDots)
  }
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  percent <- function(x, digits = 4, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  
  # debug
  # dot_scale <- 100
  # plot_path <- "sim_all_seir/"
  # raster_map <- tokyo
  # sim_output <- sim_no_public_seir
  # every_n_day <- 5
  # dot_scale <- 1000
  # dot_size <- .7
  # dot_alpha <- .7
  # output_dpi <- 200
  # day <- 1
  
  # set plot dir
  if (substrRight(plot_path, 1) != "/") {
    plot_path <- paste0(plot_path, "/")
  }
  dir.create(file.path(getwd(), plot_path), showWarnings = FALSE)
  
  # read base map of the city area
  
  map_box <- raster_map@extent %>% as.matrix() # get the coordinates of the area
  base_map <- get_stamenmap(map_box, zoom = 10, maptype = "terrain-lines", color = "bw")
  
  # generate plot
  df_pop_norm <- sim_output$curves
  infected_by_day <- sim_output$infected_by_day
  
  total_days <- df_pop_norm$day
  days <- total_days[seq(1, length(total_days), every_n_day)]
  plot_id <- 0
  
  for (day in days) {
    
    print(paste0("plotting day ", day, "/", max(df_pop_norm$day)))
    measure_traffic <- ifelse(day >= sim_output$reduced_OD_start &
                                day <= sim_output$reduced_OD_stop, "Reduced Traffic", "")
    measure_beta <- ifelse(day >= sim_output$reduced_beta_start &
                             day <= sim_output$reduced_beta_stop, "Reduced Transmission", "")
    plot_id <- plot_id + 1
    infected_vec <- infected_by_day[[day]]
    values(raster_map) <- infected_vec
    
    n_infected <- round(sum(infected_vec))
    i_p <- percent(df_pop_norm[day, "i_p"], 1)
    i_p_c <- percent(1 - df_pop_norm$s_p[day], 1)
    
    # creating dots for plots 
    # create 1 dot for every x cases in the raster cell
    dfDots <- generate_dots(raster_map, dot_scale)
    
    g_map <- ggmap(base_map, extent = "device", darken = 0) +
      annotate("text", x = 140.75100, y = 35.37100, hjust = 1,
               label = paste0("Day ", day), colour = "#484848", size = 10) +
      annotate("text", x = 140.75100, y = 35.32100, hjust = 1, 
               label = paste0("Infected: ", i_p), colour = "#484848", size = 4) +
      annotate("text", x = 140.75100, y = 35.29100, hjust = 1,
               label = paste0("Infected (cumulative): ", i_p_c), 
               colour = "#484848", size = 4) +
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = .5, lineheight = .5, size = 15, colour= "#484848"))
    
    
    if (!is.null(dfDots)) {
      g_map <- g_map + 
        geom_point(data = dfDots,
                   aes(x = x,
                       y = y),
                   size = dot_size,
                   alpha = dot_alpha,
                   color = "#bc5090",
                   shape = 19, 
                   fill = "#bc5090")
    }
    
    # creating the statistic plot at the corner
    dfLine <- df_pop_norm[1:day, ] %>% 
      reshape2::melt(., id.vars = "day", variable.name = "group") %>%
      dplyr::mutate(group = fct_relevel(group, "r_p", "s_p", "e_p", "i_p")) %>%
      dplyr::mutate(group = recode(group, `r_p` = "R", `s_p` = "S", `e_p` = "E", `i_p` = "I"))
    
    g_line <- ggplot(dfLine, aes(x = day, y = value, fill = group)) + 
      geom_area(alpha=0.6 , size=.5, position = "stack") +
      scale_fill_manual(values=c("#003f5c", "#ffa600", "#7a5195", "#ef5675")) +
      ggplot2::xlim(1, day) +
      coord_cartesian(ylim=c(0,1)) +
      
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            legend.position="bottom",
            legend.title = element_blank(),
            legend.text = element_text(colour= "#484848", 
                                       size= 10))
    
    g <- 
      g_map +
      ggmap::inset(ggplotGrob(g_line), 
                   xmin = 140.25100, xmax = 140.87516, 
                   ymin = 34.97162, ymax = 35.27100)
    
    if (measure_traffic != "" | measure_beta != "") {
      g <- g +
        theme(panel.border = element_rect(colour = "red", fill=NA, size=2)) +
        annotate("text", x = 138.9487, y = 36.20625, hjust = 0,
                 label = paste(measure_traffic, measure_beta, sep = "\n"),
                 colour = "red", size = 4)
    }
    
    ggsave(file = paste0(plot_path, str_pad(plot_id, 3, "left", "0"), ".jpg"), 
           plot = g, width = 10, height = 10, dpi = output_dpi)
  }
  print(paste("Completed. Generated", plot_id, "plot(s) in total.", sep = " "))
}

# plot_title <- "Simulation of coronavirus outbreak in the Greater Tokyo Area\n
#               (by databentobox.com; for demostration only)"

generate_sim_visual(plot_path = "sim_all_seir", 
                    raster_map = tokyo, 
                    sim_output = sim_all_seir, 
                    dot_scale = 1000, dot_size = 1.2, dot_alpha = 0.5, 
                    plot_title = "Simulation of coronavirus outbreak in the Greater Tokyo Area\n
              Scenario: All as usual; by databentobox.com")

generate_sim_visual(plot_path = "sim_no_public_seir", 
                    raster_map = tokyo, 
                    sim_output = sim_no_public_seir, 
                    dot_scale = 1000, dot_size = 1.2, dot_alpha = 0.5, 
                    plot_title = "Simulation of coronavirus outbreak in the Greater Tokyo Area\n
              Scenario: No public transport; by databentobox.com")

generate_sim_visual(plot_path = "sim_only_essential_seir", 
                    raster_map = tokyo, 
                    sim_output = sim_only_essential_seir, 
                    dot_scale = 1000, dot_size = 1.2, dot_alpha = 0.5, 
                    plot_title = "Simulation of coronavirus outbreak in the Greater Tokyo Area\n
              Scenario: Only essential travel; by databentobox.com")

generate_sim_visual(plot_path = "sim_lockdown_40_seir", 
                    raster_map = tokyo, 
                    sim_output = sim_lockdown_40_seir, 
                    dot_scale = 1000, dot_size = 1.2, dot_alpha = 0.5, 
                    plot_title = "Simulation of coronavirus outbreak in the Greater Tokyo Area\n
              Scenario: Lockdown from Day 40; by databentobox.com")

generate_sim_visual(plot_path = "sim_lockdown_60_seir", 
                    raster_map = tokyo, 
                    sim_output = sim_lockdown_60_seir, 
                    dot_scale = 1000, dot_size = 1.2, dot_alpha = 0.5, 
                    plot_title = "Simulation of coronavirus outbreak in the Greater Tokyo Area\n
              Scenario: Lockdown from Day 60; by databentobox.com")

generate_sim_visual(plot_path = "sim_lockdown_40_stop_seir", 
                    raster_map = tokyo, 
                    sim_output = sim_lockdown_40_stop_seir, 
                    dot_scale = 1000, dot_size = 1.2, dot_alpha = 0.5, 
                    plot_title = "Simulation of coronavirus outbreak in the Greater Tokyo Area\n
              Scenario: Lockdown from Day 40 to Day 60; by databentobox.com")

# Create Annimations ------------------------------------------------------

# List those Plots, Read them in, and then make animation
# list.files(path = "plots/", pattern = "*.png", full.names = T) %>% 
#   map(image_read) %>% # reads each path file
#   image_join() %>% # joins image
#   image_animate(fps=2) %>% # animates, can opt for number of loops
#   image_write("scenario.gif") # write to current dir

# use ffmepg
## ffmpeg -f image2 -framerate 9 -i %003d.jpg out.gif
## ffmpeg -framerate 9 -i %3d.jpg -y all.mp4

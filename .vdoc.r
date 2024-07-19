#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: sample-season-missing

missing_prob <- function(data, mu=NULL, alpha = 5, theta = 0.3){
                              x_ <- seq_len(nrow(data))
                              mu_ <- ifelse(is.null(mu), length(x_) / 2, mu)
                              std_ <- sqrt(mu_ / theta)
                              y_ <- abs((alpha * exp((-(x_ - mu_)^2) / std_^2)) - alpha) + alpha
                              yn_ <- y_ / (sum(y_))
                        return(yn)
                  }

sample_missing <- function(data, propMissing = 0.25){
            missing.prob <- data.table::data.table()
            for(i in data[, unique(years)]){
                  for(j in data[, unique(site_id)]){
                  missing.prob <- rbind(missing.prob, missing_prob(data[years == i & site_id == j, ]))
                  }
            }
            missing.week <- data[sample(seq_len(.N), round(propMissing * .N), prob = unlist(missing.prob)), ]  
      return(missing.week)
}
missing_prob_vector <- data.table()

for(i in btfl_week_smpl[, unique(years)]){
      for(j in btfl_week_smpl[, unique(site_id)]){
            missing_prob_vector <- rbind(missing_prob_vector, missing_prob(data  = btfl_week_smpl[years == i & site_id == j, ]))
      }
}

btfl_week_missing <- btfl_week_smpl[sample(seq_len(.N), round(0.25*.N), prob = unlist(missing_prob_vector)), ]
btfl_week_missing <- sample_missing(data = btfl_week_smpl, propMissing = 0.25)

btfl_fig3 <- ggplot() +
                geom_point(data=btfl_week_smpl, aes(x=doy, y=count, colour = "count")) + 
                geom_point(data=btfl_week_missing, aes(x=doy, y=count, colour = "missing"), 
                            shape=4, size=2, stroke=2) + 
                geom_line(data = btfl_ts,
                aes(x = doy, y = act, colour = "activity")) +
                xlim(1,365) + ylim(0, max(btfl_ts$count, btfl_ts$act)) + 
                scale_colour_manual("", 
                      breaks = c("count", "activity", "missing"),
                      values = c(cnt_col, flc_col, missing_col)) +
                theme_light() + 
                theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8)) +
                labs(title = paste0("Simulated butterfly counts (", y,")"),
                     subtitle = "- weekly visit (random resampled)",
                     x = "Day of Year",
                     y = "Count")
btfl_fig3
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

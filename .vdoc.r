#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
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

sample_missing <- function(data, quantile_breaks = c(0.20, 0.35, 0.65, 0.80), quantile_probs = c(0.25, 0.5, 0.75)){
      prob_vct <- rep(quantile_probs[1], data[,.N])
      prob_vct[round(quantile(seq_len(data[,.N]), quantile_breaks[1])):round(quantile(seq_len(data[,.N]), quantile_breaks[2]))] <- quantile_probs[2]
      prob_vct[round(quantile(seq_len(data[,.N]), quantile_breaks[3])):round(quantile(seq_len(data[,.N]), quantile_breaks[4]))] <- quantile_probs[2]
      prob_vct[1:round(quantile(seq_len(data[,.N]), quantile_breaks[1]))] <- quantile_probs[3]
      prob_vct[round(quantile(seq_len(data[,.N]), quantile_breaks[4])):data[,.N]] <- quantile_probs[3]
      btfl_week_missing <- data[sample(seq_len(.N), round(0.25*.N), prob = prob_vct/sum(prob_vct)), ]
  return(btfl_week_missing)
}

btfl_week_missing <- data.table()

for(i in btfl_week_smpl[, unique(years)]){
      for(j in btfl_week_smpl[, unique(site_id)]){
            btfl_week_missing <- rbind(btfl_week_missing, sample_missing(data  = btfl_week_smpl[years == i & site_id == j, ]))
      }
}

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

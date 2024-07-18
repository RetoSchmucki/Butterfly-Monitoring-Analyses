#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
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

btfl_week_smpl2 <- btfl_week_smpl[month(date) %in% c(4:9), ]
prob_vct <- rep(0.5,btfl_week_smpl2[,.N])
prob_vct[1:round(quantile(seq_len(btfl_week_smpl2[,.N]), c(0.25)))] <- 0.75
prob_vct[round(quantile(seq_len(btfl_week_smpl2[,.N]), c(0.75))):btfl_week_smpl2[,.N]] <- 0.75

btfl_week_missing <- btfl_week_smpl2[-sample(.N, round(0.75*.N), prob = prob_vct), ]

btfl_fig3 <- ggplot() +
                geom_point(data=btfl_week_smpl2, aes(x=doy, y=count, colour = "count")) + 
                geom_point(data=btfl_week_missing, aes(x=doy, y=count, colour = "missing"), 
                            shape=4, size=2, stroke=2) + 
                geom_line(data = btfl_dt2,
                aes(x = doy, y = act, colour = "activity")) +
                xlim(1,365) + ylim(0, max(btfl_dt1$count, btfl_dt2$act)) + 
                scale_colour_manual("", 
                      breaks = c("count", "activity", "missing"),
                      values = c(cnt_col, flc_col, missing_col)) +
                theme_light() + 
                theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8)) +
                labs(title = "Simulated butterfly counts",
                     subtitle = "- weekly visit (random resampled)")
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

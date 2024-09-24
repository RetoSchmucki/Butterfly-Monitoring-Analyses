## Generalized Additive Model approach

```{R}
#| label: setup_gam
#| results: hide
#| messages: false
#| warning: false

set.seed(13276)

if(!require("rbms")) devtools::install_github("RetoSchmucki/rbms")
if(!require("butterflyGamSims")) devtools::install_github("cbedwards/butterflyGamSims")
if(!require("mixR")) devtools::install_github("RetoSchmucki/mixR") # small fix in plot function, not yet implemented in original and CRAN version

library(rbms)
library(mixR)
library(butterflyGamSims)
library(data.table)
library(ggplot2)

flc_col <- '#ff8c00'
cnt_col <- '#008b8b'
missing_col <- '#8b0000'
GAM_col <- '#483d8b'

## local functions
sim2bms <- function(data, yearKeep = NULL, weeklySample = FALSE, weekdayKeep = NULL, monitoringSeason = NULL){
                
                  btfl_ts <- data.table::data.table(data)[, site_id := paste0("site_",sim.id)]
                  if(!is.null(yearKeep)){
                        btfl_ts <- btfl_ts[years %in% yearKeep, ]
                  }
                  btfl_ts[ , date := as.Date(doy, origin = paste0(years, "-01-01"))-1]
                  btfl_ts[ , week := isoweek(date)]
                  btfl_ts[month(date) != 1 | week < 50, weekday := rowid(week), by = .(site_id, years)]
                  if(isTRUE(weeklySample)){
                        if(!is.null(weekdayKeep)){
                              btfl_ts <- btfl_ts[weekday %in% weekdayKeep, ]
                        }
                        btfl_ts <- btfl_ts[btfl_ts[,.I[sample(.N, 1)], by = .(week, site_id, years)][["V1"]],]
                  }
                  if(!is.null(monitoringSeason)){
                  btfl_ts <- btfl_ts[month(date) %in% monitoringSeason, ]
                  }
            return(btfl_ts)
            }

missing_prob <- function(data, mu=NULL, alpha = 5, theta = 0.3){
                              x_ <- seq_len(nrow(data))
                              mu_ <- ifelse(is.null(mu), length(x_) / 2, mu)
                              std_ <- sqrt(mu_ / theta)
                              y_ <- abs((alpha * exp((-(x_ - mu_)^2) / std_^2)) - alpha) + alpha
                              yn_ <- y_ / (sum(y_))
                        return(yn_)
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
```

## Simulated count

```{R}
#| label: bivoltine
#| results: hide
#| messages: false
#| warning: false

size1 <- 500
peak1 <- 155
sd1 <- 10
size2 <- 250
peak2 <- 230
sd2 <- 10
y <- c(2023)

set.seed(13276)
btfl_data1 <- timeseries_sim(nsims=1,
               year = y,
               doy.samples = seq(from=1, to=365, by=1),
               abund.type = "exp",
               activity.type = "gauss",
               sample.type = "pois",
               sim.parms = list(growth.rate = 0,
                                init.size = size1, 
                                act.mean = peak1,
                                act.sd = sd1)
               )

set.seed(13276)
btfl_data2 <- timeseries_sim(nsims=1,
               year = y,
               doy.samples = seq(from=1, to=365, by=1),
               abund.type = "exp",
               activity.type = "gauss",
               sample.type = "pois",
               sim.parms = list(growth.rate = 0,
                                init.size = size2, 
                                act.mean = peak2,
                                act.sd = sd2)
               )

btfl_data_b2 <- rbind(btfl_data1$timeseries[,c("years", "doy", "count", "act", "sim.id")], 
                      btfl_data2$timeseries[,c("years", "doy", "count", "act", "sim.id")])
btfl_data_b2 <- unique(data.table(btfl_data_b2)[, ":="(count=sum(count), act=sum(act)), by = .(years, doy)])
btfl_data_b2[, abund.true:= sum(c(size1, size2))]

btfl_ts <- sim2bms(data = btfl_data_b2, yearKeep = y)

btfl_fig1 <- ggplot() +
                geom_point(data=btfl_ts, aes(x=doy, y=count, colour = "count")) + 
                geom_line(data = btfl_ts,
                aes(x = doy, y = act, colour = "activity")) +
                xlim(1,365) + ylim(0, max(btfl_ts$count, btfl_ts$act)) + 
                scale_colour_manual("", 
                      breaks = c("count", "activity"),
                      values = c(cnt_col, flc_col)) +
                theme_light() + 
                theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8)) +
                labs(title = paste0("Simulated butterfly counts (", y,")"),
                     subtitle = "- daily visit",
                     x = "Day of Year",
                     y = "Count")
btfl_fig1
## ===========================
## degradation1: weekly count
## ===========================

btfl_week_smpl <- sim2bms(data = btfl_data_b2, yearKeep = y, 
                  weeklySample = TRUE,
                  weekdayKeep = c(2:5),
                  monitoringSeason = c(4:9))

btfl_fig2 <- ggplot() +
                geom_point(data=btfl_week_smpl, aes(x=doy, y=count, colour = "count")) + 
                geom_line(data = btfl_ts,
                aes(x = doy, y = act, colour = "activity")) +
                xlim(1,365) + ylim(0, max(btfl_ts$count, btfl_ts$act)) + 
                scale_colour_manual("", 
                      breaks = c("count", "activity"),
                      values = c(cnt_col, flc_col)) +
                theme_light() + 
                theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8)) +
                labs(title = paste0("Simulated butterfly counts (", y,")"),
                     subtitle = "- weekly visit (random resampled)",
                     x = "Day of Year",
                     y = "Count")

btfl_fig2
## ==============================
## degradation: add missing week
## ==============================

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
```

The flight curve underlying the observed butterfly counts can be seen as a function of time, with counts increasing as adult individuals emerge and decreasing as they die. In the case of a bivoltine species, we would expect two emergence phases over the year, one for each generation. The shape of this bimodal function can be represented by a spline derived from the GAM using the `flight_curve` function in the rbms package. By organising the simulated data into a table of monitoring visits and a table of counts, we can generate a time series of counts for the 2023 season.

The `ts_dwmy_table` function builds a *time series* where every *day*, *week*, *month* and *year* is informed for the defined period (`InitYear` and `LastYear`). We can then inform the monitoring season over the time series, informing the start and end month, if we want to include additional zeros (*anchor*) at the start and end of the monitoring season, and how strong the *anchor* should be. Here we will define a monitoring season that starts on April (4) 1 and ends on the last day of September (9). We will add anchors of two zeros positioned two days before and after the monitoring season. As the time series is defined on a daily (d) basis, the anchors and the position (lag) are in the same units.

```{R}
#| label: bivoltine_flightcurve1
#| results: hide
#| messages: false
#| warning: false

visit_sim <- btfl_week_smpl[!date %in% btfl_week_missing$date, .(site_id, date, count)]
count_sim <- visit_sim[count>=1,][, species := "sp1"]

names(visit_sim) <- toupper(names(visit_sim))
names(count_sim) <- toupper(names(count_sim))

ts_date <- rbms::ts_dwmy_table(InitYear = 2023, LastYear = 2023, WeekDay1 = 'monday')

ts_season <- rbms::ts_monit_season(ts_date,
                       StartMonth = 4,
                       EndMonth = 9, 
                       StartDay = 1,
                       EndDay = NULL,
                       CompltSeason = TRUE,
                       Anchor = TRUE,
                       AnchorLength = 2,
                       AnchorLag = 2,
                       TimeUnit = 'd')

ts_season_visit <- rbms::ts_monit_site(ts_season, visit_sim)

ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, count_sim, sp = "sp1")
```

We can now use the `flight_curve` function of the rbms package to fit a GAM and get the spline corresponding to the smoothed function of the number of butterflies recorded over the time series, days. Here we will fit the model on data simulated for one site over one year (2023). A site is only included if it has been visited at least five time over the season and had 3 visit with positive occurrence. Although this is not very relevant for this dataset, these filters can be useful to discard some poorly monitored site. Remember that the purpose of this model is to inform about temporal change in abundance related to species' phenology and that this information is better in more complete time series and where the species and change in abundance has been adequately detected over the season.

When plotting the resulting flight curve, we can choose to plot the standardised form (NM), where the area under the curve is adjusted to sum to 1, or the curve converted to the count scale by multiplying the standardised curve by the total abundance. Here we plot the flight curve on the count scale.

```{R}
#| label: bivoltine_flightcurve2
#| results: hide
#| messages: false
#| warning: false
### Fitting a GAM to Butterfly Counts

ts_flight_curve <- rbms::flight_curve(ts_season_count, 
                       NbrSample = 300,
                       MinVisit = 5,
                       MinOccur = 3,
                       MinNbrSite = 1,
                       MaxTrial = 4,
                       GamFamily = 'nb',
                       SpeedGam = FALSE,
                       CompltSeason = TRUE,
                       SelectYear = NULL,
                       TimeUnit = 'd')


## plot-fitted-flight-curve

pheno <- ts_flight_curve$pheno

btfl_fig4 <- ggplot() +
                geom_point(data=ts_season_count[ANCHOR == 0 & !is.na(COUNT), ], aes(x=DAY_SINCE, y=COUNT, colour = "count")) +
                geom_point(data=btfl_week_missing, aes(x=doy, y=count, colour = "missing"), 
                            shape=4, size=2, stroke=2) + 
                geom_line(data = btfl_ts, aes(x = doy, y = act, colour = "activity")) +
                geom_line(data = pheno,
                    aes(x = trimDAYNO, y = btfl_ts[, unique(abund.true)]*NM, colour = "GAM_fit")) +
                xlim(1,365) + ylim(0, max(btfl_ts$act, 
                                          pheno$NM*btfl_ts[,unique(abund.true)], 
                                          btfl_week_missing$count, 
                                          ts_season_count[!is.na(COUNT), COUNT] )) + 
                scale_colour_manual("", 
                      breaks = c("count", "activity", "missing", "GAM_fit"),
                      values = c(cnt_col, flc_col, missing_col, GAM_col)) +
                theme_light() + 
                theme(legend.position = "inside", legend.position.inside = c(0.9, 0.8)) +
                labs(title = paste0("Simulated butterfly counts (", y,")"),
                     subtitle = "- Fitting GAM model with rbms",
                     x = "Day of Year",
                     y = "Count")
btfl_fig4
```

## GAM basis for flight curve spline

The spline that represents the flight curve is build from addition of individual basis functions. Although the individual basis functions have no direct interpretation, it may be helpful to plot them to understand the building block of the spline and the mechanics of the Generalised Additive Model (GAM). To do this, we will use the R package `gratia`, which contains useful functions for exploring, evaluating and plotting GAMs. First we will retrieve the model, which is stored in the ts_flight_curve object under the slot, model. This is a list of annual flight curve models (GAM). Here we use the 1st model in the list (we only have one year, 2023).  From this model we can retrieve the summary of the model and the number of knots used for the spline.

```{R}
#| label: gam_basis

if(!require("gratia")) install.packages("gratia")
library("gratia")
library("colorspace")

## retrieve GAM model
mod <- ts_flight_curve$model[[1]]
summary(mod)
mod$smooth[[1]]$bs.dim # knots for the spline

## retrieve basis for the smooth trimDAYNO
x2_bs <- data.table(gratia::basis(mod, select = "s(trimDAYNO)", data = ts_flight_curve$data[, c("trimDAYNO")]))
x2_spl <- x2_bs[, spline := sum(.value), by = trimDAYNO]

# plot GAM basis
ggplot() +
    geom_line(data= x2_bs, aes(x = trimDAYNO, y = .value, colour = .bf, group = .bf)) +
    geom_line(data = x2_spl, aes(x = trimDAYNO, y = spline, colour = "Spline"), lty =2) +
    labs(y = expression(f(trimDAYNO)), x = "DAYNO") +
    scale_colour_manual("Basis functions", breaks =  c(1:9,"Spline"), values = c(qualitative_hcl(9, palette = "Dark 3"), 'black')) +
    theme(legend.key = element_blank())
```

You might think that this (black dashed line) does not really look like the flight curve estimated and returned by the `rbms::flight_curve` function. One thing to remember is that the model is fitted using the log-link, which means that the spline is estimated on the log scale and so are the bases that underlie the spline. The spline is the sum of the 9 basis functions. To demonstrate this, we can extract the spline that is the sum of the bases and then transform the spline to the count scale (i.e. exp()). Since we assume that the count outside the monitoring season is zero, we will multiply the result by the binary variable (0,1) `M_SEASON`. We will then normalise the spline so that the area under the curve sums to one.

:::{.callout-note}
Note that the individual bases have no direct interpretation in the model other than as part of the building blocks of the spline, so their contribution must be considered in the context of the other bases. Also, it is the sum of the bases, the spline, that is optimised in the likelihood function. While the spline can be transformed back to the count scale, the basis should not be, as the spline is the sum of the basis in logarithmic rather than exponential scale.
:::

```{R}
#| label: gam_spline_sum_basis

# retrieve the flight curve from the spline
m_data_sp <- unique(merge(pheno, x2_spl[, .(trimDAYNO, spline)], by= "trimDAYNO"))
m_data_sp[, spline_nm := exp(spline) * M_SEASON][, spline_nm := spline_nm/sum(m_data_sp$spline_nm)]

# plot scaled spline and GAM
ggplot() +
    geom_line(data= m_data_sp, aes(x = trimDAYNO, y = spline_nm, colour = "Spline"), linewidth = 1.5) +
    geom_line(data= pheno, aes(x = trimDAYNO, y = NM, colour = "Flight Curve")) +
    labs(y = expression(f(trimDAYNO)), x = "DAYNO") 
```

#### Incidence analysis ####

# load incidence data

cattleDatInc <- read.csv("FOI_analysis/data/cattleIncidenceDataWithAedesPred.csv")
cattleDatInc <- read.csv("FOI_analysis/data/cattleIncidenceDataWithCulexPred.csv") |> 
  dplyr::select(newID, meanCulex) %>%
  left_join(cattleDatInc, by = "newID") |>
  mutate(culexAbundance.sds = scale(meanCulex),
    aedesAbundance.sds = scale(meanAedes)) 

crs.km <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
cattleDatInc <- st_as_sf(cattleDatInc, coords = c("Longitude", "Latitude"), crs = 4326) |>
  st_transform(crs.km)

#' District isn't in the incidence CSVs, so bring it in from cattle.age
#' (one District per Village, so this is a safe 1:1 join -- verified no
#' Village maps to >1 District and no animal is left unmatched)
cattleDatInc <- cattleDatInc %>%
  left_join(cattle.age %>% distinct(Village, District), by = "Village")

# make a map of where the seroconversions have happened.

villageSeroconversions <-
  cattleDatInc %>%
  mutate(X = st_coordinates(geometry)[, 1],
         Y = st_coordinates(geometry)[, 2]) %>%
  st_drop_geometry() %>%
  group_by(Village) %>%
  summarise(X = mean(X), Y = mean(Y), n = sum(exposed == 1, na.rm = TRUE)) %>%
  st_as_sf(coords = c("X", "Y"), crs = st_crs(cattleDatInc), remove = FALSE)

ggplot(villageSeroconversions) +
  annotation_map_tile(type = "cartolight", zoomin = 0) +  # You can change the tile type
  geom_sf(aes(fill = n),
          size = 3, shape = 21, colour = "darkgrey") +
  scale_fill_viridis_c( option = "G",
    direction = -1, breaks = scales::breaks_width(2)) +
  annotation_scale(location = "tl", width_hint = 0.5) +
  coord_sf(crs = crs.km, expand = TRUE) +
  theme_minimal() +
  theme(
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    plot.margin = margin(2, 2, 2, 2),
    legend.position = "bottom"
  ) +
  labs(fill = "n seroconversions")

ggsave("FOI_analysis/nseroConversions.pdf")

#'Fit serocatalytic model to incidence data for seronegative animals re-sampled in 2024
#'Does vector abundance explain the 2024 seroconversion in light of El Nino 

mods$inputs$spatial24 <- mods$inputs$spatial
mods$inputs$spatial24$obs <-
  bru_obs(
    formula = 
      exposed ~ 
      Intercept +
      culexAbundance.sds + 
      aedesAbundance.sds +
      Village + 
      hhID + 
      field + 
      offset(log(prop_year)),
    family = "binomial",
    control.family = list(link = "cloglog"),
    data = cattleDatInc)


#' Fit model
mods$outputs$spatial24 <- 
  bru(
    components = mods$inputs$spatial24$components,
    obs = mods$inputs$spatial24$obs,
    options = list(
      control.fixed = mods$priors$prior.fixed,
      control.predictor = list(compute = TRUE),
      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, 
                             config = TRUE,
                             return.marginals.predictor = TRUE)))
#summary(model2.inla)
summary(mods$outputs$spatial24)
exp(mods$outputs$spatial24$summary.fixed)

# visualise this 
fit.summary.fixed <- mods$outputs$spatial24$summary.fixed
fit.summary.fixed$ID <- rownames(fit.summary.fixed)
fit.summary.fixed$ID <- factor(fit.summary.fixed$ID, fit.summary.fixed$ID)

# Plot
vectorIncidence <- 
  fit.summary.fixed %>%
  mutate(ID = fct_recode(ID
                         ,Culex = "culexAbundance.sds"
                         ,Aedes = "aedesAbundance.sds"
  )) %>%
  ggplot(aes(y = ID, x = exp(mean))) +
  geom_pointrange(aes(xmin = exp(`0.025quant`), xmax = exp(`0.975quant`))) +
  geom_vline(xintercept = 1, linetype = 2) +
   scale_y_discrete(labels = c(
     "Intercept"
     ,expression(italic("Culex"))
     ,expression(italic("Aedes"))
   )) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")+
  labs(y = "", x = "Estimate \u00B1 95% CI")

#### spatial covariates instead of VA ####
#' 
#' mods$inputs$spatial24env <- mods$inputs$spatial
#' mods$inputs$spatial24env$obs <-
#'   bru_obs(
#'     formula = 
#'       exposed ~ 
#'       Intercept +
#'       ndvi.sds + 
#'       ndvi.sds2 + 
#'       buildingDensity500m.log10.sds +
#'       popDens500m.log10.sds +
#'       proportionIrrigatedCrop.sds +
#'       NRainDays.log10.sds +
#'       elevation500m.sds +
#'       Village + 
#'       hhID + 
#'       field + 
#'       offset(log(prop_year)),
#'     family = "binomial",
#'     control.family = list(link = "cloglog"),
#'     data = cattleDatInc)
#' 
#' 
#' #' Fit model
#' mods$outputs$spatial24env <- 
#'   bru(
#'     components = mods$inputs$spatial24env$components,
#'     obs = mods$inputs$spatial24env$obs,
#'     options = list(
#'       control.fixed = mods$priors$prior.fixed,
#'       control.predictor = list(compute = TRUE),
#'       control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, 
#'                              config = TRUE,
#'                              return.marginals.predictor = TRUE)))
#' #summary(model2.inla)
#' summary(mods$outputs$spatial24env)
#' exp(mods$outputs$spatial24env$summary.fixed)
#' 
#' # visualise this 
#' fit.summary.fixed <- mods$outputs$spatial24env$summary.fixed
#' fit.summary.fixed$ID <- rownames(fit.summary.fixed)
#' fit.summary.fixed$ID <- factor(fit.summary.fixed$ID, fit.summary.fixed$ID)
#' fit.summary.fixed <-
#'   fit.summary.fixed %>%
#'   mutate(ID = fct_recode(ID,
#'                        NDVI = "ndvi.sds",
#'                        `NDVI²` = "ndvi.sds2",
#'                        `proportion \n irrigated` = "proportionIrrigatedCrop.sds",
#'                        elevation = "elevation500m.sds",
#'                        `population density` = "popDens500m.log10.sds",
#'                        `building density` = "buildingDensity500m.log10.sds"#, 
#'                        #`rainy days` = "NRainDays.log10.sds"
#'                        )
#'          )
#' 
#' # Plot
#' spatial24env.plot <- 
#'   ggplot(fit.summary.fixed, aes(y = ID, x = exp(mean))) +
#'   geom_pointrange(aes(xmin = exp(`0.025quant`), xmax = exp(`0.975quant`))) +
#'   geom_vline(xintercept = 1, linetype = 2) +
#'   theme_minimal(base_size = 12) +
#'   labs(y = "", x = "Estimate \u00B1 95% CI")
#' 
#' sapply(mods$outputs, function(f) f$waic$waic)
#' 
#' 
#' #### check in glmer
#' 
#' library(lme4)
#' 
#' test <- glmer(exposed ~ (1|Village) + (1|hhID) + ndvi.sds + ndvi.sds2 + NRainDays.log10.sds + 
#'         elevation500m.sds + buildingDensity500m.log10.sds + popDens500m.log10.sds + 
#'         proportionIrrigatedCrop.sds + offset(log(prop_year)), 
#'         family = binomial(link = "cloglog"), data = cattleDatInc)

###### Null spatial fixed effects model #####
mods$inputs$spatial24.null <- mods$inputs$spatial
mods$inputs$spatial24.null$obs <-
  bru_obs(
    formula = 
      exposed ~ 
      Intercept +
      Village + hhID + field + offset(log(prop_year)),
    family = "binomial",
    control.family = list(link = "cloglog"),
    data = cattleDatInc)


#' Fit spatial model
mods$outputs$spatial24.null <- 
  bru(
    components = mods$inputs$spatial24.null$components,
    obs = mods$inputs$spatial24.null$obs,
    options = list(
      control.fixed = mods$priors$prior.fixed,
      control.predictor = list(compute = TRUE),
      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, 
                             config = TRUE,
                             return.marginals.predictor = TRUE)))

summary(mods$outputs$spatial24.null)

exp(mods$outputs$spatial24.null$summary.fixed)

sapply(mods$outputs, function(f) f$waic$waic)


#############################################
#####  Apply CV to incidence modelling  #####
#############################################

#' Store cluster ID in dat.all.inc for now (don't want to do anything silly to the main df)
#' Assign and plot clusters
#' cv_cluster() uses k-means internally, so fix the seed - otherwise the fold
#' assignment (and therefore the CV comparison between models) is different
#' every time this is re-run
cluster.k <- 5
set.seed(456)
cv <- cv_cluster(x = dat, k = cluster.k)
cv_plot(cv, dat)

#' Store cluster ID in dat.all for now (don't want to do anything silly to the main df)
dat.all.inc <- dat
dat.all.inc$fold <- NA
dat.all.inc$fold[1:nrow(cv$biomod_table)] <- 
  ((1 - cv$biomod_table) %*% 1:cluster.k)[, 1]
table(dat.all.inc$fold)

inla.opt <- 
  list(
    num.threads = "1:1",
    inla.mode = "experimental",
    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, 
                           config = TRUE,
                           return.marginals.predictor = TRUE,
                           openmp.strategy = NULL),
    control.family = list(link = "cloglog"))

inla.opt.cv <- inla.opt
inla.opt.cv$control.compute$cpo <- inla.opt.cv$control.compute$dic <-
  inla.opt.cv$control.compute$waic <- FALSE
#inla.opt.cv$control.compute$config <- FALSE
inla.opt.cv$control.compute$return.marginals.predictor <- FALSE
inla.opt.cv$num.threads <- detectCores()
inla.opt.cv$inla.mode <- "compact"
inla.opt.cv$control.inla = list(strategy = "simplified.laplace", int.strategy = "eb", h = 0.02)


modelSpecsIncidence <-
  list(spatial24 = "culexAbundance.sds + aedesAbundance.sds + offset(log(exposed.yrs))")
#modelSpecsIncidence$spatial24env <-
  #"ndvi.sds + ndvi.sds2 +  elevation500m.sds + buildingDensity500m.log10.sds + popDens500m.log10.sds + proportionIrrigatedCrop.sds  + offset(log(exposed.yrs))"
modelSpecsIncidence$spatial24.null <-
  "offset(log(exposed.yrs))"


#' Loop over each fold, fitting model to training folds and predicting
#' FOI on held-out fold

#' Storage for the CV diagnostic checks below: per-model out-of-fold
#' predictions and per-fold Brier scores, so we can check (1) whether models
#' actually make different predictions and (2) whether the Brier score
#' ranking is consistent across folds rather than driven by one or two.
cv.diagnostics <- list()

cv.results <-
  sapply(names(mods$outputs)[4:length(names(mods$outputs))], function(fit.name) {
    
    dat.all.inc$pred.fold <- NA
    dat.all.inc$fold.diameter <- NA
    
    print(fit.name)
    fit <- mods$outputs[[fit.name]]
    
    #' Original components
    components.cv <- mods$inputs[[fit.name]]$components 
    
    #' Loop over all k clusters
    for(i in 1:cluster.k) {
      print(i)
      print(table(!is.na(dat.all.inc$pred.fold)))
      
      # Calculate distance from centroid of fold to nearest point outside fold
      # get approximate diameter for the distance from the fold centroid to the nearest point outside the fold (scale of fold)
      fold.centroid <- colMeans(st_coordinates(dat.all.inc[dat.all.inc$fold %in% i, ]))
      dat.all.inc$fold.diameter[dat.all.inc$fold %in% i] <- 
        2 * min(st_distance(st_sfc(st_point(fold.centroid), crs = crs.km), 
                            dat.all.inc[dat.all.inc$fold %in% (1:ncol(cv$biomod_table))[-i], ]))
      print(unique(dat.all.inc$fold.diameter[dat.all.inc$fold %in% i])) 
      
      # Refit model to 4 of the 5 folds, leaving 1 fold out (fold i)
      inla.opt.cv$verbose <- inla.opt$verbose <- FALSE
      dat.all.cv <- dat.all.inc
      dat.all.cv$exposed24[dat.all.cv$fold %in% i] <- NA

      # Specify the observation model

      obs.cv <-
          bru_obs(
            formula = formula(paste0("exposed24 ~ Intercept + Village + hhID + field +", modelSpecsIncidence[which(names(modelSpecsIncidence)==fit.name)])),
            family = "binomial",
            data = dat.all.cv)
    

      # Fit the model
      cv.fit.start <- Sys.time()
      fit.cv <- 
        bru(
          components = components.cv,
          obs = obs.cv,
          options = inla.opt)    
      cv.fit.finish <- Sys.time()
      print(cv.fit.finish - cv.fit.start)
      

      dat.all.inc$pred.fold[which(!cv$biomod_table[, i])] <-
        predict(fit.cv,
              formula = paste0("Intercept + ", modelSpecsIncidence[[fit.name]]),
              newdata = dat.all.cv)$mean[which(!cv$biomod_table[, i])]
      print(table(!is.na(dat.all.inc$pred.fold)))
  }
  
    #' Back-transform from the cloglog-scale linear predictor to a probability
    dat.all.inc$probPred <- 1 - exp(-exp(dat.all.inc$pred.fold))
    plot(dat.all.inc$probPred, dat.all.inc$exposed24)

    #' exposed24 is NA for animals already seropositive in 2023, so exclude
    #' those from the Brier score (BrierScore() requires strictly binary input)
    scoredRows <- !is.na(dat.all.inc$exposed24)
    brierScore <- BrierScore(dat.all.inc$exposed24[scoredRows], pred = dat.all.inc$probPred[scoredRows])

    #' Brier score per fold, so we can check afterwards whether the overall
    #' ranking is consistent across folds or driven by just one or two
    foldBrier <-
      sapply(1:cluster.k, function(k) {
        rows <- which(dat.all.inc$fold == k & !is.na(dat.all.inc$exposed24))
        BrierScore(dat.all.inc$exposed24[rows], pred = dat.all.inc$probPred[rows])
      })
    cv.diagnostics[[fit.name]] <<- list(probPred = dat.all.inc$probPred, foldBrier = foldBrier)

    brierScore
  })

#############################################
##### CV diagnostic checks              #####
#############################################

#' 1. Do the models' out-of-fold predicted probabilities actually differ?
#' If they're near-identical, a lower Brier score for the simpler model isn't
#' meaningful - it just means the extra covariates aren't changing anything.
#' i.e. if their correlation is qutie good then the brier score doesn't mean
#' anything. 
predMatrix2024 <- sapply(cv.diagnostics, function(x) x$probPred)
cat("Correlation between models' out-of-fold predicted probabilities:\n")
print(cor(predMatrix2024, use = "pairwise.complete.obs"))

#' 2. Is the Brier score ranking consistent across all folds, or is it
#' driven by one or two folds?
#' 
foldBrierTab2024 <- sapply(cv.diagnostics, function(x) x$foldBrier)
rownames(foldBrierTab2024) <- paste0("fold", 1:cluster.k)
cat("\nBrier score by fold and model:\n")
print(round(foldBrierTab2024, 3))

cat("\nHow many folds does each model 'win' (lowest Brier)?\n")
foldWinnersInc <- factor(apply(foldBrierTab2024, 1, function(x) names(x)[which.min(x)]),
                      levels = colnames(foldBrierTab2024))
print(table(foldWinnersInc))


#' ## Tabulate and plot results for all models
#' Get WAIC for all models


WAIC.tab.incidence <- cbind(Fixed_effects = unlist(modelSpecsIncidence), IC.inla(mods$outputs[4:length(names(mods$outputs))], crit = "waic"))

covariates <- 
   c(
  #   NDVI = "ndvi", 
  #   NDVI2 = "ndvi.sds2",
  #   NRainDays = "NRainDays.log10", 
  #   PropIrrigated = "proportionIrrigatedCrop",
  #   Elevation = "elevation500m",
  #   BuildDens = "buildingDensity500m.log10", 
  #  PopDens = "popDens500m.log10", 
    Culex = "culexAbundance.sds",
    Aedes = "aedesAbundance.sds" )

for(i in 1:nrow(WAIC.tab.incidence)) {
  for(j in 1:length(covariates)) {
    WAIC.tab.incidence[i, "Fixed_effects"] <-
      gsub(covariates[j], names(covariates)[j], WAIC.tab.incidence[i, "Fixed_effects"])
  }
  WAIC.tab.incidence[i, "Fixed_effects"] <- gsub("\\.sds", "", WAIC.tab.incidence[i, "Fixed_effects"])
}
names(WAIC.tab.incidence) <- gsub("delta", "\u0394", names(WAIC.tab.incidence))

#' Add 5-fold CV results to WAIC table

WAIC.tab.incidence <- cbind(WAIC.tab.incidence, as.matrix(round(cv.results, 4))[rownames(WAIC.tab.incidence), ])
colnames(WAIC.tab.incidence) <- c("Fixed_effects"
                        ,"WAIC"
                        ,"ΔWAIC"
                        ,"Weights"
                        ,"BrierScore")

write.csv(WAIC.tab.incidence,  "modelFitTabIncidence.csv")

Xinc <- as.data.frame(model.matrix(~ prop_year + 
                                  culexAbundance.sds + 
                                    aedesAbundance.sds,
                                data = cattleDatInc))


# extract FOI using the same X object as in the v9 script with aedes and culex and age. 
#' Evaluate the spatial field (object field.est) from our model only at our sampled locations
spatialRE24 <- 
  fm_evaluate(fm_evaluator(mods$inputs$spatial24.null$mesh, loc = cattleDatInc), 
              mods$outputs$spatial24.null$summary.random$field[, point.est])

#' Sum the fixed effects predictions, the village random effect, and the spatial
#' random effect, giving village-level log FOI predictions
fixedPred24 <-
  (as.matrix(Xinc[, 1]) %*% mods$outputs$spatial24.null$summary.fixed[, point.est])
villageRE24 <- 
  mods$outputs$spatial24.null$summary.random$Village[
    match(cattleDatInc$Village, mods$outputs$spatial24.null$summary.random$Village$ID), point.est]
cattleDatInc$logvillageFOI24 <- 
  fixedPred24 +
  villageRE24 +
  spatialRE24

#' Add the herd random effect, giving herd-level log FOI predictions
herdRE24 <-
  mods$outputs$spatial24.null$summary.random$hhID$`0.5quant`[
    match(as.character(cattleDatInc$hhID), mods$outputs$spatial24.null$summary.random$hhID$ID)]
cattleDatInc$logherdFOI24 <- cattleDatInc$logvillageFOI24 + herdRE24

#' Exponentiate log FOI to FOI
cattleDatInc$villageFOI24 <- exp(cattleDatInc$logvillageFOI24)
cattleDatInc$herdFOI24 <- exp(cattleDatInc$logherdFOI24)

#' This is the RAW foi. What this is basically saying, is an animal would experienced 
#' FOI number of of infectious exposure events per year if it kept being exposed indefinitely
#' but this isn't possible because seroconverting is a one time event (in this framewokr)
#' so basically how many exposure events happen on average...

data.frame(
  quantile = c("0%", "2.5%", "50%", "97.5%", "100%"),
  villageFOI24 = quantile(cattleDatInc$villageFOI24, probs = c(0, 0.025, 0.5, 0.975, 1)),
  #villageFOI24lower = quantile(cattleDatInc$villageFOI24lower, probs = c(0, 0.025, 0.5, 0.975, 1)),
  #villageFOI24higher = quantile(cattleDatInc$villageFOI24higher, probs = c(0, 0.025, 0.5, 0.975, 1)),
  row.names = NULL)

#' But what we want to know actually is what is the chance at least one 
#' event happens before the year is out?
#' so this is where this comes in....

#' Annual % of susceptible animals expected to seroconvert (bounded 0-100%,
#' unlike the raw FOI rate above): 1 - exp(-FOI), expressed as a percentage.
#' although it is actually the probability of seroconverting per year. 
#' 1 - exp(-x) is monotonically increasing, so the lower/upper FOI bounds map
#' to the lower/upper bounds of the percentage without swapping.
cattleDatInc$annualPercentSeroconversion <- (1-exp(-cattleDatInc$villageFOI24))*100
#cattleDatInc$annualPercentSeroconversionLower <- (1 - exp(-cattleDatInc$villageFOIlower)) * 100
#cattleDatInc$annualPercentSeroconversionHigher <- (1 - exp(-cattleDatInc$villageFOIhigher)) * 100


data.frame(
  quantile = c("0%", "2.5%", "50%", "97.5%", "100%"),
  annualvillageSC = quantile(cattleDatInc$annualPercentSeroconversion, probs = c(0, 0.025, 0.5, 0.975, 1)),
  #annualvillageScLow = quantile(cattleDatInc$annualPercentSeroconversionLower, probs = c(0, 0.025, 0.5, 0.975, 1)),
  #annualvillageScHigh = quantile(cattleDatInc$annualPercentSeroconversionHigher, probs = c(0, 0.025, 0.5, 0.975, 1)),
  row.names = NULL)

#' Same summary, broken down by District
cattleDatInc %>%
  st_drop_geometry() %>%
  group_by(District) %>%
  summarise(
    n = n(),
    `0%` = quantile(annualPercentSeroconversion, 0),
    `2.5%` = quantile(annualPercentSeroconversion, 0.025),
    `50%` = quantile(annualPercentSeroconversion, 0.5),
    `97.5%` = quantile(annualPercentSeroconversion, 0.975),
    `100%` = quantile(annualPercentSeroconversion, 1)
  )



#' at low FOIs FOI ~= to the annual percent, but they diverge at higher FOIs. 
#' so an expected count (i.e. incident rate) 

#' Express FOI as a standard Poisson incidence rate (cases per 100 animal-years
#' at risk). FOI is already a per-animal-year rate because exposed.yrs was
#' used as the model offset, so this is a re-expression, not a new estimate:
#' N infections in T animal-years of exposure ~ Poisson(FOI * T)
incidenceRate24.tab <-
  data.frame(
    quantile = c("0%", "2.5%", "50%", "97.5%", "100%"),
    incidenceRate100AY = quantile(cattleDatInc$villageFOI24, probs = c(0, 0.025, 0.5, 0.975, 1), na.rm = TRUE) * 100,
    row.names = NULL)

write.csv(incidenceRate24.tab, "incidenceRatePer100AnimalYearsIncidence.csv", row.names = FALSE)

VillagePred$FOI24 <-
  tapply(cattleDatInc$villageFOI24, cattleDatInc$Village, mean)[as.character(VillagePred$Village)]
# VillagePred$obsPos24 <- 
#   tapply(cattleDatInc$exposed24, cattleDatInc$Village, sum, na.rm = TRUE)[as.character(VillagePred$Village)]
# VillagePred$NTested24 <- 
#   tapply(!is.na(cattleDatInc$exposed24), cattleDatInc$Village, sum)[as.character(VillagePred$Village)]

VillagePred$annualPercentSeroconversion2024 <- (1 - exp(-VillagePred$FOI24)) * 100
VillagePred$timesHigher <- 
  round(VillagePred$annualPercentSeroconversion2024/VillagePred$annualPercentSeroconversion)

write.csv(VillagePred, "FOI_analysis/villagePred.csv", row.names = F)

# # make FOI factor
# cutpoints24 <- seq(0, 0.95, 0.05)
# VillagePred$FOI24_level <-
#   cut(VillagePred$FOI24, cutpoints24, 
#       labels = paste0(100 * cutpoints24[-length(cutpoints24)], "-", 100 * cutpoints24[-1], "%"))
# table(VillagePred$FOI24_level)
# 

#hcl_palettes("diverging", n = nlevels(VillagePred$FOI_level), plot = TRUE)
# my_palette24 <- diverging_hcl(nlevels(VillagePred$FOI24_level), palette = "Blue-Red 2")


# incidencePlot <- ggplot() +
#   annotation_map_tile(type = "cartolight", zoomin = 0) +  # You can change the tile type
#   geom_sf(data = VillagePred, 
#           aes(colour = FOI24*100,
#               shape = District),
#           size = 3) +
#   annotation_scale(location = "tl", width_hint = 0.5) +
#   coord_sf(crs = crs.km, expand = T) +
#   scale_colour_viridis_c(breaks = seq(from = 0, to = 95, by = 10)) +
#   guides(colour = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"),
#          shape = guide_legend(title.position = "top")) +
#   theme_minimal() +
#   theme(
#     plot.margin = margin(2, 2, 2, 2),
#     legend.key.size = unit(0.3, "cm"),
#     legend.text = element_text(size = 8),
#     legend.title = element_text(size = 9),
#     legend.spacing.y = unit(0.1, "cm"),
#     legend.position = "bottom",
#     legend.box = "horizontal"
#   ) +
#   labs(colour = "% Exposed \n2023-2024") 
# 
# 
# 
# 
# mean(VillagePred$FOI24)
# var(VillagePred$FOI24)
# range(VillagePred$FOI24)
# sd(VillagePred$FOI24)
# 
# mean(VillagePred$FOI)
# sd(VillagePred$FOI)

##### how many times higher #####

(changePlot <- 
  ggplot() +
  annotation_map_tile(type = "cartolight", zoomin = 0) +  # You can change the tile type
  geom_sf(data = VillagePred,
          aes(fill = timesHigher,  
              shape = District),
          colour = "darkgrey",
          size = 4) +
  scale_fill_viridis_c(option = "D",
                       direction = -1,
                       breaks = scales::breaks_width(2)) +
  scale_shape_manual(values = c(Babati = 21, Hai = 22, Moshi = 24)) +
  #scale_color_distiller(palette = "Set3")+
  annotation_scale(location = "tl", width_hint = 0.5) +
  coord_sf(crs = crs.km) +
  theme_minimal() +
   theme(
     legend.key.size = unit(0.8, "cm"),
     legend.spacing.y = unit(0.3, "cm"),
     plot.margin = margin(2, 2, 2, 2),
     legend.position = "bottom",
     legend.title = element_text(size = 12)
   )+
  labs(fill = "x-fold increase in \nannual seroconversion \n2023-2024")
 # + 
 #  guides(
 #    colour = guide_legend(
 #      keywidth = 1,
 #      keyheight = 1.2
 #    ),
 #    shape = guide_legend(
 #      keywidth = 1,
 #      keyheight = 1.2
 #    )
  )


# villageSeroconversionRate <-
#   cattleDatInc %>%
#   mutate(X = st_coordinates(geometry)[, 1],
#          Y = st_coordinates(geometry)[, 2]) %>%
#   st_drop_geometry() %>%
#   group_by(Village) %>%
#   summarise(X = mean(X), Y = mean(Y), `annual seroconversion` = mean(annualPercentSeroconversion)) %>%
#   st_as_sf(coords = c("X", "Y"), crs = st_crs(cattleDatInc), remove = FALSE)
# 
# changePlot <- 
#   cattleDatInc %>%
#   group_by(Village, District) %>%
#   summarise(`annual seroconversion` = mean(annualPercentSeroconversion)) %>%
#   ggplot() +
#   annotation_map_tile(type = "cartolight", zoomin = 0) +  # You can change the tile type
#   geom_sf(data = VillagePred,
#           aes(colour = timesHigher,
#               shape = District),
#           size = 3) +
#   annotation_scale(location = "tl", width_hint = 0.5) +
#   coord_sf(crs = crs.km, expand = T) +
#   scale_colour_viridis_c(direction = -1) +
#   guides(colour = guide_colourbar(title.position = "top"),
#          shape = guide_legend(title.position = "top")) +
#   theme_minimal() +
#   theme(
#     plot.margin = margin(2, 2, 2, 2),
#     legend.key.size = unit(.8, "cm"),
#     legend.text = element_text(size = 12),
#     legend.title = element_text(size = 12),
#     legend.spacing.y = unit(0.1, "cm"),
#     legend.position = "bottom",
#     legend.box = "horizontal"
#   ) +
#   labs(color = "FOI x-fold increase \n2023-2024")



# plot_grid(changePlot, ggdraw() + draw_plot(vectorIncidence, scale = 0.7),
#           labels = c("A", "B"),
#           rel_widths = c(1.4, 1), 
#           rel_heights = c(8,6), 
#           label_x = 0, 
#           label_y = 0.86
#           )
# ggsave("incidencePanel.pdf")

plot_grid(changePlot, 
          vectorIncidence, 
          nrow = 1, ncol = 2, labels = c("A", "B")
)
ggsave("FOI_analysis/FOIpanelIncidence.pdf")


# ------------------------------------------------------------------------------

#### ignore #####

######################
### w/ log ratio #####
######################

modsVecs$inputs$spatial24Ratio <- mods$inputs$spatial
modsVecs$inputs$spatial24Ratio$obs <-
  bru_obs(
    formula = 
      exposed24 ~ 
      Intercept +
      difference +
      abundanceAedes.sds +
      Village + 
      hhID + 
      field + 
      offset(log(exposed.yrs)),
    family = "binomial",
    control.family = list(link = "cloglog"),
    data = cattleDatInc)

#' Fit model
modsVecs$outputs$spatial24Ratio <- 
  bru(
    components = modsVecs$inputs$spatial$components,
    obs = modsVecs$inputs$spatial24Ratio$obs,
    options = list(
      control.fixed = mods$priors$prior.fixed,
      control.predictor = list(compute = TRUE),
      control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, 
                             config = TRUE,
                             return.marginals.predictor = TRUE)))
#summary(model2.inla)
summary(modsVecs$outputs$spatial24Ratio)
exp(modsVecs$outputs$spatial24Ratio$summary.fixed)

# visualise this 
fit.summary.fixed <- modsVecs$outputs$spatial24Ratio$summary.fixed
fit.summary.fixed$ID <- rownames(fit.summary.fixed)
fit.summary.fixed$ID <- factor(fit.summary.fixed$ID, fit.summary.fixed$ID)

# Plot
ratioIncidence <- 
  fit.summary.fixed %>%
  # mutate(ID = fct_recode(ID
  #                        ,Culex = "culexAbundance.sds"
  #                        ,Aedes = "abundanceAedes.sds"
  # )) %>%
  ggplot(aes(y = ID, x = exp(mean))) +
  geom_pointrange(aes(xmin = exp(`0.025quant`), xmax = exp(`0.975quant`))) +
  geom_vline(xintercept = 1, linetype = 2) +
  # scale_y_discrete(labels = c(
  #   "Intercept"
  #   ,expression(italic("Culex"))
  #   ,expression(italic("Aedes"))
  # )) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")+
  labs(y = "", x = "Estimate \u00B1 95% CI")

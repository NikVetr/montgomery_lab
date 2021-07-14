#### describe animation ####

# open on frame w/ no axes, just title "Inference of Wiggly Functions from Data"
# bayes can pop in (using speech bubbles and xkcd font popping in one at a time) all "Hi folks! Today we'll be exploring how soon into a non-linear temporal process it's end behavior can be well understood. We'll be using a linear combination of logistic functions as our non-linear model of choice, featuring both growth and decay components (draw these) whose sum determine (sum them) the system's overall behavior. 
# we simulate sampling times in that most natural of ways, by drawing from a rescaled dirichlet(5,5,5...) to yield an expected rate at each time, which we use to parameterize a negative-binomial distribution to realize rates by which poisson-distributed counts may be drawn (do the lightning bolt thing). We'll fit the true, perfectly specified model to these simulated data one observation at a time to show you the true power of my ABS: Applied Bayesian Sta-"
# then Fisher rises from the bottom of the frame saying "Not so fast, BARNIE! Your paltry method of inverse probability will not avail you here! Truly, it is only the Method of Maximum Likelihood that provides a principled means of estimating model parameters."
# then he'll pause and say "Wait, how do I express the logistic MLE in closed form again? Little help here boys!" and newton and efron and turing will appear saying "You have my compute!" "And my method!" "And MY bootstraps!"
# have them comment throughout, e.g. "incredible" pun on credible intervals, "so uninformed" on flat priors, "A MAP by any other name" for penalized likelihood, etc.


# phase 1a -- fullscreen plot, visualize curve growing, visualize exponentially distributed sampling events 
# visualize gamma distribution moving, visualize poisson distribution forming, shoot arrows from poisson distribution to points that grow and then shrink
# like the lightning bolts of the mighty zeus
# phase 1b -- shrink plot to the upper left corner, then make three copies into the other three corners, specifying each copy's nature before shrinking it (and having bodybuilders poke out)
# phase 2a -- fade in the different model specifications, as well as the prior distributions for each model's parameters in the appropriate subplots, zooming out as needed
# phase 2b -- pan a line across the plot, w/ either Bayes' head (+ the Stan logo) or Newton's head (+ the R logo) above (maybe give newton googly eyes?)
# ORRRR put arnold bayes and... 
# shifting prior of model parameters into posterior, and plot 89% HPDIs and 89% posterior predictive intervals in shaded region
# shade area to the left of the line with a light grey, and color open points according to the PSIS / LPPD score
# in total, have 4 models: 1) true model, flat priors 2) true model, shrinkage priors 3) cubic splines, shrinkage priors 4) true model, bfgs MLE w/ bootstrap proportions
# fade in gold star on winner, fade everything out

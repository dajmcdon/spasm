#'@export
particleFilter <- function(Yin, stateTransition, 
               f = dnorm,            #density of Yt given Xt
               nParticles = 1000, generateParticles = rnorm, verbose = FALSE){
    Y = as.matrix(Yin)
    particles = as.matrix(generateParticles(nParticles))
    output = array(dim = c(nrow(Y), ncol(particles), nParticles))
    means = array(dim = c(nrow(Y), ncol(particles)))
    for(t in 1:nrow(Y)){
        if(verbose){
            print(t)
        }
        weights = f(Y[t,], particles)
        indices = sample(1:nrow(particles), prob = weights, replace = TRUE)
        output[t,,] = particles[indices,]
        means[t,] = colMeans(particles[indices,])
        particles = as.matrix(stateTransition(as.matrix(particles[indices,])))
    }
    list(particles = output, means = means)
}

ilogit <- function(x){
    1/(1 + exp(-1*x))
}

stateTransition <- function(x){
    0.9*x + rnorm(length(x), mean = 0, sd = sqrt(0.1))
}

observationFunction <- function(x){
    4*ilogit(x)*(1 - ilogit(x)) + rnorm(length(x), mean = 0, sd = sqrt(0.1))
}

generateDat <- function(TT){
    x = double(TT)
    y = double(TT)
    x[1] = rnorm(1)
    for(t in 2:TT){
        x[t] = stateTransition(x[t - 1])
    }
    y = observationFunction(x)
    list(x = x, y = y)
}

yGivenX <- function(y, x){
    dnorm(y, mean = 4*ilogit(x)*(1-ilogit(x)), sd = sqrt(0.1))
}

gp <- function(n){
    matrix(rnorm(n), nrow = n, ncol = 1)
}

#' @export
filterTest <- function(TT, np){
    dat =  generateDat(200)
    filtered = particleFilter(dat$y, stateTransition = stateTransition, f = yGivenX, generateParticles = gp, nParticles = 200000)
    plot(dat$x ~ c(1:200))
    lines(filtered$means)
    remove(dat)
    remove(filtered)
    
}
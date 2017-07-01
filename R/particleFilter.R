#'@export
particleFilter <- function(Yin, stateTransition, 
               f = dnorm,            #density of Yt given Xt
               nParticles = 1000, generateParticles = rnorm, verbose = FALSE){
    Y = as.matrix(Yin)
    particles = as.matrix(generateParticles(nParticles))
    output = array(dim = c(nrow(Y), ncol(particles), nParticles))
    for(t in 1:nrow(Y)){
        if(verbose){
            print(t)
        }
        weights = f(Y[t,], particles)
        weights = weights/sum(weights)
        indices = sample(1:nrow(particles), prob = weights, replace = TRUE)
        output[t,,] = particles[indices,]
        particles = as.matrix(stateTransition(as.matrix(particles[indices,])))
    }
    output
}

dat =  generateSPASM(1000, 3, 2)
stateTransition <- function(x){
    (0.9 * x) + (0.1 * rnorm(length(x)))
}

stateT <- function(X){
    t(apply(X, 1, stateTransition))
}

generateParticles <- function(n){
    matrix(rep(0, n*3), nrow = n, ncol = 3)
}

yGivenX <- function(y, x){
    dmvnorm(y, mean = as.vector(dat$Zt %*% sin((1:length(x))/(2*pi) * x)), sigma = diag(0.1, length(y)))
}

ygx <- function(y, x){
    apply(x, MARGIN = 1, function(z){yGivenX(y, z)})
}

gp <- function(n){
    matrix(0, nrow = n, ncol = 3)
}


#' @export
filterTest <- function(){
    particleFilter(dat$y, stateTransition = stateT, f = ygx, generateParticles = gp, verbose = TRUE)
}
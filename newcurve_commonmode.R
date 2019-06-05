#### Clear workspace ####
rm(list=ls())

#### Include libraries ####
library(ggplot2)

#### Declaration of functions ####
# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  sqrt(mean(error^2))
}

# Function that returns Mean Absolute Error
mae <- function(error)
{
  mean(abs(error))
}

#### Initialize parameters ####

benchmarks <- c("Gold-rader","Blowfish","SHA")
frequency_MHz <- c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500)
voltage <- c(0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98,1.09,1.23,1.23,1.23,1.23)


#### Fitting execution time over frequency ####
# Read data from the file and assign coloumn names
benchmark_time <- read.table("data_time_kamesh.txt")
colnames(benchmark_time) <- c("Gold-rader","Blowfish","SHA")

my_benchmark_time <- benchmark_time
prediction_time <- benchmark_time

# Build a matrix for the constants and initilize them to '0'
# and assign row names and coloumn names to the matrix.
constant_parameters <- matrix(0,length(benchmarks),3)
colnames(constant_parameters) <- c("ccb","cck","beta")
rownames(constant_parameters) <- c("Gold-rader","Blowfish","SHA")

starting_values <- rbind(c(3.0,0.072,1.27),
                         c(10.0,0.072,1.27),
                         c(10.0,0.072,1.27))

x <- frequency_MHz
y <- benchmark_time

# Run the regression analysis (fitting) for each benchmark
for(size in 1:length(benchmarks)) {
  cat("Benchmark:",benchmarks[size],"\n")
  
  x_fit <- frequency_MHz/1e3
  # regression analysis
  regression <- nls(y[,size] ~ (constant1 * ((1 / (x_fit - constant2))  + constant3)),
                    control=list(maxiter = 1000),
                    start = list(constant1=starting_values[size,1],
                                 constant2=starting_values[size,2],
                                 constant3=starting_values[size,3]),
                    trace = T,
                    algorithm = "port",
                    lower = c(0,0,-100))
  
  constant_parameters[size,"ccb"] <- coef(regression)[1]
  constant_parameters[size,"cck"] <- coef(regression)[2]
  constant_parameters[size,"beta"] <- coef(regression)[3]
  
  prediction_time[,size] <- (constant_parameters[size,"ccb"] / (x_fit-constant_parameters[size,"cck"]) + constant_parameters[size,"beta"])
  
}

  #### Plot the execution time vs frequency ####
plot(x=1,
     y=1,
     main=paste("Execution time vs Frequency"),
     xaxp  = c(0, 1500, 15),
     yaxp  = c(0, 800, 8),
     xlab= "Frequency (MHz)",
     ylab= "Execution time (s)",
     type="n",
     ylim=range(my_benchmark_time[,1:3],prediction_time[,1:3]),
     xlim=range(frequency_MHz))

abline(v=seq(0,max(x),100), h=seq(0,max(y),100), col="gray", lty="dotted")

legend("topright",
       #title=expression(paste("Benchmarks")),
       paste("",benchmarks,sep=""),
       col=(1:8)+1,
       lty=1,
       pch=1:8,
       ncol=1,
       bg="white",
       cex=0.75)

for (size in 1:length(benchmarks)) {
  lines(x=frequency_MHz,
        y=my_benchmark_time[,size],
        col="black",
        type="o",
        lty="dashed",
        pch=size)
  
  lines(x=frequency_MHz,
        y=prediction_time[,size],
        col=size+1,
        type="o",
        lty="solid",
        pch=size)
}

#### Print summary and constant parameter from the regression ####
print(summary(regression))
print(constant_parameters)

#### Fitting power over frequency ####

my_power_estimates <- read.table("power_frequency_kamesh.txt",sep="")
colnames(my_power_estimates) <- c("Gold-rader","Blowfish","SHA")

my_benchmark_power <- my_power_estimates
prediction_power <- my_power_estimates

par(mfrow=c(1,1))

#Matrix to hold the fitted parameters
my_parameters_power <- matrix(0,length(benchmarks),3)
colnames(my_parameters_power) <- c("Pstatic", "gamma", "naC")
rownames(my_parameters_power) <- c("Gold-rader","Blowfish","SHA")
  
for (size in 1:length(benchmarks)) {
  
  cat("Benchmark:",benchmarks[size],"\n")
  
  f <- frequency_MHz/1e3
  V <- voltage
  P <- my_benchmark_power[,size]
  
  cat("\tplotting power vs frequency curves...\n")
  
  plot(x=frequency_MHz,
       y=P,
       col="red",
       main=paste("Average power vs Frequency:", benchmarks[size]),
       xlab="Frequency (MHz)",
       ylab="Average power (mW)",
       type="n",
       lty="dashed",
       xlim=range(frequency_MHz),
       ylim=c(0,max(my_benchmark_power[,size]+5))
  )
  
  abline(v=(seq(0,1600,100)), h=(seq(0,4500,500)), col="lightgray", lty="dotted")
  
  # FITTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  cat("\tfitting default ...\n")
  
  regression_power <- nls(P ~ Pstatic + (1+V*gamma)*naC*f^2*V^2,
                          control=list(maxiter = 1000,
                                       #minFactor=1e-5,
                                       warnOnly=T),
                          start = list(Pstatic = 1.5,
                                       gamma = 1.2,
                                       naC   = 10),
                          trace = T,
                          algorithm = "port",
                          lower = c(-100,-100,-100),
                          upper = c(500,100,900))
  
  my_parameters_power[size, "Pstatic"] <- Pstatic <- coef(regression_power)[1]
  my_parameters_power[size, "gamma"]   <- gamma   <- coef(regression_power)[2]
  my_parameters_power[size, "naC"]     <- naC     <- coef(regression_power)[3]
  
  Vfit <- voltage
  P <- my_benchmark_power[,size]
  
  lines(x=frequency_MHz,
        y=P,
        col="blue",
        type="l")
  
  lines(x=frequency_MHz,
        y=(Pstatic+(1+Vfit*gamma)*naC*f^2*Vfit^2),
        col="red",
        type="o",
        pch=3,
        lty="dashed")
  
  prediction_power[,size] <- (Pstatic+(1+Vfit*gamma)*naC*f^2*Vfit^2)
  
}

#### Print summary and constant parameter from the regression ####
print(summary(regression_power))
print(my_parameters_power)

#### Energy modeling and fitting ####

experimental_energy <- my_benchmark_power * my_benchmark_time

theoretical_energy <- prediction_power * prediction_time

#### Plotting of Energy vs Frequency curves ####
plot(x=1,
     y=1,
     main=paste("Total energy consumption vs Frequency"),
     xaxp  = c(0, 1600, 16),
     yaxp  = c(0, 400000, 10),
     xlab= "Frequency (MHz)",
     ylab= "Total energy consumption (mJ)",
     type="n",
     ylim=range(experimental_energy[,1:3],theoretical_energy[,1:3]),
     xlim=range(frequency_MHz))

abline(v=(seq(0,1600,100)), h=(seq(0,400000,50000)), col="lightgray", lty="dotted")

for (size in 1:length(benchmarks)) {
  lines(x=frequency_MHz,
        y=theoretical_energy[,size],
        col="black",
        type="o",
        lty="dashed",
        pch=size)
  
  lines(x=frequency_MHz,
        y=experimental_energy[,size],
        col=size+1,
        type="o",
        lty="solid",
        pch=size)
}


df2 <- data.frame(type=rep(c("experimental", "theoretical"), each=15),
                  frequency_sha=rep(c(frequency_MHz),2),
                  benchmark_sha=c(experimental_energy$SHA,theoretical_energy$SHA))

meow <- ggplot(data=df2, aes(x=frequency_sha, y=benchmark_sha, group=type)) +
  geom_line(aes(linetype=type, color=type))+
  geom_point(aes(shape=type, color=type))

print(meow)

#### Printing results ####

error_goldrader <- experimental_energy$`Gold-rader` - theoretical_energy$`Gold-rader`
error_Blowfish <- experimental_energy$`Blowfish` - theoretical_energy$`Blowfish`
error_sha <- experimental_energy$`SHA` - theoretical_energy$`SHA`

cat("\t \n")
cat("\t Goodness-of-Fit Statistics \n")
cat("\t Printing root mean square error (rmse) and mean absolute error (mae) \n")
cat("\t ==================================================================== \n")
cat("\t \n")
cat("\t Gold-rader: \n")
print(rmse(error_goldrader))
print(mae(error_goldrader))
cat("\t \n")
cat("\t Blowfish: \n")
print(rmse(error_Blowfish))
print(mae(error_Blowfish))
cat("\t \n")
cat("\t SHA: \n")
print(rmse(error_sha))
print(mae(error_sha))

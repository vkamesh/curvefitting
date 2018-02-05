fdat <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
pdat <- c(2.02, 3.19, 4.89, 6.43, 8.84, 11.65, 15.04, 18.69, 22.81, 27.85, 31.09, 38.13, 45.74, 49.87, 54.15)
pdata <- data.frame(y=fdat, t=pdat)

selfstart1 <- nls(y ~ SSlogis(t, pstatic, gamma, naC), data=pdata)
print(summary(selfstart1))

fdat <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
pdat <- c(2.52, 3.52, 4.6, 5.75, 7.2, 8.81, 10.34, 12.35, 14, 15.67, 17.66, 20, 23.09, 26.54, 29.58)
pdata <- data.frame(y=fdat, t=pdat)

selfstart2 <- nls(y ~ SSlogis(t, pstatic, gamma, naC), data=pdata)
print(summary(selfstart2))

fdat <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
pdat <- c(2.8, 3.53, 4.48, 5.57, 6.76, 8.11, 9.49, 11.25, 12.9, 14.52, 16.35, 19.17, 21.7, 23.67, 25.75)
pdata <- data.frame(y=fdat, t=pdat)

selfstart3 <- nls(y ~ SSlogis(t, pstatic, gamma, naC), data=pdata)
print(summary(selfstart3))
temp = read.csv("/home/morikie/Documents/Forschungspraktikum/C++/polar/scripts/notFoundPAS.csv", header=FALSE)

x = seq(0, 491)
y = rep(0, 492)

for (i in 1:length(y)) {
  y[i] = mean(temp[, i])
}

plot(x, y, main="notFoundPAS data set", type="h")

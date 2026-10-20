


library("yaml")

n <- "Normal"

str(x <- read_yaml(sprintf("config-%s.yaml", n)))
names(x$parameters)

names(x$params)



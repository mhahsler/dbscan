# Source: http://glaros.dtc.umn.edu/gkhome/cluto/cluto/download

chameleon_ds4 <- read.table("t4.8k.dat")
chameleon_ds5 <- read.table("t5.8k.dat")
chameleon_ds7 <- read.table("t7.10k.dat")
chameleon_ds8 <- read.table("t8.8k.dat")

colnames(chameleon_ds4) <- colnames(chameleon_ds5) <- colnames(chameleon_ds7) <- colnames(chameleon_ds8) <- c("x", "y")

plot(chameleon_ds4)
plot(chameleon_ds5)
plot(chameleon_ds7)
plot(chameleon_ds8)

save(chameleon_ds4, chameleon_ds5, chameleon_ds7, chameleon_ds8, 
     file="Chameleon.rda")

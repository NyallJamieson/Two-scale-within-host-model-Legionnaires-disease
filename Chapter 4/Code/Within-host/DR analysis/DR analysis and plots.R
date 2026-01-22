# Read in results
Erlang_averaged <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Simulation results/Erlang/Results.Erlang.Averaged.csv")[,2:3]
Erlang_intermediate <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Simulation results/Erlang/Erlang.Intermediate.2(64).csv")[,2:3]
Erlang_pois <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Simulation results/Erlang/Results.Erlang.Pois.csv")[,2:3]

# Read in results
Burr_averaged <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Simulation results/Burr/Results.Burr.Averaged.csv")[,2:3]
Burr_intermediate <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Simulation results/Burr/Burr.Intermediate.csv")[,2:3]
Burr_pois <- read.csv("~/PhD-work/Chapter 4/Data/Within-host/Simulation results/Burr/Results.Burr.Pois.csv")[,2:3]

# Now go through models and find ID50 and uncertainty
dose <- seq(1,500,by=1)

prob.1 <- unname(table(Erlang_averaged[,1]))/1000
prob.2 <- unname(table(Erlang_intermediate[,1]))/1000
prob.3 <- unname(table(Erlang_pois[,1]))/1000
prob.4 <- unname(table(Burr_averaged[,1]))/1000
prob.5 <- unname(table(Burr_intermediate[,1]))/1000
prob.6 <- unname(table(Burr_pois[,1]))/1000

# Fit exponential DR models now
dr1 <- minpack.lm::nlsLM(prob.1~1-exp(-log(2)*dose/d50),start=list(d50=1),control = nls.control(maxiter = 1000))
dr2 <- minpack.lm::nlsLM(prob.2~1-exp(-log(2)*dose/d50),start=list(d50=1),control = nls.control(maxiter = 1000))
dr3 <- minpack.lm::nlsLM(prob.3~1-exp(-log(2)*dose/d50),start=list(d50=1),control = nls.control(maxiter = 1000))
dr4 <- minpack.lm::nlsLM(prob.4~1-exp(-log(2)*dose/d50),start=list(d50=1),control = nls.control(maxiter = 1000))
dr5 <- minpack.lm::nlsLM(prob.5~1-exp(-log(2)*dose/d50),start=list(d50=1),control = nls.control(maxiter = 1000))
dr6 <- minpack.lm::nlsLM(prob.6~1-exp(-log(2)*dose/d50),start=list(d50=1),control = nls.control(maxiter = 1000))

# Estimate ID50 and standard error
A1 <- summary(dr1)$coefficients["d50", c("Estimate","Std. Error")]
A2 <- summary(dr2)$coefficients["d50", c("Estimate","Std. Error")]
A3 <- summary(dr3)$coefficients["d50", c("Estimate","Std. Error")]
A4 <- summary(dr4)$coefficients["d50", c("Estimate","Std. Error")]
A5 <- summary(dr5)$coefficients["d50", c("Estimate","Std. Error")]
A6 <- summary(dr6)$coefficients["d50", c("Estimate","Std. Error")]

# Now want to compare to P2 work
source("~/PhD-work/Chapter 2/Code/DR and TDR analysis/DR plot for paper.R")
R1 <- append(0, prob.1)
R2 <- append(0, prob.2)
R3 <- append(0, prob.3)
R4 <- append(0, prob.4)
R5 <- append(0, prob.5)
R6 <- append(0, prob.6)

# Do plot similar to P2

# Save plot with 600 DPI resolution
png("DR1.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(dose_M,response_M,main="Dose-response for Legionnaires' disease",xlab="Deposited dose (number of Legionella)",ylab="Probability of symptom onset",lwd=2,ylim=c(0,1),xlim=c(0,100))
lines(seq(0, 500, by = 1), R1, col=2, lwd=2)
lines(seq(0, 500, by = 1), R2, col=3, lwd=2)
lines(seq(0, 500, by = 1), R3, col=4, lwd=2)
lines(seq(0, 500, by = 1), R4, col=5, lwd=2)
lines(seq(0, 500, by = 1), R5, col=6, lwd=2)
lines(seq(0, 500, by = 1), R6, col=7, lwd=2)
lines(doses,1-exp(-beta_hat*doses),lwd=2,col="black")
lines(doses,1-exp(-ci_upper*doses),col="black",lty=2)
lines(doses,1-exp(-ci_lower*doses),col="black",lty=2)
for (i in 1:4){
  arrows(x0=dose_M[i], y0=point_lower[i], x1=dose_M[i], y1=point_upper[i], length=0.05, angle=90, code=3)
}
legend("bottomright",legend=c("Erlang intermediate", "Erlang Poisson", "Erlang constant", "Burr intermediate", "Burr Poisson", "Burr constant","Binomial likelihood \n to experimental results"), col = c(2, 3, 4, 5, 6, 7, 1), lty = c(1, 1, 1, 1, 1, 1, 1), lwd = c(2, 2, 1, 2, 2, 2), bty = "n")

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()


cols <- RColorBrewer::brewer.pal(9, "Set1")
cols

# Save plot with 600 DPI resolution
png("DR2.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(seq(8.75, 9.05, by = 0.001), dnorm(seq(8.75, 9.05, by = 0.001), mean = A1[1], sd = A1[2]), type = "l", lwd = 2, col = cols[1], ylim = c(0, 33), xlab = expression(paste(ID[50], " (number of Legionella)")), ylab = "Probability density")
lines(seq(8.75, 9.05, by = 0.001), dnorm(seq(8.75, 9.05, by = 0.001), mean = A2[1], sd = A2[2]), lwd = 2, col = cols[2])
lines(seq(8.75, 9.05, by = 0.001), dnorm(seq(8.75, 9.05, by = 0.001), mean = A3[1], sd = A3[2]), lwd = 2, col = cols[3])
lines(seq(8.75, 9.05, by = 0.001), dnorm(seq(8.75, 9.05, by = 0.001), mean = A4[1], sd = A4[2]), lwd = 2, col = cols[4])
lines(seq(8.75, 9.05, by = 0.001), dnorm(seq(8.75, 9.05, by = 0.001), mean = A5[1], sd = A5[2]), lwd = 2, col = cols[5])
lines(seq(8.75, 9.05, by = 0.001), dnorm(seq(8.75, 9.05, by = 0.001), mean = A6[1], sd = A6[2]), lwd = 2, col = cols[6])
lines(seq(8.75, 9.05, by = 0.001), dnorm(seq(8.75, 9.05, by = 0.001), mean = 8.83699058, sd = 0.01779069), lwd = 2, col = cols[7], lty = 2)
lines(seq(8.75, 9.05, by = 0.001), dnorm(seq(8.75, 9.05, by = 0.001), mean = 8.94816397, sd = 0.01544307), lwd = 2, col = cols[8], lty = 2)
lines(seq(8.75, 9.05, by = 0.001), dnorm(seq(8.75, 9.05, by = 0.001), mean = 8.9706873, sd = 0.0200399), lwd = 2, col = cols[9], lty = 2)
legend("top", legend = c("Erlang averaged", "Erlang Poisson", "Erlang intermediate", "Burr averaged", "Burr Poisson", "Burr intermediate", "[7] model A", "[7] model B", "[7] model C"), col = c(cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8], cols[9]), lwd = 2, lty = c(1, 1, 1, 1, 1, 1, 2, 2, 2), ncol = 3, horiz = FALSE)

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()


# Now TDR analysis
# Get TDR data in correct form
d <- seq(1,500,by=1)
t <- seq(0,200,by=0.5)

val <- c()
for (i in t){
  print(i)
  for (j in d){
    val <- append(val,length(Erlang_averaged[which(Erlang_averaged[,1]==j&Erlang_averaged[,2]<=i),2])/1000)}}
data1 <- expand.grid(X=d,Y=t)
data1$Z <- val
data1[length(d)*length(t)+1,] <- c(1000,1000,1)

val <- c()
for (i in t){
  print(i)
  for (j in d){
    val <- append(val,length(Erlang_pois[which(Erlang_pois[,1]==j&Erlang_pois[,2]<=i),2])/1000)}}
data2 <- expand.grid(X=d,Y=t)
data2$Z <- val
data2[length(d)*length(t)+1,] <- c(1000,1000,1)

val <- c()
for (i in t){
  print(i)
  for (j in d){
    val <- append(val,length(Erlang_intermediate[which(Erlang_intermediate[,1]==j&Erlang_intermediate[,2]<=i),2])/1000)}}
data3 <- expand.grid(X=d,Y=t)
data3$Z <- val
data3[length(d)*length(t)+1,] <- c(1000,1000,1)

val <- c()
for (i in t){
  print(i)
  for (j in d){
    val <- append(val,length(Burr_averaged[which(Burr_averaged[,1]==j&Burr_averaged[,2]<=i),2])/1000)}}
data4 <- expand.grid(X=d,Y=t)
data4$Z <- val
data4[length(d)*length(t)+1,] <- c(1000,1000,1)

val <- c()
for (i in t){
  print(i)
  for (j in d){
    val <- append(val,length(Burr_pois[which(Burr_pois[,1]==j&Burr_pois[,2]<=i),2])/1000)}}
data5 <- expand.grid(X=d,Y=t)
data5$Z <- val
data5[length(d)*length(t)+1,] <- c(1000,1000,1)

val <- c()
for (i in t){
  print(i)
  for (j in d){
    val <- append(val,length(Burr_intermediate[which(Burr_intermediate[,1]==j&Burr_intermediate[,2]<=i),2])/1000)}}
data6 <- expand.grid(X=d,Y=t)
data6$Z <- val
data6[length(d)*length(t)+1,] <- c(1000,1000,1)


# Save plot with 600 DPI resolution
tiff("Fig1.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# TDR heatmap
ggplot(data1, aes(X, Y, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colours = rainbow(5)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(50, 200), expand = c(0, 0)) +
  labs(x = "Dose", y = "Time in Hours", fill = "Probability") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 20),      # Axis labels
    axis.text = element_text(size = 20),       # Axis tick labels
    legend.title = element_text(size = 20, margin = margin(b = 10)),    # Legend title
    legend.text = element_text(size = 20),      # Legend labels
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5) 
  )

# Close device
dev.off()

# Save plot with 600 DPI resolution
tiff("Fig2.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# TDR heatmap
ggplot(data2, aes(X, Y, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colours = rainbow(5)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(50, 200), expand = c(0, 0)) +
  labs(x = "Dose", y = "Time in Hours", fill = "Probability") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 20),      # Axis labels
    axis.text = element_text(size = 20),       # Axis tick labels
    legend.title = element_text(size = 20, margin = margin(b = 10)),    # Legend title
    legend.text = element_text(size = 20),      # Legend labels
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5) 
  )

# Close device
dev.off()

# Save plot with 600 DPI resolution
tiff("Fig3.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# TDR heatmap
ggplot(data3, aes(X, Y, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colours = rainbow(5)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(50, 200), expand = c(0, 0)) +
  labs(x = "Dose", y = "Time in Hours", fill = "Probability") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 20),      # Axis labels
    axis.text = element_text(size = 20),       # Axis tick labels
    legend.title = element_text(size = 20, margin = margin(b = 10)),    # Legend title
    legend.text = element_text(size = 20),      # Legend labels
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5) 
  )

# Close device
dev.off()

# Save plot with 600 DPI resolution
tiff("Fig4.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# TDR heatmap
ggplot(data4, aes(X, Y, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colours = rainbow(5)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(50, 200), expand = c(0, 0)) +
  labs(x = "Dose", y = "Time in Hours", fill = "Probability") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 20),      # Axis labels
    axis.text = element_text(size = 20),       # Axis tick labels
    legend.title = element_text(size = 20, margin = margin(b = 10)),    # Legend title
    legend.text = element_text(size = 20),      # Legend labels
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5) 
  )

# Close device
dev.off()

# Save plot with 600 DPI resolution
tiff("Fig5.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# TDR heatmap
ggplot(data5, aes(X, Y, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colours = rainbow(5)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(50, 200), expand = c(0, 0)) +
  labs(x = "Dose", y = "Time in Hours", fill = "Probability") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 20),      # Axis labels
    axis.text = element_text(size = 20),       # Axis tick labels
    legend.title = element_text(size = 20, margin = margin(b = 10)),    # Legend title
    legend.text = element_text(size = 20),      # Legend labels
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5) 
  )

# Close device
dev.off()

# Save plot with 600 DPI resolution
tiff("Fig6.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# TDR heatmap
ggplot(data6, aes(X, Y, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colours = rainbow(5)) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(50, 200), expand = c(0, 0)) +
  labs(x = "Dose", y = "Time in Hours", fill = "Probability") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 20),      # Axis labels
    axis.text = element_text(size = 20),       # Axis tick labels
    legend.title = element_text(size = 20, margin = margin(b = 10)),    # Legend title
    legend.text = element_text(size = 20),      # Legend labels
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5) 
  )

# Close device
dev.off()



# Get moments of incubation periods
m1_1 <- c()
m1_2 <- c()
m1_3 <- c()
m1_4 <- c()
m2_1 <- c()
m2_2 <- c()
m2_3 <- c()
m2_4 <- c()
m3_1 <- c()
m3_2 <- c()
m3_3 <- c()
m3_4 <- c()
m4_1 <- c()
m4_2 <- c()
m4_3 <- c()
m4_4 <- c()
m5_1 <- c()
m5_2 <- c()
m5_3 <- c()
m5_4 <- c()
m6_1 <- c()
m6_2 <- c()
m6_3 <- c()
m6_4 <- c()

m1_5 <- c()
m2_5 <- c()
m3_5 <- c()
m4_5 <- c()
m5_5 <- c()
m6_5 <- c()

m1_6 <- c()
m2_6 <- c()
m3_6 <- c()
m4_6 <- c()
m5_6 <- c()
m6_6 <- c()

for (i in 1:500){
  times1 <- Erlang_averaged[which(Erlang_averaged[,1] == i), 2]
  m1_1 <- append(m1_1,mean(times1))
  m1_2 <- append(m1_2,var(times1))
  m1_3 <- append(m1_3,PerformanceAnalytics::skewness(times1))
  m1_4 <- append(m1_4,PerformanceAnalytics::kurtosis(times1,method="moment"))
  times2 <- Erlang_pois[which(Erlang_pois[,1] == i), 2]
  m2_1 <- append(m2_1,mean(times2))
  m2_2 <- append(m2_2,var(times2))
  m2_3 <- append(m2_3,PerformanceAnalytics::skewness(times2))
  m2_4 <- append(m2_4,PerformanceAnalytics::kurtosis(times2,method="moment"))
  times3 <- Erlang_intermediate[which(Erlang_intermediate[,1] == i), 2]
  m3_1 <- append(m3_1,mean(times3))
  m3_2 <- append(m3_2,var(times3))
  m3_3 <- append(m3_3,PerformanceAnalytics::skewness(times3))
  m3_4 <- append(m3_4,PerformanceAnalytics::kurtosis(times3,method="moment"))
  times4 <- Burr_averaged[which(Burr_averaged[,1] == i), 2]
  m4_1 <- append(m4_1,mean(times4))
  m4_2 <- append(m4_2,var(times4))
  m4_3 <- append(m4_3,PerformanceAnalytics::skewness(times4))
  m4_4 <- append(m4_4,PerformanceAnalytics::kurtosis(times4,method="moment"))
  times5 <- Burr_pois[which(Burr_pois[,1] == i), 2]
  m5_1 <- append(m5_1,mean(times5))
  m5_2 <- append(m5_2,var(times5))
  m5_3 <- append(m5_3,PerformanceAnalytics::skewness(times5))
  m5_4 <- append(m5_4,PerformanceAnalytics::kurtosis(times5,method="moment"))
  times6 <- Burr_intermediate[which(Burr_intermediate[,1] == i), 2]
  m6_1 <- append(m6_1,mean(times6))
  m6_2 <- append(m6_2,var(times6))
  m6_3 <- append(m6_3,PerformanceAnalytics::skewness(times6))
  m6_4 <- append(m6_4,PerformanceAnalytics::kurtosis(times6,method="moment"))
  m1_5 <- append(m1_5, median(times1))
  m2_5 <- append(m2_5, median(times2))
  m3_5 <- append(m3_5, median(times3))
  m4_5 <- append(m4_5, median(times4))
  m5_5 <- append(m5_5, median(times5))
  m6_5 <- append(m6_5, median(times6))
  m1_6 <- append(m1_6, IQR(times1))
  m2_6 <- append(m2_6, IQR(times2))
  m3_6 <- append(m3_6, IQR(times3))
  m4_6 <- append(m4_6, IQR(times4))
  m5_6 <- append(m5_6, IQR(times5))
  m6_6 <- append(m6_6, IQR(times6))
}

R1 <- read.csv("~/PhD-work/Chapter 2/Data/Simulated CTMC results/Results.Constant.csv")[,2:3]
R2 <- read.csv("~/PhD-work/Chapter 2/Data/Simulated CTMC results/Results.NegBin.csv")[,2:3]
R3 <- read.csv("~/PhD-work/Chapter 2/Data/Simulated CTMC results/Results.Hetero.Exact.1to500.csv")[,2:3]

M1_1 <- c()
M1_2 <- c()
M1_3 <- c()
M1_4 <- c()
M1_5 <- c()
M2_1 <- c()
M2_2 <- c()
M2_3 <- c()
M2_4 <- c()
M2_5 <- c()
M3_1 <- c()
M3_2 <- c()
M3_3 <- c()
M3_4 <- c()
M3_5 <- c()

M1_6 <- c()
M2_6 <- c()
M3_6 <- c()

for (i in 1:500){
  times1 <- R1[which(R1[,1] == i), 2]
  M1_1 <- append(M1_1,mean(times1))
  M1_2 <- append(M1_2,var(times1))
  M1_3 <- append(M1_3,PerformanceAnalytics::skewness(times1))
  M1_4 <- append(M1_4,PerformanceAnalytics::kurtosis(times1,method="moment"))
  times2 <- R2[which(R2[,1] == i), 2]
  M2_1 <- append(M2_1,mean(times2))
  M2_2 <- append(M2_2,var(times2))
  M2_3 <- append(M2_3,PerformanceAnalytics::skewness(times2))
  M2_4 <- append(M2_4,PerformanceAnalytics::kurtosis(times2,method="moment"))
  times3 <- R3[which(R3[,1] == i), 2]
  M3_1 <- append(M3_1,mean(times3))
  M3_2 <- append(M3_2,var(times3))
  M3_3 <- append(M3_3,PerformanceAnalytics::skewness(times3))
  M3_4 <- append(M3_4,PerformanceAnalytics::kurtosis(times3,method="moment"))
  M1_5 <- append(M1_5, median(times1))
  M2_5 <- append(M2_5, median(times2))
  M3_5 <- append(M3_5, median(times3))
  M1_6 <- append(M1_6, IQR(times1))
  M2_6 <- append(M2_6, IQR(times2))
  M3_6 <- append(M3_6, IQR(times3))
}
  
# COLOUR BLIND COLOUR SCHEME
# Okabe-Ito colorblind safe palette
cb_palette <- c(
  "#000000", # black
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999"  # grey (optional)
)

# Now plot the legend separately
# Example color palette (replace with your actual cb_palette)
cb_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                "#FFFF33", "#A65628", "#F781BF", "#999999")

# Save plot with 600 DPI resolution
png("Fig1.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(M1_1, col = cb_palette[1], type = "l", lty = 2, lwd = 2, ylim = c(50, 160), xlab = "Deposited dose (number of Legionella)", ylab = "Time (hours)", main = "Mean incubation period")
lines(M2_1, col = cb_palette[2], lwd = 2, lty = 2)
lines(M3_1, col = cb_palette[3], lwd = 2, lty = 2)
lines(m1_1, col = cb_palette[4], lwd = 2, lty = 1)
lines(m2_1, col = cb_palette[5], lwd = 2)
lines(m3_1, col = cb_palette[6], lwd = 2)
lines(m4_1, col = cb_palette[7], lwd = 2)
lines(m5_1, col = cb_palette[8], lwd = 2)
lines(m6_1, col = cb_palette[9], lwd = 2)
lines(M1_1, col = cb_palette[1], lwd = 2, lty = 2)
lines(M2_1, col = cb_palette[2], lwd = 2, lty = 2)
lines(M3_1, col = cb_palette[3], lwd = 2, lty = 2)

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# Save plot with 600 DPI resolution
png("Fig2.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(M1_2^0.5, col = cb_palette[1], type = "l", lty = 2, lwd = 2, ylim = c(0, 30), xlab = "Deposited dose", ylab = "Standard deviation", main = "Standard deviation")
lines(M2_2^0.5, col = cb_palette[2], lwd = 2, lty = 2)
lines(M3_2^0.5, col = cb_palette[3], lwd = 2, lty = 2)
lines(m1_2^0.5, col = cb_palette[4], lwd = 2, lty = 1)
lines(m2_2^0.5, col = cb_palette[5], lwd = 2)
lines(m3_2^0.5, col = cb_palette[6], lwd = 2)
lines(m4_2^0.5, col = cb_palette[7], lwd = 2)
lines(m5_2^0.5, col = cb_palette[8], lwd = 2)
lines(m6_2^0.5, col = cb_palette[9], lwd = 2)

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# Save plot with 600 DPI resolution
png("Fig3.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(M1_3, col = cb_palette[1], type = "l", lwd = 2, lty = 2, ylim = c(0, 6), xlab = "Deposited dose (number of Legionella)", ylab = "Skewness", main = "Skewness")
lines(M2_3, col = cb_palette[2], lwd = 2, lty = 2)
# lines(M3_3, col = cb_palette[3], lwd = 2, lty = 2)
lines(m1_3, col = cb_palette[4], lwd = 2, lty = 1)
lines(m2_3, col = cb_palette[5], lwd = 2)
lines(m3_3, col = cb_palette[6], lwd = 2)
lines(m4_3, col = cb_palette[7], lwd = 2)
lines(m5_3, col = cb_palette[8], lwd = 2)
lines(m6_3, col = cb_palette[9], lwd = 2)

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# Save plot with 600 DPI resolution
png("Fig4.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(M1_4^0.5, col = cb_palette[1], type = "l", lwd = 2, lty = 2, ylim = c(0, 10), xlab = "Deposited dose (number of Legionella)", ylab = "Kurtosis", main = "Kurtosis")
lines(M2_4^0.5, col = cb_palette[2], lwd = 2, lty = 2)
# lines(M3_4^0.5, col = cb_palette[3], lwd = 2, lty = 2)
lines(m1_4^0.5, col = cb_palette[4], lwd = 2, lty = 1)
lines(m2_4^0.5, col = cb_palette[5], lwd = 2)
lines(m3_4^0.5, col = cb_palette[6], lwd = 2)
lines(m4_4^0.5, col = cb_palette[7], lwd = 2)
lines(m5_4^0.5, col = cb_palette[8], lwd = 2)
lines(m6_4^0.5, col = cb_palette[9], lwd = 2)

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# Save plot with 600 DPI resolution
png("Fig5.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(M1_5, col = cb_palette[1], type = "l", lty = 2, lwd = 2, ylim = c(50, 160), xlab = "Deposited dose (number of Legionella)", ylab = "Time (hours)", main = "Median incubation period")
lines(M2_5, col = cb_palette[2], lwd = 2, lty = 2)
lines(M3_5, col = cb_palette[3], lwd = 2, lty = 2)
lines(m1_5, col = cb_palette[4], lwd = 2, lty = 1)
lines(m2_5, col = cb_palette[5], lwd = 2)
lines(m3_5, col = cb_palette[6], lwd = 2)
lines(m4_5, col = cb_palette[7], lwd = 2)
lines(m5_5, col = cb_palette[8], lwd = 2)
lines(m6_5, col = cb_palette[9], lwd = 2)

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# Save plot with 600 DPI resolution
png("Fig6.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(M1_6, col = cb_palette[1], type = "l", lty = 2, lwd = 2, ylim = c(0, 50), xlab = "Deposited dose (number of Legionella)", ylab = "Time (hours)", main = "Interquartile range")
lines(M2_6, col = cb_palette[2], lwd = 2, lty = 2)
lines(M3_6, col = cb_palette[3], lwd = 2, lty = 2)
lines(m1_6, col = cb_palette[4], lwd = 2, lty = 1)
lines(m2_6, col = cb_palette[5], lwd = 2)
lines(m3_6, col = cb_palette[6], lwd = 2)
lines(m4_6, col = cb_palette[7], lwd = 2)
lines(m5_6, col = cb_palette[8], lwd = 2)
lines(m6_6, col = cb_palette[9], lwd = 2)

# Add minor ticks (nx, ny = number of intervals between major ticks)
minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)

dev.off()

# Save legend as PNG
png(
  filename = "Legend.png",
  width = 8, height = 3, units = "in", res = 600
)

# Create a completely empty plot
plot.new()  # blank canvas (no box, no axes)

# Add only the legend
legend(
  "center",   # position
  legend = c(
    "Erlang averaged", "Erlang Poisson", "Erlang intermediate",
    "Burr averaged", "Burr Poisson", "Burr intermediate",
    "[7] model A", "[7] model B", "[7] model C"
  ),
  col = cb_palette[c(4, 5, 6, 7, 8, 9, 1, 2, 3)],
  lwd = 2,
  lty = c(1, 1, 1, 1, 1, 1, 2, 2, 2),
  ncol = 3,
  horiz = FALSE,
  bty = "n"  # removes the legend box
)

# Close the PNG device
dev.off()

# Now do interquartile range
iqr_1 <- c()
iqr_2 <- c()
iqr_3 <- c()

# Now do fancy SA type stuff
samples_1 <- round(rbeta(n = 1000, shape1 = 0.5, shape2 = 5) * 500)
samples_2 <- round(rbeta(n = 1000, shape1 = 0.7, shape2 = 5) * 500)
samples_3 <- round(rbeta(n = 1000, shape1 = 1, shape2 = 5) * 500)
samples_4 <- round(rbeta(n = 1000, shape1 = 2, shape2 = 5) * 500)

SAMPLES_1 <- c()
SAMPLES_2 <- c()
SAMPLES_3 <- c()
SAMPLES_4 <- c()

for (i in samples_1) {
  print(i)
  subset_vals <- Erlang_pois[Erlang_pois[,1] == i, 2]
  if(length(subset_vals) > 0){
    SAMPLES_1 <- c(SAMPLES_1, sample(subset_vals, size = 1000, replace = TRUE))
  }
}

for (i in samples_2) {
  print(i)
  subset_vals <- Erlang_pois[Erlang_pois[,1] == i, 2]
  if(length(subset_vals) > 0){
    SAMPLES_2 <- c(SAMPLES_2, sample(subset_vals, size = 1000, replace = TRUE))
  }
}

for (i in samples_3) {
  print(i)
  subset_vals <- Erlang_pois[Erlang_pois[,1] == i, 2]
  if(length(subset_vals) > 0){
    SAMPLES_3 <- c(SAMPLES_3, sample(subset_vals, size = 1000, replace = TRUE))
  }
}

for (i in samples_4) {
  print(i)
  subset_vals <- Erlang_pois[Erlang_pois[,1] == i, 2]
  if(length(subset_vals) > 0){
    SAMPLES_4 <- c(SAMPLES_4, sample(subset_vals, size = 1000, replace = TRUE))
  }
}

# Save plot with 600 DPI resolution
png("Fig10.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(seq(0, 1, by = 0.001) * 500, dbeta(seq(0, 1, by = 0.001), shape1 = 0.5, shape2 = 5), type = "l", lwd = 2, xlab = "Deposited dose", ylab = "Density", main = "Deposited-dose distributions")
lines(seq(0, 1, by = 0.001) * 500, dbeta(seq(0, 1, by = 0.001), shape1 = 0.7, shape2 = 5), lwd = 2, col = 2)
lines(seq(0, 1, by = 0.001) * 500, dbeta(seq(0, 1, by = 0.001), shape1 = 1, shape2 = 5), lwd = 2, col = 3)
lines(seq(0, 1, by = 0.001) * 500, dbeta(seq(0, 1, by = 0.001), shape1 = 2, shape2 = 5), lwd = 2, col = 4)
legend("topright",
       legend = c(
         expression(paste(alpha, " = 0.5, ", beta, " = 5")),
         expression(paste(alpha, " = 0.7, ", beta, " = 5")),
         expression(paste(alpha, " = 1, ", beta, " = 5")),
         expression(paste(alpha, " = 2, ", beta, " = 5"))
       ),
       col = c(1, 2, 3, 4))

dev.off()


DENS1 <- density(SAMPLES_1)
DENS2 <- density(SAMPLES_2)
DENS3 <- density(SAMPLES_3)
DENS4 <- density(SAMPLES_4)


# Save plot with 600 DPI resolution
png("Fig11.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(DENS1, type = "l", lwd = 2, xlab = "Time in hours", ylab = "Density", main = "Erlang incubation-period distributions", ylim = c(0, 0.04))
lines(DENS2, lwd = 2, col = 2)
lines(DENS3, lwd = 2, col = 3)
lines(DENS4, lwd = 2, col = 4)
legend("topright",
       legend = c(
         expression(paste(alpha, " = 0.5, ", beta, " = 5")),
         expression(paste(alpha, " = 0.7, ", beta, " = 5")),
         expression(paste(alpha, " = 1, ", beta, " = 5")),
         expression(paste(alpha, " = 2, ", beta, " = 5"))
       ),
       col = c(1, 2, 3, 4), lwd = 2)

dev.off()

SAMPLES_1 <- c()
SAMPLES_2 <- c()
SAMPLES_3 <- c()
SAMPLES_4 <- c()

for (i in samples_1) {
  print(i)
  subset_vals <- Burr_pois[Burr_pois[,1] == i, 2]
  if(length(subset_vals) > 0){
    SAMPLES_1 <- c(SAMPLES_1, sample(subset_vals, size = 1000, replace = TRUE))
  }
}

for (i in samples_2) {
  print(i)
  subset_vals <- Burr_pois[Burr_pois[,1] == i, 2]
  if(length(subset_vals) > 0){
    SAMPLES_2 <- c(SAMPLES_2, sample(subset_vals, size = 1000, replace = TRUE))
  }
}

for (i in samples_3) {
  print(i)
  subset_vals <- Burr_pois[Burr_pois[,1] == i, 2]
  if(length(subset_vals) > 0){
    SAMPLES_3 <- c(SAMPLES_3, sample(subset_vals, size = 1000, replace = TRUE))
  }
}

for (i in samples_4) {
  print(i)
  subset_vals <- Burr_pois[Burr_pois[,1] == i, 2]
  if(length(subset_vals) > 0){
    SAMPLES_4 <- c(SAMPLES_4, sample(subset_vals, size = 1000, replace = TRUE))
  }
}

DENS1 <- density(SAMPLES_1)
DENS2 <- density(SAMPLES_2)
DENS3 <- density(SAMPLES_3)
DENS4 <- density(SAMPLES_4)


# Save plot with 600 DPI resolution
png("Fig12.png", width = 8, height = 6, units = "in", res = 600)

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

plot(DENS1, type = "l", lwd = 2, xlab = "Time in hours", ylab = "Density", main = "Burr incubation-period distributions", ylim = c(0, 0.04))
lines(DENS2, lwd = 2, col = 2)
lines(DENS3, lwd = 2, col = 3)
lines(DENS4, lwd = 2, col = 4)
legend("topright",
       legend = c(
         expression(paste(alpha, " = 0.5, ", beta, " = 5")),
         expression(paste(alpha, " = 0.7, ", beta, " = 5")),
         expression(paste(alpha, " = 1, ", beta, " = 5")),
         expression(paste(alpha, " = 2, ", beta, " = 5"))
       ),
       col = c(1, 2, 3, 4), lwd = 2)

dev.off()


# Time - dose correlation analysis
cor(Erlang_averaged$L, Erlang_averaged$time, method = "spearman")
cor(Erlang_pois$L, Erlang_pois$time, method = "spearman")
cor(Erlang_intermediate$L, Erlang_intermediate$time, method = "spearman")
cor(Burr_averaged$L, Burr_averaged$time, method = "spearman")
cor(Burr_pois$L, Burr_pois$time, method = "spearman")
cor(Burr_intermediate$L, Burr_intermediate$time, method = "spearman")

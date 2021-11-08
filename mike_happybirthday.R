invlogit <- function(x) exp(x) / (1+exp(x))

w1 <- "HAPPY BIRTHDAY!"
w2 <- "MIKE GLOUDEMANS!"

current_members <- c("Nathan Abell",
                     "Olivia de Goede",
                     "Marianne DeGorter",
                     "Tiffany Eulalio",
                     "Nicole Ersaro",
                     "Nicole Gay",
                     "Mike Gloudemans",
                     "Pagé Goddard",
                     "Emily Greenwald",
                     "Tanner Jensen",
                     "Stephen Montgomery",
                     "Daniel Nachun",
                     "Kameron Rodrigues",
                     "Jarod Rutledge",
                     "Kevin Smith",
                     "Nikki Teran",
                     "Rachel Ungar",
                     "Nikolai Gates Vetr",
                     "Andrew Marderstein"
)

first_names <- do.call(rbind, strsplit(current_members, " "))[,1]
first_names <- c("Mike", "Jarod", "Stephen")

h = 250
w = 500
png(filename = "w1.png", width = w, height = h)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0,1))
text(0,0.5, pos = 4, label = w1, cex = 3, font  = 2)
dev.off()

png(filename = "w2.png", width = w, height = h)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0,1))
text(0,0.5, pos = 4, label = w2, cex = 3, font  = 2)
dev.off()

w1c <- png::readPNG("w1.png")
w2c <- png::readPNG("w2.png")

file.remove("w1.png")
file.remove("w2.png")

w1c <- which(w1c == 0, arr.ind = T)[,2:1]
w1c[,2] <- h-w1c[,2]
w2c <- which(w2c == 0, arr.ind = T)[,2:1]
w2c[,2] <- h-w2c[,2]

#jitter points
w1c[,1] <- w1c[,1] + rnorm(nrow(w1c))
w2c[,1] <- w2c[,1] + rnorm(nrow(w2c))
w1c[,2] <- w1c[,2] + rnorm(nrow(w1c))
w2c[,2] <- w2c[,2] + rnorm(nrow(w2c))

#rescale so total range is the same
w1c[,1] <- w1c[,1] - min(w1c[,1])
w1c[,1] <- w1c[,1] / max(w1c[,1]) * 365.25
w2c[,1] <- w2c[,1] - min(w2c[,1])
w2c[,1] <- w2c[,1] / max(w2c[,1]) * max(w1c[,1])

#downsample longer word
w2c <- w2c[sample(1:nrow(w2c), size = nrow(w1c), replace = F),]

#sort both words and swap in xlocs from word 1
w1c <- w1c[order(w1c[,1]),]
w2c[,1] <- w2c[,1]

cols <- c('blue4','skyblue','darkgreen','orange','red')
alpha = 0.9
cols_w1 <- cols[ceiling(invlogit((w1c[,2] - mean(w1c[,2]) * 1.078) / sd(w1c[,2]) + rnorm(nrow(w1c))) * length(cols))]
cols_w2 <- cols[ceiling(invlogit((w2c[,2] - mean(w2c[,2]) * 1.05) / sd(w2c[,2]) + rnorm(nrow(w2c))) * length(cols))]

r <- 0.999
w1w2r <- cbind(w1c[,2], w2c[,2]) %*% (chol(diag(2) + r - diag(2)*r))
# w1w2r[,2] <- (w1w2r[,2] - mean(w1w2r[,2])) * 0.5 + mean(w1w2r[,2])
# w1w2r[,1] <- w1w2r[,1] / sd(w1w2r[,1])
# w1w2r[,2] <- w1w2r[,2] / sd(w1w2r[,2])
w1w2r[,2] <- w1w2r[,2] + rnorm(nrow(w1w2r), 0, sd = abs(w1w2r[,1] - max(w1w2r[,1])) / 5)
w1w2r[,1] <- w1w2r[,1] + rnorm(nrow(w1w2r), 0, sd = abs(w1w2r[,2] - max(w1w2r[,2])) / 5)
w1w2r <- exp(w1w2r)
w1w2r[,1] <- (w1w2r[,1] - mean(w1w2r[,1])) / sd(w1w2r[,1])
w1w2r[,2] <- (w1w2r[,2] - mean(w1w2r[,2])) / sd(w1w2r[,2])
cols_w1w2r <- cols[ceiling(invlogit((apply(w1w2r,1,prod) - mean(apply(w1w2r,1,prod)) * 50) / 
                              sd(apply(w1w2r,1,prod)) + rnorm(nrow(w1w2r))) * length(cols))]

pchs_w1w2r <- c(rep(21, which.max(apply(w1w2r,1,prod)) - 1), 
                23, 
                rep(21, nrow(w1w2r) - which.max(apply(w1w2r,1,prod))))
cex_w1w2r <- c(rep(1.25, which.max(apply(w1w2r,1,prod)) - 1), 
                2, 
                rep(1.25, nrow(w1w2r) - which.max(apply(w1w2r,1,prod))))


w1c[,2] <- (w1c[,2] - min(w1c[,2])) / diff(range(w1c[,2])) * max(w1w2r[,1])
w2c[,2] <- (w2c[,2] - min(w2c[,2])) / diff(range(w2c[,2])) * max(w1w2r[,2])

#now do the actual plotting
layout( matrix(c(1,1,3,2,2,3), ncol=3, byrow = T))
par(mar = c(4,5,2,4))

plot(w1c, pch = 21, col = adjustcolor(1, alpha), bg = adjustcolor(cols_w1, alpha),
     xlab = "DAY OF THE YEAR (DAYS)",
     ylab = latex2exp::TeX(paste0("$-log_{10}(", w1, ")$")))
plot(w2c, pch = 21, col = adjustcolor(1, alpha), bg = adjustcolor(cols_w2, alpha),
     xlab = "DAY OF THE YEAR (DAYS)",
     ylab = latex2exp::TeX(paste0("$-log_{10}(", w2, ")$")))
plot(w1w2r, pch = pchs_w1w2r, col = adjustcolor(1, alpha), bg = adjustcolor(cols_w1w2r, alpha), 
     cex = cex_w1w2r, xlab = latex2exp::TeX(paste0("$-log_{10}(", w1, ")$")),
     ylab = latex2exp::TeX(paste0("$-log_{10}(", w2, ")$")))
arrows(y0 = mean(par("usr")[3:4]), y1 = mean(par("usr")[3:4]),
       x0 = -diff(par("usr")[1:2]) / 2, x1 = -diff(par("usr")[1:2]) / 3.5, xpd = NA, lwd = 5)

xl <- par("usr")[1] + diff(par("usr")[1:2]) * 0.85
yb <- par("usr")[3] + diff(par("usr")[3:4]) * 0.05
xr <- par("usr")[1] + diff(par("usr")[1:2]) * 0.9
yt <- par("usr")[3] + diff(par("usr")[3:4]) * 0.275

rect(
  xl,
  head(seq(yb,yt,(yt-yb)/length(cols)),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/length(cols)),-1),
  col=cols
)
text(1:4*0.2, y = seq(yb,yt,(yt-yb)/length(cols))[2:5], x = xr, pos = 4)
text("BIRTHDAY POWER", y = seq(yb,yt,(yt-yb)/length(cols))[3], srt = 90, x = xl - (xr-xl)*0.35, pos = 3, cex = 0.7)
text(first_names, x = w1w2r[order(apply(w1w2r, 1, prod), decreasing = T)[1:3],][,1],
     y = w1w2r[order(apply(w1w2r, 1, prod), decreasing = T)[1:3],][,2], pos = 2, cex = 0.85)

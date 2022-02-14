# Date: 21/12/4
#---- Setting E.1 ----
p = 400
b = 70 * 2
K = 3
n = K * b

mu_1 = rep(0, p)
mu_2 = rep(0, p)
mu_3 = rep(0, p)
mu_1[1:3] <- -2
mu_2[4:6] <- 2
mu = cbind(mu_1, mu_2, mu_3)

rou = 0.8
Sigma = matrix(rep(1, p * p), nrow = p)
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j] = rou ^ abs(i - j)
  }
}

U <- mu.to.U(mu)
W.oracle = solve(Sigma) %*% U

#---- model fit ----
s = 40
M.list <- c(20, 30, 40, 50, 60)
mcr.oracle <- matrix(rep(0, 5 * s), nrow = s)
mcr.cen <- matrix(rep(0, 5 * s), nrow = s)
mcr.local <- matrix(rep(0, 5 * s), nrow = s)
mcr.dmSLDA <- matrix(rep(0, 5 * s), nrow = s)

error.cen <- matrix(rep(0, 5 * s), nrow = s)
error.dmSLDA <- matrix(rep(0, 5 * s), nrow = s)
error.local <- matrix(rep(0, 5 * s), nrow = s)

for(j in 1:5){
  M = M.list[j]
  for(i in 1:s){
    data <- data.generation(M, b, mu, Sigma)
    X = matrix(rep(0, K * b * M * p), ncol = p)
    for(k in 1:K){
      for(m in 1:M){
        X[((k-1)*M*b + (m-1)*b + 1):((k-1)*M*b + m*b),] = data[[3]][[m]][((k-1)*b + 1):(k*b),]
      }
    }
    data.test <- data.generation(M = 1, b = 300, mu, Sigma)[[3]][[1]]
    
    W.cen <- cen.SLDA(data, K = 3)
    result <- dmSLDA(data, K = 3, whole.path = TRUE)
    W.local = result[[1]][[1]]
    score = result[[4]]
    W.dmSLDA = result[[1]][[which.min(score)]]
    
    mcr.oracle[i, j] <- classifier(W.oracle, X, data.test)
    mcr.cen[i, j] <- classifier(W.cen, X, data.test)
    mcr.local[i, j] <- classifier(W.local, X, data.test)
    mcr.dmSLDA[i, j] <- classifier(W.dmSLDA, X, data.test)
    
    error.cen[i, j] <- sqrt(sum((W.cen - W.oracle)^2))
    error.dmSLDA[i, j] <- sqrt(sum((W.dmSLDA - W.oracle)^2))
    error.local[i, j] <- sqrt(sum((W.local - W.oracle)^2))
  }
}


#---- df.supp.1 ----
mcr.mean.oracle <- apply(mcr.oracle, 2, mean)
mcr.mean.cen <- apply(mcr.cen, 2, mean)
mcr.mean.dmSLDA <- apply(mcr.dmSLDA, 2, mean)
mcr.mean.local <- apply(mcr.local, 2, mean)

mcr.sd.oracle <- apply(mcr.oracle, 2, sd)
mcr.sd.cen <- apply(mcr.cen, 2, sd) 
mcr.sd.local <- apply(mcr.local, 2, sd) 
mcr.sd.dmSLDA <- apply(mcr.dmSLDA, 2, sd)

mcr.mean <- rep(0, 20)
mcr.sd <- rep(0, 20)
mcr.mean[1:5] <- mcr.mean.local
mcr.mean[6:10] <- mcr.mean.dmSLDA
mcr.mean[11:15] <- mcr.mean.cen
mcr.mean[16:20] <- mcr.mean.oracle
mcr.sd[1:5] <- mcr.sd.local
mcr.sd[6:10] <- mcr.sd.dmSLDA
mcr.sd[11:15] <- mcr.sd.cen
mcr.sd[16:20] <- mcr.sd.oracle

methods <- factor(c(rep('Local', 5), rep('dmSLDA', 5), rep('Centralized', 5), rep('Oracle', 5)), levels = c('Local', 'dmSLDA', 'Centralized', 'Oracle'))
mcr.M <- rep(M.list, 4)

df.plot.supp.1 <- data.frame(methods,
                             mcr.mean,
                             mcr.sd,
                             M = mcr.M)
save(df.plot.supp.1, file = "Numerical/data/df_supp_1.RData")

#---- df.supp.2 ----
error.mean.cen <- apply(error.cen, 2, mean)
error.mean.dmSLDA <- apply(error.dmSLDA, 2, mean)
error.mean.local <- apply(error.local, 2, mean)

error.sd.cen <- apply(error.cen, 2, sd)
error.sd.dmSLDA <- apply(error.dmSLDA, 2, sd) 
error.sd.local <- apply(error.local, 2, sd)

error.mean <- rep(0, 15)
error.sd <- rep(0, 15)
error.mean[1:5] <- error.mean.cen
error.mean[6:10] <- error.mean.dmSLDA
error.mean[11:15] <- error.mean.local
error.sd[1:5] <- error.sd.cen
error.sd[6:10] <- error.sd.dmSLDA
error.sd[11:15] <- error.sd.local

error.methods <- factor(c(rep('Centralized', 5), rep('dmSLDA', 5), rep('Local', 5)), levels = c('Local', 'dmSLDA', 'Centralized', 'Oracle'))
error.M <- rep(M.list, 3)

df.plot.supp.2 <- data.frame(error.mean,
                             error.sd,
                             methods = error.methods,
                             M = error.M)
save(df.plot.supp.2, file = "Numerical/data/df_supp_2.RData")

#---- plot supp.1 ----
p.supp.1 <- ggplot(data = df.plot.supp.1, aes(x = M, y = mcr.mean)) +
  geom_point(aes(color = methods, shape = methods, fill = methods)) +
  geom_line(aes(color = methods)) +
  theme_article() +
  scale_color_manual(values = c("dmSLDA" = "red", "Local" = "green", "Centralized" = "blue", "Oracle" = "purple"),labels = c("Local","dmSLDA","Centralized","Oracle")) +
  scale_fill_manual(values = c("dmSLDA" = "red", "Local" = "green", "Centralized" = "blue", "Oracle" = "purple"),labels = c("Local","dmSLDA","Centralized","Oracle")) +
  scale_shape_manual(values = c("dmSLDA" = 21, "Local" = 22, "Centralized" = 24, "Oracle" = 25),labels = c("Local","dmSLDA","Centralized","Oracle")) +
  ggtitle(TeX("$\\sigma = 0.8$")) +
  theme(legend.position = "right",
        legend.title=element_blank(),
        legend.background = element_rect(color = "#666666", size = 0.1)
  )+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        title = element_text(size = 13),
        legend.text = element_text(size = 13))

#---- plot supp.2 ----
p.supp.2 <- ggplot(data = df.plot.supp.2, aes(x = M, y = error.mean)) +
  geom_point(aes(color = methods, shape = methods, fill = methods)) +
  geom_line(aes(color = methods)) +
  geom_ribbon(aes(ymin = error.mean - error.sd,
                  ymax = error.mean + error.sd,
                  fill = methods), alpha = 0.1) +
  ylab("estimation error") +
  theme_article() +
  scale_color_manual(values = c("dmSLDA" = "red", "Local" = "green", "Centralized" = "blue", "Oracle" = "purple"),labels = c("Local","dmSLDA","Centralized")) +
  scale_fill_manual(values = c("dmSLDA" = "red", "Local" = "green", "Centralized" = "blue", "Oracle" = "purple"),labels = c("Local","dmSLDA","Centralized")) +
  scale_shape_manual(values = c("dmSLDA" = 21, "Local" = 22, "Centralized" = 24, "Oracle" = 25),labels = c("Local","dmSLDA","Centralized")) +
  ggtitle(TeX("$\\sigma = 0.8$")) +
  theme(legend.position = "right",
        legend.title=element_blank(),
        legend.background = element_rect(color = "#666666", size = 0.1)
  )+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        title = element_text(size = 13),
        legend.text = element_text(size = 13))

#---- E.2 ----
# df.plot.supp.2.2
K = 2
# setting
p = 200
b = 200
n = K * b
mu.1 = c(rep(0.5, 10), rep(0, p - 10))
mu.2 = - mu.1
mu = cbind(mu.1, mu.2)
Sigma = matrix(rep(1, p * p), nrow = p)
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j] = 0.8 ^ abs(i - j)
  }
}

U <- mu.to.U(mu)
W.oracle = solve(Sigma) %*% U

# experiments
s = 40
M.list <- c(5, 10, 15, 20, 30)

mcr.cen = matrix(rep(0, 5 * s), nrow = s)
mcr.dmSLDA = matrix(rep(0, 5 * s), nrow = s)
mcr.dSLDA = matrix(rep(0, 5 * s), nrow = s)
mcr.dsvm = matrix(rep(0, 5 * s), nrow = s)
for(j in 1:5){
  M = M.list[j]
  for(i in 1:s){
    data = data.generation(M, b, mu, Sigma)
    # X: data.aggregated
    X = matrix(rep(0, K * b * M * p), ncol = p)
    for(k in 1:K){
      for(m in 1:M){
        X[((k-1)*M*b + (m-1)*b + 1):((k-1)*M*b + m*b),] = data[[3]][[m]][((k-1)*b + 1):(k*b),]
      }
    }
    data.test <- data.generation(M = 1, b = 300, mu, Sigma)[[3]]
    mu.bar = apply(X, 2, mean)
    
    W.cen <- cen.SLDA(data, K = 2)
    W.dmSLDA <- dmSLDA(data, K = 2)
    beta.dSLDA = dSLDA.bin(data)
    beta.dSVM = dSVM(data)
    
    mcr.dmSLDA[i,j] = classifier(W.dmSLDA, X, data.test[[1]])
    mcr.cen[i,j] = classifier(W.cen, X, data.test[[1]])
    mcr.dSLDA[i,j] = classifier.bi(beta.dSLDA, mu.bar, data.test[[1]]) / nrow(data.test[[1]])
    mcr.dsvm[i,j] = classifier.svm(beta.dSVM, data.test[[1]])/nrow(data.test[[1]])
  }
}

mcr.mean.cen <- apply(mcr.cen, 2, mean)
mcr.mean.dmSLDA <- apply(mcr.dmSLDA, 2, mean)
mcr.mean.dSLDA <- apply(mcr.dSLDA, 2, mean)
mcr.mean.dsvm <- apply(mcr.dsvm, 2, mean)
mcr.sd.cen <- apply(mcr.cen, 2, sd) 
mcr.sd.dSLDA <- apply(mcr.dSLDA, 2, sd) 
mcr.sd.dmSLDA <- apply(mcr.dmSLDA, 2, sd)
mcr.sd.dsvm <- apply(mcr.dsvm, 2, sd)

mcr.mean <- rep(0, 20)
mcr.sd <- rep(0, 20)

mcr.mean[1:5] <- mcr.mean.dsvm
mcr.mean[6:10] <- mcr.mean.dSLDA
mcr.mean[11:15] <- mcr.mean.dmSLDA
mcr.mean[16:20] <- mcr.mean.cen
mcr.sd[1:5] <- mcr.sd.dsvm
mcr.sd[6:10] <- mcr.sd.dSLDA
mcr.sd[11:15] <- mcr.sd.dmSLDA
mcr.sd[16:20] <- mcr.sd.cen

methods <- factor(c(rep("dSVM",5), rep('dSLDA', 5), rep('dmSLDA', 5), rep('Cen.SLDA', 5)), levels = c("dSVM", 'dSLDA', 'dmSLDA', 'Cen.SLDA'))
mcr.M <- rep(M.list, 4)

df.plot.supp.5 <- data.frame(methods,
                             mcr.mean,
                             mcr.sd,
                             M = mcr.M)
save(df.plot.supp.5, file = "Numerical/data/df_supp_5_1.RData")

p.supp.5 <- ggplot(data = df.plot.supp.5, aes(x = M, y = mcr.mean)) +
  geom_point(aes(color = methods, shape = methods, fill = methods)) +
  geom_line(aes(color = methods)) +
  theme_article() +
  scale_color_manual(values = c("dSVM" = "#003300", "dmSLDA" = "red", "Cen.SLDA" = "blue", "dSLDA" = "#663366"),labels = c("dSVM","dSLDA","dmSLDA","Cen.SLDA")) +
  scale_fill_manual(values = c("dSVM" = "#003300", "dmSLDA" = "red", "Cen.SLDA" = "blue", "dSLDA" = "#663366"),labels = c("dSVM","dSLDA","dmSLDA","Cen.SLDA")) +
  scale_shape_manual(values = c("dSVM" = 4, "dmSLDA" = 21, "Cen.SLDA" = 24, "dSLDA" = 23),labels = c("dSVM","dSLDA","dmSLDA","Cen.SLDA")) +
  ggtitle("binary") +
  theme(legend.position = "right",
        legend.title=element_blank(),
        legend.background = element_rect(color = "#666666", size = 0.3)
  )+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        title = element_text(size = 27),
        legend.text = element_text(size = 27)) +
  scale_x_continuous(limits=c(5,30), breaks=seq(0,30,5))

ggsave(p.supp.5, file = "fig.supp.5.pdf", width = 10, height = 6)


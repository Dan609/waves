### PCA

calcium2 <- cbind(calcium$ROI553,
                  calcium$ROI461,
                  calcium$ROI418,
                  calcium$ROI333,
                  calcium$ROI222,
                  calcium$ROI148)

pairs(calcium2)

arc.pca1 <- princomp(calcium2, scores = TRUE, cor = TRUE)

arc.pca2 <- prcomp(calcium2)

summary(arc.pca2)

par(mfrow=c(1,1))

plot(arc.pca2)

biplot(arc.pca2)

arc.pca1$loadings

arc.pca1$scores

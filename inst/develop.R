data("wang2019")

res = deft_do(wang2019, group_level = c("Male", "Female"), method = "FE")

forest(res$all$model,
       slab = res$all$data$entry, atransf = exp, xlab="Hazard ratio", showweights = TRUE,
       mlab = "overall")
op = par(cex = 0.75, font = 2)
text(-13, 7.5, "Trial(s) and subgroup", pos = 4)
text(6, 15, "Relative Risk [95% CI]", pos = 2)
par(op)

forest(res$subgroup$model, showweights = TRUE, atransf = exp)


res3 = rma(yi ~ trial + 0, sei = sei,
           data = res$subgroup$data)

summary(res)

forest(res)

funnel(res)


res2 <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg,
           slab=paste(author, year, sep=", "))
forest(res2)
op <- par(xpd=TRUE)
text(x=-8.5, y=-1:16, -1:16, pos=4, cex=.5)
par(op)


res$subgroup$data
fit <- rma(yi = yi, sei=sei, data = res$subgroup$data)
forest(fit, atransf = exp)

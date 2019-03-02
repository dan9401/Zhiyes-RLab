require(survival)
data(pbc)

surObj = Surv(pbc$time, pbc$status)

model = survfit(surObj ~ 1)
summary(model)
model

plot(model, ylab = "Survival Rate", xlab = "Days Elapsed")

survdiff(surObj ~ pbc$trt)

coxModel = coxph(surObj ~ pbc$bili)
zphModel = cox.zph(coxModel)
coxModel
zphModel

plot(zphModel)

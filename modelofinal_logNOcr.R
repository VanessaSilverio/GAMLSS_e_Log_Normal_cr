m1 <-gamlss(Surv(diassobrev, ULTINFO) ~ 1,
            family= cens('LOGNO2cr'),data=dados1)
nC <- detectCores(); nC
m1 <- stepGAIC(m1, scope=expl, parallel = 'multicore', ncpus = nC, direction = 'forward')
m1 <- stepGAIC(m1, scope=expl, parallel = 'multicore', ncpus = nC, direction = 'forward', what = 'nu') 
~QUIMIO + ECGRUP + pb(IDADE) + CIRURGIA + RADIO + pb(datatrat) 
m1 <- stepGAIC(m1, scope=expl, parallel = 'multicore', ncpus = nC, direction = 'forward', what = 'sigma')
~QUIMIO + RADIO + CIRURGIA + pb(IDADE)

mLOGNO <- gamlss(Surv(diassobrev, ULTINFO) ~ ECGRUP + pb(datatrat) + CIRURGIA + 
                   QUIMIO + RADIO + pb(IDADE) + SEXO + M, # removi T
                 sigma.formula = ~QUIMIO + RADIO + CIRURGIA, # removi idade
                 nu.formula = ~QUIMIO + ECGRUP + IDADE +  
                   datatrat, # removi radio e cirurgia
                 family = cens("LOGNOcr"), 
                 data = dados1)
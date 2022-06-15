
#PARÂMETROS DA SIMULAÇÃO
Tamanho.amostras = 50
Densidade.y = "Assimetrico" # "simetrico"
Volume.dos.dados = "Big" # "Nao.Big"
Presenca.de.Multicolinearidade = "X.indep" # X.multicol





gera.dados = function(n, tipo.y, tipo.x, relacao){
  library(car)
  
X.indep = data.frame(x1= rnorm(n,4,0.2),
             x2= (rexp(n,1)),
             x3= rbinom(n,4,0.2),
             x4= rpois(n,1),
             x5= sample(c("0","1"),n, prob = c(0.5,0.5), replace = T),
             x6= rnorm(n,7,0.5),
             x7= (rexp(n,2)),
             x8= rbinom(n,10,0.5),
             x9=  sample(c("0","1"),n, prob = c(0.3,0.7), replace = T),
             x10= rchisq(n,1),
             x11= rnorm(n,11,1),
             x12= (rnorm(n,1,0.1)),
             x13= rbinom(n,2,0.5),
             x14= rpois(n,2),
             x15= rchisq(n,1))

x1= rnorm(n,4,0.2)
x2= x1*(rexp(n,1))
x3= x2*rbinom(n,4,0.2)
x4= x3*rpois(n,1)
x5= ifelse(x4 > mean(x4),"1","0")
x6= x4*rnorm(n,7,0.5)
x7= x6*(rexp(n,2))
x8= x7*rbinom(n,10,0.5)
x9= ifelse(x8 > quantile(x8, prob= 0.75),"1","0")
x10= x8*rchisq(n,1)
x11= x10*rnorm(n,11,1)
x12= x11*(rnorm(n,1,0.1))
x13= x12*rbinom(n,2,0.5)
x14= x13*rpois(n,2)
x15= x14*rchisq(n,1)
X.multicol = data.frame(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15)

X = list(X.indep = X.indep, X.multicol = X.indep)
X = X[[relacao]]

C5 = ifelse(X[,5]=='1', 1, (-1))
C9 = ifelse(X[,9]=='1', (-3), 3)
y1 = 2*(X[,1]) + (1/10)*X[,2] + X[,3] -X[,4] + C5 + rnorm(n,2,0.5)
y2 = 2*(X[,1]) + (1/10)*X[,2] + X[,3] -X[,4] + C5 + rgamma(n,2,4)
y3 = 3*(X[,1]) + (1/10)*X[,2] + X[,3] -X[,4] + C5 + X[,6]*X[,7] + X[,8] + C9 - 
  X[,10] + (X[,11]^-1) + X[,12] + X[,13]^2 - X[,14] + X[,15]*0.2 + rnorm(n,20,0.05)
y4 = 2*(X[,1]) + (1/10)*X[,2] + X[,3] -X[,4] + C5 + X[,6]*X[,7] + X[,8] + C9 - 
  X[,10]*(X[,11]^-1) + X[,12] + X[,13]^2 - X[,14] + X[,15]*0.2 + rgamma(n,2,4)
y = ifelse(tipo.y=="simetrico" & tipo.x=="Big", data.frame(y3), ifelse(tipo.y=="Assimetrico" & tipo.x=="Big",data.frame(y4),
           ifelse(tipo.y=="Assimetrico" & tipo.x=="Nao.Big",data.frame(y2),data.frame(y1))))

y= y[[1]]

titulo = ifelse(tipo.y=="simetrico","Densidade de Y - N~(10,1)","Densidade de Y - 10(Exp~(2)")
  plot(density(y), main= titulo,xlab="",ylab="") ## Verificar assimetria
  abline(v= median(y), col= "red", lwd = 2)
  abline(v= quantile(y,probs=0.25), col= "red", lwd = 2)
  abline(v= quantile(y,probs=0.75), col= "red", lwd = 2)
  
  refe = ifelse(tipo.x=="Big", data.frame(1:15), data.frame(1:5))
  refe= refe[[1]]
  X= X[,refe]
  dados = X
  dados$y = y
  modelo.vif= lm(y~.,dados) 
  vif(modelo.vif)
  return(dados)
}
dados = gera.dados(n=Tamanho.amostras, tipo.y = Densidade.y,
                   tipo.x = Volume.dos.dados, relacao= Presenca.de.Multicolinearidade)


################################################
varResp = 'y'
Medidas = function(residuos,pred,teste, formula){
  ind = which(colnames(teste) == as.character(formula)[2])
  y = teste[, ind]
  estimados.medios= mean(pred)
  observados.medios= mean(y)
  estimados.desvios = sd(pred)
  observados.desvios = sd(y)
  r = cor(pred,y)
  EQM = mean((residuos)^2)
  numerador.Var = (estimados.medios - observados.medios)^2
  numerador.Vies = (estimados.desvios - observados.desvios)^2
  numerador.cov = 2*(1-r)*estimados.desvios*observados.desvios
  Um = numerador.Vies/ EQM
  Ur = numerador.Var/ EQM
  Uc = numerador.cov/EQM
  return(list(Um = Um,Ur = Ur,Uc = Uc))
}

formula = as.formula(paste(varResp, '~ .')) 

AuxSaida = function(pred, teste, formula, alg){
  ind = which(colnames(teste) == as.character(formula)[2])
  residuos = pred-teste[, ind]
  rs = Medidas(residuos,pred,teste,formula)
  rs[['Classificador']] = alg
  return(do.call(data.frame, rs))
}

ValidCruzada = function(cls, dados, formula, kfolds = 10, ...){
  library(dismo)
  id = kfold(1:nrow(dados), kfolds)
  vcs = data.frame()
  for(i in 1:kfolds){
    treino = dados[id != i, ]
    teste = dados[id == i, ]
    kcls = cls(treino, teste, formula) 
    vcs = rbind(vcs, kcls)
  }
  vcs = vcs[, c(ncol(vcs), 1:(ncol(vcs) - 1))]
  avg = aggregate(. ~ Classificador, data = vcs, FUN = mean)
  std = aggregate(. ~ Classificador, data = vcs, FUN = sd)
  return(list(Media = avg, Desvio = avg, Modelos = vcs))
}

ValidCruzadarep = function(cls, dados, formula, kfolds = 10, reps = 10){
  x = data.frame()
  for(i in 1:reps) x = rbind(x, ValidCruzada(cls, dados, formula, kfolds)$Media) 
  avg = aggregate(. ~ Classificador, data = x, FUN = mean)
  std = aggregate(. ~ Classificador, data = x, FUN = sd)
  return(list(Media = avg, Desvio = std, Modelos = x))
}
###### FASE I   
#SVR 
library(kernlab)
KERNEL1= "polydot"
SVR1 = function(treino, teste, formula){
  cls = ksvm(formula, data = treino, kernel = KERNEL1)
  pred = predict(cls, newdata = data.frame(teste))
  return(pred)
}
medSVR1 = function(treino, teste, formula){
  alg = 'SVR - polydot' 
  pred = SVR1(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}
KERNEL2= "rbfdot"
SVR2 = function(treino, teste, formula){
  cls = ksvm(formula, data = treino, kernel = KERNEL2)
  pred = predict(cls, newdata = data.frame(teste))
  return(pred)
}
medSVR2 = function(treino, teste, formula){
  alg = 'SVR - rbfdot' 
  pred = SVR2(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}
KERNEL3= "vanilladot"
SVR3 = function(treino, teste, formula){
  cls = ksvm(formula, data = treino, kernel = KERNEL3)
  pred = predict(cls, newdata = data.frame(teste))
  return(pred)
}
medSVR3 = function(treino, teste, formula){
  alg = 'SVR - vanilladot' 
  pred = SVR3(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}

#SpAM

SPAM = function(treino, teste, formula){
  library(SAM)
  ind.treino = which(colnames(treino) == as.character(formula)[2])
  ind.teste = which(colnames(teste) == as.character(formula)[2])
  dados.x.treino= (treino[,-ind.treino])
  dados.y.treino= (treino[,ind.treino])
  dados.x.teste = (teste[,-ind.teste])
  for(i in 1:ncol(dados.x.treino)) {
    dados.x.treino[,i] <- as.numeric(dados.x.treino[,i])
  }
  for(i in 1:ncol(dados.x.teste)) {
    dados.x.teste[,i] <- as.numeric(dados.x.teste[,i])
  }
  cls = samQL(dados.x.treino,dados.y.treino,p=1)
  pred = predict(object=cls, dados.x.teste)
  pred= matrix(pred$values,nrow= length(teste[,ind.teste]),
               ncol=(length(pred$values)/length(teste[,ind.teste])))
  pred = apply(pred,1, mean)
  return(pred)
}

medSPAM = function(treino, teste, formula){
  alg = 'SpAM'
  pred = SPAM(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}
#Ridge
library(glmnet)
RIDGE= function(treino, teste, formula){
  ind.treino = which(colnames(treino) == as.character(formula)[2])
  ind.teste = which(colnames(teste) == as.character(formula)[2])
  dados.x.treino= data.matrix(treino[,-ind.treino])
  dados.y.treino= data.matrix(treino[,ind.treino])
  dados.x.teste = data.matrix(teste[,-ind.teste])
  cls =  cv.glmnet(dados.x.treino,dados.y.treino,alpha = 0,type.measure = "mse")
  lambda= as.numeric(cls$lambda.1se)
  pred = predict(object=cls, s =lambda, newx =dados.x.teste)
  return(pred)
}
medRIDGE = function(treino, teste, formula){
  alg = 'RIDGE' 
  pred = RIDGE(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}

#LASSO
library(glmnet)
LASSO= function(treino, teste, formula){
  ind.treino = which(colnames(treino) == as.character(formula)[2])
  ind.teste = which(colnames(teste) == as.character(formula)[2])
  dados.x.treino= data.matrix(treino[,-ind.treino])
  dados.y.treino= data.matrix(treino[,ind.treino])
  dados.x.teste = data.matrix(teste[,-ind.teste])
  cls =  cv.glmnet(dados.x.treino,dados.y.treino,alpha = 1,type.measure = "mse")
  lambda= as.numeric(cls$lambda.1se)
  pred = predict(object=cls, s =lambda, newx =dados.x.teste)
  return(pred)
}
medLASSO = function(treino, teste, formula){
  alg = 'LASSO' 
  pred = LASSO(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}

#ELASTIC NET
library(glmnet)
E.NET= function(treino, teste, formula){
  ind.treino = which(colnames(treino) == as.character(formula)[2])
  ind.teste = which(colnames(teste) == as.character(formula)[2])
  dados.x.treino= data.matrix(treino[,-ind.treino])
  dados.y.treino= data.matrix(treino[,ind.treino])
  dados.x.teste = data.matrix(teste[,-ind.teste])
  cls =  cv.glmnet(dados.x.treino,dados.y.treino,alpha = 0.5,type.measure = "mse")
  lambda= as.numeric(cls$lambda.1se)
  pred = predict(object=cls, s =lambda, newx =dados.x.teste)
  return(pred)
}
medE.NET = function(treino, teste, formula){
  alg = 'ELASTIC NET' 
  pred = E.NET(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}


#RQ
library(quantreg)
RQ= function(treino, teste, formula){
  library(quantreg)
  cls =  rq(formula = formula, tau = 0.5, data = data.frame(treino))
  pred = predict(object=cls, newdata = data.frame(teste))
  return(pred)
}
medRQ = function(treino, teste, formula){
  alg = 'RQ' 
  pred = RQ(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}
#MLG
familia1 = gaussian
MLG1= function(treino, teste, formula,familia){
  cls =  glm(y~., data = data.frame(treino),family = familia1)
  pred =predict(object=cls, newdata = data.frame(teste))
  return(pred)
}
medMLG1 = function(treino, teste, formula){
  alg = 'MLG - Gaussian Family' 
  pred = MLG1(treino, teste, formula,familia)
  return(AuxSaida(pred, teste, formula, alg))
}

familia2 = Gamma
MLG2= function(treino, teste, formula,familia){
  cls =  glm(y~., data = data.frame(treino),family = familia2)
  pred =predict(object=cls, newdata = data.frame(teste))
  pred = 1/pred  
  return(pred)
}
medMLG2 = function(treino, teste, formula){
  alg = 'MLG - Gamma Family' 
  pred = MLG2(treino, teste, formula,familia)
  return(AuxSaida(pred, teste, formula, alg))
}

#POLINOMIAL ORTOGONAL 
formula.polinomial = function(dados){
  ref1 = c(rep('x', (length(dados)- 1)))
  ref2 = c(1:(length(dados)-1))
  ref2 = as.character(ref2)
  ref8 = paste0( ref1,ref2)
 return(ref8)
}
ref8 = formula.polinomial(dados)


Grau1 = 2
POLIORT1= function(treino, teste, formula,familia,ref8){
  treino= data.matrix(treino)
  teste = data.matrix(teste)
  formula.poli = as.formula(paste('y~poly(', paste(ref8, collapse=" + "),
                                  paste(", degree =", Grau1,")",sep="")))
  cls = lm(formula=formula.poli,data = data.frame(treino))
  pred = predict(cls, data.frame(teste))
  return(pred)
  
}
medPOLIORT1 = function(treino, teste, formula){
  alg = 'POLINOMIAL ORTOGONAL - Grau 2' 
  pred = POLIORT1(treino, teste, formula,familia,ref8)
  return(AuxSaida(pred, teste, formula, alg))
}
Grau2 = 3
POLIORT2= function(treino, teste, formula,familia,ref8){
  treino= data.matrix(treino)
  teste = data.matrix(teste)
  formula.poli = as.formula(paste('y~poly(', paste(ref8, collapse=" + "),
                                  paste(", degree =", Grau2,")",sep="")))
  cls = lm(formula=formula.poli,data = data.frame(treino))
  pred = predict(cls, data.frame(teste))
  return(pred)
}
medPOLIORT2 = function(treino, teste, formula){
  alg = 'POLINOMIAL ORTOGONAL- Grau 3' 
  pred = POLIORT2(treino, teste, formula,familia,ref8)
  return(AuxSaida(pred, teste, formula, alg))
}
Grau3 = 4
POLIORT3= function(treino, teste, formula,familia,ref8){
  treino= data.matrix(treino)
  teste = data.matrix(teste)
  formula.poli = as.formula(paste('y~poly(', paste(ref8, collapse=" + "),
                                  paste(", degree =", Grau3,")",sep="")))
  cls = lm(formula=formula.poli,data = data.frame(treino))
  pred = predict(cls, data.frame(teste))
  return(pred)
}
medPOLIORT3 = function(treino, teste, formula){
  alg = 'POLINOMIAL ORTOGONAL- Grau 4' 
  pred = POLIORT3(treino, teste, formula,familia,ref8)
  return(AuxSaida(pred, teste, formula, alg))
}

tabela.de.desempenho = matrix(NA, ncol = 5, nrow = 13)
colnames(tabela.de.desempenho) = c("Modelos", "Um","Ur", "Uc", "Medida")
tabela.de.desempenho = data.frame(tabela.de.desempenho)


medida  = medSVR1
nome.medida = "medSVR1"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[1,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3), 
                             round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)
medida  = medSVR2
nome.medida = "medSVR2"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[2,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3), 
                             round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)
medida  = medSVR3
nome.medida = "medSVR3"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[3,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3), 
                             round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)
medida  = medSPAM
nome.medida = "medSPAM"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[4,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3), 
                             round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)
medida  = medRIDGE
nome.medida = "medRIDGE"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[5,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3),
                             round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)
medida  = medLASSO
nome.medida = "medLASSO"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[6,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3), 
                             round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)

medida  = medE.NET
nome.medida = "medE.NET"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[7,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3), 
                             round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)

medida  = medRQ
nome.medida = "medRQ"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[8,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3), 
                             round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)
medida  = medMLG1
nome.medida = "medMLG1"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[9,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]], 3),
                             round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)
medida  = medMLG2
nome.medida = "medMLG2"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[10,] = c(paste(saida$Media[[1]]), saida$Media[[2]], 
                             saida$Media[[3]],saida$Media[[4]], nome.medida)
medida  = medPOLIORT1
nome.medida = "medPOLIORT1"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[11,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3), 
                              round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)
medida  = medPOLIORT2
nome.medida = "medPOLIORT2"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[12,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3), 
                              round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)
medida  = medPOLIORT3
nome.medida = "medPOLIORT3"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[13,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3), 
                              round(saida$Media[[3]],3),round(saida$Media[[4]],3), nome.medida)

tabela.de.desempenho = tabela.de.desempenho[order(tabela.de.desempenho$Uc), ]
tabela.de.desempenho[nrow(tabela.de.desempenho),1]
medida = tabela.de.desempenho[nrow(tabela.de.desempenho),5]
tabela.de.desempenho
######################################################################################
Medidas = function(residuos,pred,teste,formula){
  medianas = median(pred)
  maximos = max(pred)
  minimos = min(pred)
  Q3 = quantile(pred, probs= 0.75)[[1]]
  Q1 = quantile(pred, probs= 0.25)[[1]]
  interquartil= Q3-Q1
  return(list("Medianas"= medianas, "IQ"=interquartil, "Max"=maximos, "Min"= minimos, "Q3" = Q3,"Q1" = Q3 ))
}
J.k.f = function(dados,cls, formula,...){
  vcs = data.frame()
  n = nrow(dados)
  for (j in 1:n){
    reamostra.jack = dados[-j,]
    treino = reamostra.jack
    teste  = reamostra.jack
    kcls = cls(treino, teste, formula) 
    vcs = rbind(vcs, kcls)
  }
  vcs = vcs[, c(ncol(vcs), 1:(ncol(vcs) - 1))]
  avg = aggregate(. ~ Classificador, data = vcs, FUN = mean)
  bvg = aggregate(. ~ Classificador, data = vcs, FUN = median)
  std = aggregate(. ~ Classificador, data = vcs, FUN = sd)
  return(list("Média" = avg, "Desvio" = std, Mediana= bvg, saidas = vcs))
}

##
LIMITES.DE.CONTROLE = J.k.f(dados,get(medida),formula)
LC =  LIMITES.DE.CONTROLE$Mediana[[2]]
IQ =  LIMITES.DE.CONTROLE$Mediana[[3]]
LSC = LIMITES.DE.CONTROLE$Mediana[[4]]
LIC = LIMITES.DE.CONTROLE$Mediana[[5]]
Q3 =  LIMITES.DE.CONTROLE$Mediana[[6]]
Q1 =  LIMITES.DE.CONTROLE$Mediana[[7]]
Gamma1 = (Q1 - LIC)/IQ
Gamma2 = (LSC - Q3)/IQ

LC1 =  LIMITES.DE.CONTROLE$Média[[2]]
IQ1=  LIMITES.DE.CONTROLE$Média[[3]]
LSC1 = LIMITES.DE.CONTROLE$Média[[4]]
LIC1 = LIMITES.DE.CONTROLE$Média[[5]]
Q31 =  LIMITES.DE.CONTROLE$Média[[6]]
Q11 =  LIMITES.DE.CONTROLE$Média[[7]]
Gamma11 = (Q11 - LIC1)/IQ1
Gamma21 = (LSC1 - Q31)/IQ1


medianas.preditas = LIMITES.DE.CONTROLE$saidas[[2]]
mediana.REFERENCIA = median(dados$y)
mediana.REFERENCIA
paste(LIC,LC,LSC)
paste(LIC1,LC1,LSC1)
paste(min(dados$y), median(dados$y), max(dados$y))
tabela.de.desempenho




## 2.1 

cls = get(medida) # Modelo Vencedor

Holdout = function(dados, p = 0.7){
  n_treino = ceiling(dim(dados)[1]*p)
  ind = c(rep('treino', n_treino), rep('teste', dim(dados)[1] - n_treino) )
  ind = sample(ind)
  treino = dados[ind == 'treino', ]
  teste  = dados[ind == 'teste' , ]
  return(list(treino = treino, teste = teste))
}
treino = Holdout(dados)$treino
teste = Holdout(dados)$teste

Medidas = function(pred){
  predicoes = pred
  return(predicoes)
}
AuxSaida = function(pred, teste, formula, alg){
  ind = which(colnames(teste) == as.character(formula)[2])
  rs = Medidas(pred)
  return(rs)
}

REGULA= function(treino, teste, formula){
  pred = cls(treino,teste,formula)
  MEDIANA= median(pred)
  contador = 1
  while(any(pred< LIC) | any(pred> LSC)){
    Holdout = function(dados, p = 0.7){
      n_treino = ceiling(dim(dados)[1]*p)
      ind = c(rep('treino', n_treino), rep('teste', dim(dados)[1] - n_treino) )
      ind = sample(ind)
      treino = dados[ind == 'treino', ]
      teste  = dados[ind == 'teste' , ]
      return(list(treino = treino, teste = teste))
    }
    treino = Holdout(dados)$treino
    teste = Holdout(dados)$teste
    pred = cls(treino,teste,formula)
    MEDIANA= median(pred)
    contador = contador +1
    print(contador)
  }
  return(list(modelo = cls, mediana = MEDIANA, contador = contador))
}

saida=REGULA(treino, teste, formula)
mediana = saida$mediana
cls= saida$modelo

## 2.2 

## 2.2 

TUNAGEM = function(cls, alpha,dados,Gamma1,Gamma2, Q3,Q1,IQ){
  rep = 10
  contador = 1
  Gamma1.corrigido = Gamma1
  Gamma2.corrigido = Gamma2
  ARL0s = data.frame()
  Gammas1 = data.frame()
  Gammas2 = data.frame()
  NOMINAL= 0
  n = nrow(dados)
  MATRIZ.DE.ERROS = data.frame()
  reamostra = data.frame(matrix(NA,nrow = n*0.7,ncol = length(dados)))
  for (j in 1:rep) {
    NOMINAL= 0
    Gamma1.corrigido = Gamma1
    Gamma2.corrigido = Gamma2
    while (NOMINAL>(1/(alpha/2)) | NOMINAL <( 1/alpha)) { 
      Holdout.reamostra = function(dados, p = 0.7){
        n_treino = ceiling(dim(dados)[1]*p)
        ind = c(rep('reamostra', n_treino), rep('descartada', dim(dados)[1] - n_treino) )
        ind = sample(ind)
        reamostra = dados[ind == 'reamostra', ]
        descarte  = dados[ind == 'descartada' , ]
        return(reamostra = reamostra)
      }
      reamostra = Holdout.reamostra(dados)
      pred = cls(reamostra,reamostra,formula)
      erros = ifelse(pred > (Q3 + Gamma2.corrigido*IQ) |
                       pred< (Q1 - Gamma1.corrigido*IQ),1,0)
      prop.dos.erros = sum(erros)/n
      MATRIZ.DE.ERROS = rbind(erros,MATRIZ.DE.ERROS)
      #################################################################
      NOMINAL= 1/prop.dos.erros
      Gamma1.corrigido= ifelse((NOMINAL > (1/(alpha/2)) & all(pred >= (Q1 - Gamma1.corrigido*IQ))), Gamma1.corrigido - 0.001,
                               ifelse((NOMINAL < (1/(alpha)) & all(pred >= (Q1 - Gamma1.corrigido*IQ))), Gamma1.corrigido,
                                      ifelse((NOMINAL < (1/(alpha)) & (sum(pred >= (Q1 - Gamma1.corrigido*IQ))>=2)), Gamma1.corrigido +0.001,
                                             Gamma1.corrigido)))
      
      Gamma2.corrigido= ifelse((NOMINAL > (1/(alpha/2)) & all(pred < (Q3 + Gamma2.corrigido*IQ))), Gamma2.corrigido - 0.001,
                               ifelse((NOMINAL < (1/(alpha)) & all(pred >= (Q3 + Gamma2.corrigido*IQ))), Gamma2.corrigido,
                                      ifelse((NOMINAL < (1/(alpha)) & (sum(pred >= (Q3 + Gamma2.corrigido*IQ))>=2)), Gamma2.corrigido +0.001,
                                             Gamma2.corrigido)))
      #################################################################
      contador= contador+1
      ARL0s= rbind(NOMINAL,ARL0s)
      A= paste(contador,"-ésima iteração", "com Gamma1 =", 
               Gamma1.corrigido,"com Gamma2 =", Gamma2.corrigido,
               "e ARL_ZERO = ", NOMINAL)
      print(A)
    }
    Gammas1 = rbind(Gamma1.corrigido,Gammas1)
    Gammas2 = rbind(Gamma2.corrigido,Gammas2)
  }
  return(list("ARL ZERO" = ARL0s, "Parâmetros de Tunagem Gammas" =  
                c(Gamma1 = Gamma1.corrigido,Gamma2 = Gamma2.corrigido), 
              "Número de iterações"= contador,
              "Gammas1"=Gammas1, "Gammas2"=Gammas2, "Matriz de Erros"= MATRIZ.DE.ERROS))
}

    
   
simulacao = TUNAGEM(cls, alpha= 0.05 ,dados,Gamma1,Gamma2, Q3,Q1,IQ) # Neste exemplo só foi possível rodar para alpha 0.05
par(mfrow=c(1,2))
boxplot(simulacao$Gammas1, main = "Gamma1's Calibrados", horizontal = F, col = "lightblue")
boxplot(simulacao$Gammas2, main = "Gamma2's Calibrados", horizontal = F, col = "tomato")

Gamma1.c = median(simulacao$Gammas1[[1]])
Gamma2.c = median(simulacao$Gammas2[[1]])

LSCc= Q3 + Gamma2.c*IQ
LICc= Q1 - Gamma1.c*IQ
paste(LICc,LC,LSCc)
paste(LIC,LC,LSC)
paste(min(dados$y),median(dados$y),max(dados$y))
Gamma1.c
Gamma2.c
######### Fase III


Desempenho = function(rep=500, cls, LICc,LSCc){
  ARL = data.frame()
  p= data.frame()
  for (j in 1:rep) {
    new.dados = gera.dados(n=Tamanho.amostras, tipo.y = Densidade.y,
                           tipo.x = Volume.dos.dados, relacao= Presenca.de.Multicolinearidade)
    new.pred = cls(new.dados,new.dados,formula)
  erros = ifelse(new.pred> LSCc | new.pred< LICc,1,0)
  prop.dos.erros = sum(erros)/nrow(new.dados)
  ARL_zero = 1/prop.dos.erros
  #################################################################
  ARL = rbind(ARL_zero,ARL)
  p = rbind(p,erros)
  print(j)
    }
  return(list("ARL's"= ARL,"Erros" = p))
}  

 
des = Desempenho(rep=500, cls, LICc,LSCc)

median(des$`ARL's`)
desempenhos = des$`ARL's`[which(des$`ARL's`[[1]] !=Inf),]
box = boxplot(desempenhos, horizontal = TRUE, main = "ARL0", col = 'lightblue')
abline(v= mean(desempenhos), col = "red", lwd = 2)
median(desempenhos)
mean(desempenhos)
hist(apply(des$Erros, 1, sum), main = "Histograma dos Erros", 
     ylab = "Frequência nas 500 repetições", xlab= "Soma dos Erros das n Amostras Unitárias", 
     col = "tomato")

plot(dados$y, ylim = c(min(dados$y,LICc-2) ,max(dados$y,LSCc+2)),pch=20, xlim= c(1,length(dados$y)+1), cex.main=1.2,
  xlab= "Unidades Amostrais", ylab=expression(paste("Valores Preditos"," ", hat(y))), col=ifelse(dados$y>LSCc | dados$y<LICc,"red","black"),
     main="Gráfico de Controle dos Rankings Preditos",  col.main="blue",col.lab="blue",
     font.main=2, font.lab=2, cex.lab=1, cex.sub=1)
lines(c(rep(LSCc,length(dados$y))),col ="red", lwd= 2)
lines(c(rep(LICc,length(dados$y))),col ="red", lwd= 2)
lines(c(rep(LC,length(dados$y))),col ="blue", lwd= 2)
segments(c(1:length(dados$y)),LC, c(1:length(dados$y)),dados$y ,lty= 2)
text(x=1,y=LSCc+1.5,"LSC", col = "red", cex=0.8)
text(x=1,y=LICc-1.5,"LIC", col = "red", cex=0.8)


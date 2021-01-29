
#PARÂMETROS DA SIMULAÇÃO
Tamanho.amostras = 250
Densidade.y =  "simetrico" # "Assimetrico"
Volume.dos.dados = "Nao.Big" # "Big" "
Presenca.de.Multicolinearidade = "X.indep" # X.multicol
q = c(0,0.1,0.4,0.6,1) # Criterio de Discretização


gera.dados = function(n, tipo.y, tipo.x, relacao){
library(gtools)
  
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
  refe = ifelse(tipo.x=="Big", data.frame(1:15), data.frame(1:5))
  refe= refe[[1]]
  X= X[,refe]
  dados = X
  dados$y = y
  dados$y =  quantcut(dados$y,q)
  dados$y = (as.numeric(dados$y))-1
  dados$y = as.factor(dados$y)
  return(dados)
}
dados = gera.dados(n=Tamanho.amostras, tipo.y = Densidade.y,
                   tipo.x = Volume.dos.dados, relacao= Presenca.de.Multicolinearidade)


################################################
varResp = 'y'
Medidas = function(mc){
  PE = (sum(mc[1,])*sum(mc[,1])+sum(mc[2,])*sum(mc[,2])+
          sum(mc[3,])*sum(mc[,3])+sum(mc[4,])*sum(mc[,4]))/((sum(mc[1,])+sum(mc[2,])+sum(mc[3,])+sum(mc[4,]))^2)
  PA = sum(diag(mc))/(sum(mc[1,])+sum(mc[2,])+sum(mc[3,])+sum(mc[4,]))
  Kappa = (PA-PE)/(1-PE)
  return(list(Kappa = Kappa))
}
formula = as.formula(paste(varResp, '~ .')) 

FormFac = function(formula){
  aux = as.character(formula)[2]
  return(as.formula(paste('factor(',aux, ')', '~ .')))
}
AuxSaida = function(pred, teste, formula, alg){
  ind = which(colnames(teste) == as.character(formula)[2])
  mc = matrix(NA,4,4)
  mc[1,1] = sum(teste[, ind]==0 & pred ==0)
  mc[1,2] = sum(teste[, ind]==1 & pred ==0)
  mc[1,3] = sum(teste[, ind]==2 & pred ==0)
  mc[1,4] = sum(teste[, ind]==3 & pred ==0)
  mc[2,1] = sum(teste[, ind]==0 & pred ==1)
  mc[2,2] = sum(teste[, ind]==1 & pred ==1)
  mc[2,3] = sum(teste[, ind]==2 & pred ==1)
  mc[2,4] = sum(teste[, ind]==3 & pred ==1)
  mc[3,1] = sum(teste[, ind]==0 & pred ==2)
  mc[3,2] = sum(teste[, ind]==1 & pred ==2)
  mc[3,3] = sum(teste[, ind]==2 & pred ==2)
  mc[3,4] = sum(teste[, ind]==3 & pred ==2)
  mc[4,1] = sum(teste[, ind]==0 & pred ==3)
  mc[4,2] = sum(teste[, ind]==1 & pred ==3)
  mc[4,3] = sum(teste[, ind]==2 & pred ==3)
  mc[4,4] = sum(teste[, ind]==3 & pred ==3)
  rs = Medidas(mc)
  rs[['Classificador']] = alg    
  return(do.call(data.frame, rs))  
}
Holdout = function(dados, p = 0.7){
  n_treino = ceiling(dim(dados)[1]*p)
  ind = c(rep('treino', n_treino), rep('teste', dim(dados)[1] - n_treino) )
  ind = sample(ind)
  treino = dados[ind == 'treino', ]
  teste  = dados[ind == 'teste' , ]
  return(list(treino = treino, teste = teste))
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
  avg = aggregate(. ~ Classificador, data = x, FUN = mean)### valor guardado que entra no gráfico da figura 5.3
  std = aggregate(. ~ Classificador, data = x, FUN = sd)
  return(list("Media Global"= avg, "Erro Padrão Validação Cruzada" = std, Modelos = x))
}
###### FASE I   
############################################################# SVM
library(kernlab)
SVM1 = function(treino, teste, formula){
  cls = ksvm(formula, data = data.frame(treino),type='C-svc',kernel= 'rbfdot')
  pred = predict(cls, newdata = data.frame(teste))
  return(pred)
}
medSVM1 = function(treino, teste, formula){
  alg = 'SVM - Kernel Gaussiano' 
  pred = SVM1(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}

SVM2 = function(treino, teste, formula){
  cls = ksvm(formula, data = data.frame(treino),type='C-svc',kernel= 'polydot')
  pred = predict(cls, newdata = data.frame(teste))
  return(pred)
}
medSVM2 = function(treino, teste, formula){
  alg = 'SVM - Kernel Polinomial' 
  pred = SVM2(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}

SVM3 = function(treino, teste, formula){
  cls = ksvm(formula, data = data.frame(treino),type='C-svc',kernel= 'vanilladot')
  pred = predict(cls, newdata = data.frame(teste))
  return(pred)
}
medSVM3 = function(treino, teste, formula){
  alg = 'SVM - Kernel Linear' 
  pred = SVM3(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}
####################################################ARVORE DE DECISÃO
Arvore = function(treino, teste, formula){
  library(rpart) 
  cls = rpart(formula, data = data.frame(treino), method = 'class')
  pred = predict(cls, newdata = data.frame(teste), type = 'class')
  return(pred)
}

medArvore = function(treino, teste, formula){ 
  alg = 'Árvore de Decisão' 
  pred = Arvore(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}
###########################################################KNN
KNN = function(treino, teste, formula, K = (nrow(treino)/3)){
  library(kknn) 
  cls = kknn(FormFac(formula), train = data.frame(treino), test = data.frame(teste), k = K) 
  pred = cls$fitted.values 
  return(pred)
}

medKNN = function(treino, teste, formula, K = (nrow(treino)/3)){ 
  alg = paste(K, '- Vizinhos mais Próximos')
  pred = KNN(treino, teste, formula, K)
  return(AuxSaida(pred, teste, formula, alg))
}
###########################################################LDA
LDA = function(treino, teste, formula){
  library(MASS) 
  cls = lda(formula, data = data.frame(treino))
  pred = predict(cls, newdata = data.frame(teste))$class 
  return(pred)
}

medLDA = function(treino, teste, formula){ 
  alg = 'Análise Discriminante Linear' 
  pred = LDA(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}
###########################################################QDA
QDA = function(treino, teste, formula){
  library(MASS) 
  cls = qda(formula, data = data.frame(treino))
  pred = predict(cls, newdata = data.frame(teste))$class 
  return(pred)
}

medQDA = function(treino, teste, formula){
  alg = 'Análise Discriminante Quadrática' 
  pred = QDA(treino, teste, formula) 
  return(AuxSaida(pred, teste, formula, alg))
}
#######################################################REDE NEURAL 
library(nnet)
ANN = function(treino, teste, formula, size = 3){
  classifier = nnet(formula, data = data.frame(treino), size = 3)
  pred  = predict(classifier, newdata =data.frame(teste))
  pred = apply(pred,1,which.max)-1
  pred= as.vector(pred)
  pred= as.factor(pred)
  return(pred)
}
medANN = function(treino, teste, formula, size =3){
  alg = 'Redes Neurais' 
  pred = ANN(treino, teste, formula, size = 3 )
  return(AuxSaida(pred, teste, formula, alg))
}
#######################################################MULTINOMIAL
library(nnet)
Multinomial = function(treino, teste, formula){ 
  cls =multinom(formula, data = data.frame(treino))
  pred = predict(cls, newdata = teste, type = 'class') 
  return(pred)
}
medMultinomial = function(treino, teste, formula){
  alg = 'Regressão Multinomial' 
  pred = Multinomial(treino, teste, formula)
  return(AuxSaida(pred, teste, formula, alg))
}
##################################################FLORESTAS ALEATÓRIAS 
Floresta = function(treino, teste, formula){
  library(randomForest) 
  cls = randomForest(as.factor(y)~., data = data.frame(treino),importance= TRUE,
                     proximity=TRUE)
  pred = predict(cls, newdata = data.frame(teste), type = "class")
  return(pred)
}
medFloresta = function(treino, teste, formula){ 
  alg = 'Floresta Aleatória' 
  pred = Floresta(data.frame(treino), data.frame(teste), formula)
  return(AuxSaida(pred, teste, formula, alg))
}
##############################################################


tabela.de.desempenho = matrix(NA, ncol = 2, nrow = 10)
colnames(tabela.de.desempenho) = c("Modelos", "kappa")
tabela.de.desempenho = data.frame(tabela.de.desempenho)


medida  = medSVM1
nome.medida = "medSVM1"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[1,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3))
medida  = medSVM2
nome.medida = "medSVM2"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[2,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3))
medida  = medSVM3
nome.medida = "medSVM3"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[3,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3))
medida  = medArvore
nome.medida = "medArvore"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[4,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3))
medida  = medKNN
nome.medida = "medKNN"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[5,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3))
medida  = medLDA
nome.medida = "medLDA"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[6,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3))
medida  = medQDA
nome.medida = "medQDA"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[7,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3))
medida  = medANN
nome.medida = "medANN"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[8,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3))
medida  = medMultinomial
nome.medida = "medMultinomial"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[9,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3))
medida  = medFloresta
nome.medida = "medFloresta"
saida = ValidCruzadarep(medida, dados, formula, kfolds = 10, reps = 10)
tabela.de.desempenho[10,] = c(paste(saida$Media[[1]]), round(saida$Media[[2]],3))

tabela.de.desempenho = tabela.de.desempenho[order(tabela.de.desempenho$kappa), ]
medida = tabela.de.desempenho[nrow(tabela.de.desempenho),1]
medida = medMultinomial
tabela.de.desempenho
######################################################################################
Medidas = function(pred){
  proporcoes = table(pred)/sum(table(pred))
  prop.0 = proporcoes[[1]]
  prop.1 = proporcoes[[2]]
  prop.2 = proporcoes[[3]]
  prop.3 = proporcoes[[4]]
  
  return(list("Prop 0"= prop.0, "Prop 1"=prop.1, "Prop 2"=prop.2, "Prop 3"=prop.3))
}
AuxSaida = function(pred, teste, formula, alg){
  rs = Medidas(pred)
  rs[['Classificador']] = alg    
  return(do.call(data.frame, rs))  
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



##LIMITES.DE.CONTROLE = J.k.f(dados,get(medida),formula)

LIMITES.DE.CONTROLE = J.k.f(dados,medida,formula)
LC.0 =  LIMITES.DE.CONTROLE$Média[[2]]
LC.1 =  LIMITES.DE.CONTROLE$Média[[3]]
LC.2 = LIMITES.DE.CONTROLE$Média[[4]]
LC.3 = LIMITES.DE.CONTROLE$Média[[5]]

s.0 =  LIMITES.DE.CONTROLE$Desvio[[2]]
s.1 =  LIMITES.DE.CONTROLE$Desvio[[3]]
s.2 = LIMITES.DE.CONTROLE$Desvio[[4]]
s.3 = LIMITES.DE.CONTROLE$Desvio[[5]]

LSC.0 = LC.0 + 3*s.0
LSC.1 = LC.1 + 3*s.1
LSC.2 = LC.2 + 3*s.2
LSC.3 = LC.3 + 3*s.3
LIC.0 = LC.0 - 3*s.0
LIC.1 = LC.1 - 3*s.1
LIC.2 = LC.2 - 3*s.2
LIC.3 = LC.3 - 3*s.3

lim.0 = paste0('[',round(LIC.0,3) ,';',round(LC.0,3) ,';',round(LSC.0,3) ,']')
lim.1 = paste0('[',round(LIC.1,3) ,';',round(LC.1,3) ,';',round(LSC.1,3) ,']')
lim.2 = paste0('[',round(LIC.2,3) ,';',round(LC.2,3) ,';',round(LSC.2,3) ,']')
lim.3 = paste0('[',round(LIC.3,3) ,';',round(LC.3,3) ,';',round(LSC.3,3) ,']')

verificacao.dos.intervalos = table(dados$y)/sum(table(dados$y))
lim.0;lim.1;lim.2;lim.3 ;verificacao.dos.intervalos
## 2.1 

# cls = get(medida) # Modelo Vencedor
cls = medida

##########
TUNAGEM = function(cls, alpha,dados){
  rep = 10
  contador = 1
  Gamma = 1
  Gammas = data.frame()
  ARL0s = data.frame()
  NOMINAL= 0
  n = nrow(dados)
  MATRIZ.DE.ERROS = data.frame()
  reamostra = data.frame(matrix(NA,nrow = n*0.7,ncol = length(dados)))
  for (j in 1:rep) {
    NOMINAL= 0
    Gamma.corrigido = Gamma
    while (NOMINAL>(1/(alpha/2)) | NOMINAL <( 1/alpha)) { 
      erros= c(rep(NA,nrow(dados)))
      for (k in 1:nrow(dados)) {
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
      pred$Prop.0
      pred$Prop.1
      pred$Prop.2
      pred$Prop.3
      erros[k] = ifelse(LIC.0 > pred$Prop.0 | pred$Prop.0 > LSC.0  |
                        LIC.1 > pred$Prop.1 | pred$Prop.1 > LSC.1  | 
                         LIC.2 > pred$Prop.2 | pred$Prop.2 > LSC.2  |
                          LIC.3 > pred$Prop.3 | pred$Prop.3 > LSC.3, 1,0)
      
      }
      prop.dos.erros = sum(erros)/n
      MATRIZ.DE.ERROS = rbind(erros,MATRIZ.DE.ERROS)
################################ PAREI AQUI #################################
      NOMINAL= 1/prop.dos.erros
      Gamma.corrigido= ifelse((NOMINAL > (1/(alpha/2))),Gamma.corrigido- 0.5,  
       ifelse((NOMINAL < (1/(alpha))),Gamma.corrigido+ 0.5,Gamma.corrigido))
      #################################################################

      LSC.0 = LC.0 + 3*Gamma.corrigido*s.0
      LSC.1 = LC.1 + 3*Gamma.corrigido*s.1
      LSC.2 = LC.2 + 3*Gamma.corrigido*s.2
      LSC.3 = LC.3 + 3*Gamma.corrigido*s.3
      LIC.0 = LC.0 - 3*Gamma.corrigido*s.0
      LIC.1 = LC.1 - 3*Gamma.corrigido*s.1
      LIC.2 = LC.2 - 3*Gamma.corrigido*s.2
      LIC.3 = LC.3 - 3*Gamma.corrigido*s.3
      
      
      contador= contador+1
      ARL0s= rbind(NOMINAL,ARL0s)
      A= paste(contador,"-ésima iteração", "com Gamma =", 
               Gamma.corrigido, NOMINAL)
      print(A)
    }
    Gammas = rbind(Gamma.corrigido,Gammas)
  }
  return(list("ARL ZERO" = ARL0s, "Parâmetros de Tunagem Gammas" =  
               Gammas, 
              "Número de iterações"= contador,
              "Gammas"=Gammas, "Matriz de Erros"= MATRIZ.DE.ERROS))
}

simulacao = TUNAGEM(cls, alpha=0.05,dados) # Neste exemplo só foi possível rodar para alpha 0.05
par(mfrow=c(2,1))
ref = bquote("Gamma's Calibrados para" ~  alpha == 0.05)
boxplot(simulacao$`Parâmetros de Tunagem Gammas`, main = ref, horizontal = T, col = "lightblue")
points(mean(simulacao$`Parâmetros de Tunagem Gammas`),1,'*',col = 'red')
simulacao$`Número de iterações`
Gamma.c = median(simulacao$Gammas[[1]])

LSC.0 = LC.0 + 3*Gamma.c*s.0
LSC.1 = LC.1 + 3*Gamma.c*s.1
LSC.2 = LC.2 + 3*Gamma.c*s.2
LSC.3 = LC.3 + 3*Gamma.c*s.3
LIC.0 = LC.0 - 3*Gamma.c*s.0
LIC.1 = LC.1 - 3*Gamma.c*s.1
LIC.2 = LC.2 - 3*Gamma.c*s.2
LIC.3 = LC.3 - 3*Gamma.c*s.3

lim.0.c = paste0('[',round(LIC.0,3) ,';',round(LC.0,3) ,';',round(LSC.0,3) ,']')
lim.1.c = paste0('[',round(LIC.1,3) ,';',round(LC.1,3) ,';',round(LSC.1,3) ,']')
lim.2.c = paste0('[',round(LIC.2,3) ,';',round(LC.2,3) ,';',round(LSC.2,3) ,']')
lim.3.c = paste0('[',round(LIC.3,3) ,';',round(LC.3,3) ,';',round(LSC.3,3) ,']')

lim.0; lim.1;lim.2;lim.3; lim.0.c; lim.1.c;lim.2.c;lim.3.c

verificacao.dos.intervalos
######### Fase III



Desempenho = function(rep=500, cls){
  ARL = data.frame()
  p= data.frame()
  for (j in 1:rep) {
    for (k in 1:nrow(dados)) {
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
    pred$Prop.0
    pred$Prop.1
    pred$Prop.2
    pred$Prop.3
    erros[k] = ifelse(LIC.0 > pred$Prop.0 | pred$Prop.0 > LSC.0  |
                        LIC.1 > pred$Prop.1 | pred$Prop.1 > LSC.1  | 
                        LIC.2 > pred$Prop.2 | pred$Prop.2 > LSC.2  |
                        LIC.3 > pred$Prop.3 | pred$Prop.3 > LSC.3, 1,0)
    
  }
  prop.dos.erros = sum(erros)/nrow(dados)
  ARL_zero = 1/prop.dos.erros
  #################################################################
  ARL = rbind(ARL_zero,ARL)
  p = rbind(p,erros)
  print(j)
    }
  return(list("ARL's"= ARL,"Erros" = p))
}  

 
des = Desempenho(rep=500, cls)



median(des$`ARL's`,na.rm = FALSE)
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



a = sample(1:8,size = 600,replace = TRUE)
boxplot(a, main = ref, horizontal = T, col = "lightblue")
points(mean(a),1,pch= '*',col = 'red')




library(data.table)
library(Hmisc)
library(pls)
library(MASS)

path2data<-file.path("c:","Users","makch","STUDIES",
                     "s3","ADD","RegPLSR")


DM.dt <- fread(file.path(path2data, "Data_Cornell.csv"))
# Regression multiple Y f(X1..X7)
lm<-lm(formula = Y ~ 0 + X1 + X2+ X3 + X4 + X5 +X6 + X7, 
       data=DM.dt)
#print 
print(lm) 
#summary
print(summary(lm)) 

DM_Matrix<-as.matrix(DM.dt)
corrmat = rcorr(DM_Matrix, type=c("pearson","spearman"))$P
n_features = ncol(DM.dt)
corrmat[n_features,1:n_features]
Scale_DM<-scale(DM.dt)

colnames(Scale_DM)<- c("X1_cn","X2_cn","X3_cn","X4_cn","X5_cn",
                       "X6_cn","X7_cn","Y_cn")
DT_scale<- cbind(DM.dt,Scale_DM)
#################Construction de T1 ###########################
DT_scale<-DT_scale[,':='(T1=(1/sqrt(0.84^2+0.84^2+0.71^2
                                    +0.99^2+0.74^2))
                         *((-0.84*X1_cn)+(-0.84* X3_cn)+(-0.71*X4_cn)
                          +(0.99* X6_cn)+(-0.74* X7_cn)))]
#################Construction de T2 ###########################

# Regression y sur T1 et Xj j=1..7
#pour chercher les variables contribuant de mani�re sifnificatice
# � la construction de T2

lm11<-lm(formula = Y_cn ~ 0 + T1 + X1_cn, data=DT_scale)
print(summary(lm11))
lm12<-lm(formula = Y_cn ~ 0 + T1 + X2_cn, data=DT_scale)
print(summary(lm12))
lm13<-lm(formula = Y_cn ~ 0 + T1 + X3_cn, data=DT_scale)
print(summary(lm13))
lm14<-lm(formula = Y_cn ~ 0 + T1 + X4_cn, data=DT_scale)
print(summary(lm14))
lm15<-lm(formula = Y_cn ~ 0 + T1 + X5_cn, data=DT_scale)
print(summary(lm15))
lm16<-lm(formula = Y_cn ~ 0 + T1 + X6_cn, data=DT_scale)
print(summary(lm16))
lm17<-lm(formula = Y_cn ~ 0 + T1 + X7_cn, data=DT_scale)
print(summary(lm17))

# Seules les variables X2 et X6 sont significatives au risque 0.05

# On calcule les r�sidus X12 X16 des regressions de X2_nc, X6_nc sur T1 
lm_R12<-lm(formula = X2_cn ~ 0 + T1  , data=DT_scale)
#### Extraction des r�sidus
X12<-resid(lm_R12)
lm_R16<-lm(formula = X6_cn ~ 0 + T1  , data=DT_scale)
#### Extraction des r�sidus
X16<-resid(lm_R16)
DT_scale<- cbind(DT_scale,X12,X16)
X12n<-X12/var(X12)
X16n<-X16/var(X16)
DT_scale<- cbind(DT_scale,X12n,X16n)

# Puis on effectue les deux r�gressions multiples 
# Y_cn sur T1 et X1jn = x1j/var(x1j) j=2 et 6
lm_Y12<-lm(formula = Y_cn ~ 0 + T1 + X12n  , data=DT_scale)
lm_Y16<-lm(formula = Y_cn ~ 0 + T1 + X16n  , data=DT_scale)
print(summary(lm_Y12))
print(summary(lm_Y16))

DT_scale<-DT_scale[,':='(T2=(1/sqrt(0.19361^2+0.10067^2))*
                           ((-0.19361*X12)+(0.10067*X16)))]

#################Construction de T3 ###########################

# Regressions y sur T1, T2 et Xj j=1..7
#pour chercher les variables contribuant de mani�re sifnificatice
# � la construction de T3

lm21<-lm(formula = Y_cn ~ 0 + T1 + T2 + X1_cn, data=DT_scale)
print(summary(lm21))
lm22<-lm(formula = Y_cn ~ 0 + T1 + T2 + X2_cn, data=DT_scale)
print(summary(lm22))
lm23<-lm(formula = Y_cn ~ 0 + T1 + T2 + X3_cn, data=DT_scale)
print(summary(lm23))
lm24<-lm(formula = Y_cn ~ 0 + T1 + T2 + X4_cn, data=DT_scale)
print(summary(lm24))
lm25<-lm(formula = Y_cn ~ 0 + T1 + T2 + X5_cn, data=DT_scale)
print(summary(lm25))
lm26<-lm(formula = Y_cn ~ 0 + T1 + T2 + X6_cn, data=DT_scale)
print(summary(lm26))
lm27<-lm(formula = Y_cn ~ 0 + T1 + T2 + X7_cn, data=DT_scale)
print(summary(lm27))

# les variables X1 X2 X3 X4 et X6 sont significatives au risque 0.05
# Il faut rechercher une troisi�me composante PLS T3

# On calcule les r�sidus X21 X22 X23 X24 et X26 
# des regressions de X1_nc, X2_nc, X3_nc, X4_nc et X6_nc sur T1 T2 
lm_R21<-lm(formula = X1_cn ~ 0 + T1 + T2  , data=DT_scale)
#### Extraction des r�sidus
X21<-resid(lm_R21)
lm_R22<-lm(formula = X2_cn ~ 0 + T1 + T2  , data=DT_scale)
#### Extraction des r�sidus
X22<-resid(lm_R22)
lm_R23<-lm(formula = X3_cn ~ 0 + T1 + T2  , data=DT_scale)
#### Extraction des r�sidus
X23<-resid(lm_R23)
lm_R24<-lm(formula = X4_cn ~ 0 + T1 + T2  , data=DT_scale)
#### Extraction des r�sidus
X24<-resid(lm_R24)
lm_R26<-lm(formula = X6_cn ~ 0 + T1 + T2 , data=DT_scale)
#### Extraction des r�sidus
X26<-resid(lm_R26)
DT_scale<- cbind(DT_scale,X21,X22,X23,X24,X26)
X21n<-X21/var(X21)
X22n<-X22/var(X22)
X23n<-X23/var(X23)
X24n<-X24/var(X24)
X26n<-X26/var(X26)
DT_scale<- cbind(DT_scale,X21n,X22n,X23n,X24n,X26n)

# Puis on effectue les  r�gressions multiples 
# Y_cn sur T1 T2 et X1jn = x1j/var(x1j) j=1,2,3,4 et 6
lm_Y21<-lm(formula = Y_cn ~ 0 + T1 + T2 + X21n  , data=DT_scale)
lm_Y22<-lm(formula = Y_cn ~ 0 + T1 + T2 + X22n  , data=DT_scale)
lm_Y23<-lm(formula = Y_cn ~ 0 + T1 + T2 + X23n  , data=DT_scale)
lm_Y24<-lm(formula = Y_cn ~ 0 + T1 + T2 + X24n  , data=DT_scale)
lm_Y26<-lm(formula = Y_cn ~ 0 + T1 + T2 + X26n  , data=DT_scale)
print(summary(lm_Y21))
print(summary(lm_Y22))
print(summary(lm_Y23))
print(summary(lm_Y24))
print(summary(lm_Y26))
DT_scale<-DT_scale[,':='(T3=(1/sqrt(0.05830^2+0.015440^2+0.05864^2
                                    +0.07722^2+0.029694^2))
                         *((0.05830*X21)+(0.015440*X22)+(0.05864*X23)
                           +(-0.07722*X24)+(0.029694*X26)))]

#################Construction de T4!! ###########################

# Regressions y sur T1, T2, T3 et Xj j=1..7
#pour chercher les variables contribuant de mani�re sifnificatice
# � la construction de T3

lm31<-lm(formula = Y_cn ~ 0 + T1 + T2 + T3 + X1_cn, data=DT_scale)
print(summary(lm31))
lm32<-lm(formula = Y_cn ~ 0 + T1 + T2 + T3 + X2_cn, data=DT_scale)
print(summary(lm32))
lm33<-lm(formula = Y_cn ~ 0 + T1 + T2 + T3 + X3_cn, data=DT_scale)
print(summary(lm33))
lm34<-lm(formula = Y_cn ~ 0 + T1 + T2 + T3 + X4_cn, data=DT_scale)
print(summary(lm34))
lm35<-lm(formula = Y_cn ~ 0 + T1 + T2 + T3 + X5_cn, data=DT_scale)
print(summary(lm35))
lm36<-lm(formula = Y_cn ~ 0 + T1 + T2 + T3 + X6_cn, data=DT_scale)
print(summary(lm36))
lm37<-lm(formula = Y_cn ~ 0 + T1 + T2 + T3 + X7_cn, data=DT_scale)
print(summary(lm37))

# Aucune des variables X1...X7 sont significatives au risque 0.05
# Il faut retenir que les trois composantes PLS T1 T2 T3

##########Construction de l'�quation de r�gression PLS #########
##########           � trois composantes               #########

# Regressions y sur T1, T2 T3 


lm_PLS<-lm(formula = Y ~   T1 + T2 + T3, data=DT_scale)
print(summary(lm_PLS))

# Avec les variables d'origine

Y = 93.317 - 8.755*x1 - 7.782*x2 - 14.969*x3 - 8.434*x4 + 9.488*x6 
- 44.978*x7

Y=84.562*x1 + 85.535*x2 + 78.348*x3 + 84.883*x4 + 93.317*x5 
+ 102.805*x6 + 48.339*x7

#######Regression lineaire generalise PLS#############
####### Logistique##############################

DM_logistic.dt <- fread(file.path(path2data, "QUALITE_RAISIN.csv"))
corr<-DM_logistic.dt[,3:6]

corr_Matrix<-as.matrix(corr)
rcorr(corr_Matrix, type=c("pearson","spearman"))



Scale_DM_logistic<-scale(DM_logistic.dt)
Scale_DM_logistic<-as.data.table(Scale_DM_logistic)
a<-Scale_DM_logistic[,3:6]
colnames(a)<- c("TEMPERATURE_cn","SOLEIL_cn","CHALEUR_cn","PLUIE_cn")
DT_scale_logistic<- cbind(DM_logistic.dt,a)
donnees<-DT_scale_logistic[,7:11]
donnees$QUALITE<-as.factor(donnees$QUALITE)
model <- polr(QUALITE~.,data=donnees,Hess=TRUE )
print(summary(model))
#############La regression logistique ordinale PLS#######################
model_1 <- polr(QUALITE~TEMPERATURE_cn,data=donnees)
print(summary(model_1))
model_2 <- polr(QUALITE~SOLEIL_cn,data=donnees)
print(summary(model_2))
model_3 <- polr(QUALITE~CHALEUR_cn,data=donnees)
print(summary(model_3))
model_4 <- polr(QUALITE~PLUIE_cn,data=donnees)
print(summary(model_4))
#################Construction de T1 ###########################
donnees<-donnees[,':='(T1=(1/sqrt(3.012^2+3.34^2+2.145^2+1.791^2))*((-3.012*TEMPERATURE_cn)
                                                       +(-3.34*SOLEIL_cn)+(-2.145*CHALEUR_cn)+(1.791*PLUIE_cn)))]
#################Construction de T2 !!###########################                                                          
####### On construit les r�gression logistiques de la qualit�
####### sur T1 et chaque pr�dicteur centr� reduit##############
model_11 <- polr(QUALITE~T1+TEMPERATURE_cn,Hess=TRUE ,data=donnees)
print(summary(model_11))
model_21 <- polr(QUALITE~T1+SOLEIL_cn,Hess=TRUE ,data=donnees)
print(summary(model_21))
model_31 <- polr(QUALITE~T1+CHALEUR_cn,Hess=TRUE ,data=donnees)
print(summary(model_31))
model_41 <- polr(QUALITE~T1+PLUIE_cn,Hess=TRUE ,data=donnees)
print(summary(model_41))
# Aucune des variables sont significatives au risque 0.05
# Il faut retenir que la premi�re composante PLS T1
####### Regression Logistique PLS sur T1######################## 
model_LOG_PLS <- polr(QUALITE~T1,Hess=TRUE ,data=donnees)
print(summary(model_LOG_PLS))






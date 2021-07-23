#  IMPLEMENTATION DE LA PROCEDURE AIC 

install.packages("powerSet") # Ensemble de parties
library(rje)
library(powerSet)

# Paramètres :


n <- 20 # nombre d'observations

p <- 20 # nombre de variables explicatives / covariables

d <- 7 # au choix

sigma2 <- 1

f_star <- rep(0, n) # le signal f*

Y <- rnorm(n, 0, sigma2) # Y est la variable à expliquer. Ici, Y est un bruit blanc gaussien (le signal f* vaut exactement zéro)

mod <- powerSet(1:p, d) # mod est de type liste, c'est l'ensemble des parties de [1,p] à 0 à d éléments ; correspond à tous les modèles possibles de taille 0 à d.
head(mod) # aperçu


# Nous trions la liste mod :

m_sorted <- list()

cpt <- 1

for (dim in 0:d) {
  mod_dimension_dim <- which(lengths(mod) == dim)

  for (i in mod_dimension_dim) {
    m_sorted[[cpt]] <- mod[[i]]

    cpt <- cpt + 1
  }
}

mod <- m_sorted

# mod est maintenant "trié" (d'abord les modeles de taille 1, puis les modeles de taille 2, etc).


#_________________________________________________________________________________________________

obs_values <- c() #contient \hat{r}_{m}


for (m in mod) {
  
  dim <- length(m)

  Proj_Sm_Y <- rep(0, n) 
  
  # Projection du vecteur Y sur l'espace engendré par les colonnes I^j de l'identité, j \in m

  Proj_Sm_Y[m] <- Y[m] #\hat{f}_{m}

  v <- Y - Proj_Sm_Y

  norme <- norm(v, type = "2")^2

  obs_values <- c(obs_values, norme + (2 * dim - n) * sigma2)
}

res <- list(obs_values[which(lengths(mod) == 0)], obs_values[which(lengths(mod) == 1)], obs_values[which(lengths(mod) == 2)], obs_values[which(lengths(mod) == 3)], obs_values[which(lengths(mod) == 4)], obs_values[which(lengths(mod) == 5)], obs_values[which(lengths(mod) == 6)], obs_values[which(lengths(mod) == 7)])

boxplot(res, horizontal = FALSE, names = c("dim=0", "dim=1", "dim=2", "dim=3", "dim=4", "dim=5", "dim=6", "dim=7"), col = c("cadetblue2"))

points(0:7, col = 2, pch = 19) #risque réel pour chaque dimension 

# _________________________________________________________________________________________________


# le risque min est
min(obs_values)

# et le modèle sélectionné est :
(modele_AIC <- mod[[which.min(obs_values)]])

# ref modele_AIC: which.min(obs_values)

# L'estimateur \hat{f}_{m_AIC} correspond est : 

f_hat <- rep(0, n)
(f_hat[modele_AIC] <- Y[modele_AIC])  #Projection du vecteur Y sur le modèle sélectionné par la procédure AIC 


(risk <- norm(f_star - f_hat, type = "2")^2) #risque associé 



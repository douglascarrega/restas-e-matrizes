library(emmeans)

#RETAS EX 1
###############################################
coordenadas_das_retas <- function(x1, y1, x2, y2, x3, y3, x4, y4) {
  # coeficientes angulares das retas
  m1 <- (y2 - y1) / (x2 - x1)
  m2 <- (y4 - y3) / (x4 - x3)
  
  # verificando se são diferentes
  if (m1 != m2) {
    # coordenadas do ponto de interseção
    x <- (m1*x1 - m2*x3 + y3 - y1) / (m1 - m2)
    y <- m1*(x - x1) + y1
    
    cat("As retas são concorrentes e se cruzam no ponto (", x, ", ", y, ").\n")
  } else {
    cat("As retas são paralelas.\n")
  }
}

#insira os valores das coordenadas na 
#seguinte sequência(x1, y1, x2, y2, x3, y3, x4, y4)
#exemplo
coordenadas_das_retas(-1, 1, 0, 1, 0, -1, 1, 0)



#######################################################################
#EXERCÍCIO 2 MATRIZ LU

data=c(1,2,1,2,6,6,1,6,3)
print(data)

A=matrix(data, 3, 3)
print(A)
det(A)
inversa_A = solve(A)
inversa_A



# Função para decomposição LU
# A é a matriz do sistema e n é a ordem da matriz
decomposicao_lu <- function(A, n) {
  L <- matrix(0, nrow = n, ncol = n)
  U <- matrix(0, nrow = n, ncol = n)
  
  # Algoritmo de eliminação de Gauss com pivotação parcial
  for (k in 1:n) {
    # Pivotamento parcial
    p <- which(abs(A[k:n, k]) == max(abs(A[k:n, k])))[1] + k - 1
    if (p != k) {
      A[c(k, p), ] <- A[c(p, k), ]
      L[c(k, p), ] <- L[c(p, k), ]
      if (k > 1) {
        U[c(k, p), 1:(k-1)] <- U[c(p, k), 1:(k-1)]
      }
    }
    
    # Fatoração LU
    L[k, k] <- 1
    for (i in (k+1):n) {
      L[i, k] <- A[i, k] / A[k, k]
      U[k, i] <- A[k, i]
      for (j in (k+1):n) {
        A[i, j] <- A[i, j] - L[i, k]*U[k, j]
      }
    }
    U[k, k] <- A[k, k]
  }
  
  # Retorna L e U
  return(list(L = L, U = U))
}

# Função para resolver sistemas de equações lineares por decomposição LU
# A é a matriz do sistema, Y é o vetor coluna do lado direito e n é a ordem da matriz
lu_sistema <- function(A, Y, n) {
  # Decomposição LU
  LU <- decomposicao_lu(A, n)
  L <- LU$L
  U <- LU$U
  
  # Substituição progressiva
  Z <- rep(0, n)
  for (i in 1:n) {
    soma <- 0.0
    for (j in 1:(i-1)) {
      soma <- soma + L[i,j]*Z[j]
    }
    Z[i] <- Y[i] - soma
  }
  
  # Substituição retroativa
  X <- rep(0, n)
  for (i in n:1) {
    soma <- 0.0
    for (j in (i+1):n) {
      soma <- soma + U[i,j]*X[j]
    }
    X[i] <- (Z[i] - soma) / U[i,i]
  }
  
  # Retorna a solução
  return(X)
}

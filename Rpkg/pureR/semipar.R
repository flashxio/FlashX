
semipar = function(A1,A2,dhat)
{
  Xhat1 <- ase(A1,dhat)
  Xhat2 <- ase(A2,dhat)
  procrustes(Xhat1,Xhat2)$error
}


#include <Rcpp.h>

// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>

#include <list>

// [[Rcpp::export]]
Rcpp::DataFrame BigMatBH(Rcpp::XPtr<BigMatrix> pBigMat, double saveMax) {

  double *p = (double *)pBigMat->data_ptr();
  const size_t numCol = pBigMat->ncol();
  const size_t numRow = pBigMat->nrow();
  const size_t n = numCol * numRow;

  std::vector<double> sp;
  
  cout << "Copying values less than " << saveMax << "..." << endl;
  double *cp = p;
  for (size_t i = 0; i < n; i++) {
  	if ( *cp < saveMax )
  		sp.push_back(*cp);

  	cp++;
  }

  cout << "Sorting " << sp.size() << " values..." << endl;
  std::sort(sp.begin(), sp.end());
  
  
  std::list<size_t> indices;
  std::list<double> pvals;
  std::list<double> adj;
  
  double currentMin = 1.0;

  cout << "Finding BH points..." << endl;
  for (size_t i = sp.size(); i > 0; i--) {
    double pval = sp[i-1];
    double adjVal = pval*n/i;
	  if ( adjVal < currentMin ) {
      indices.push_front(i);
      pvals.push_front(pval);
      adj.push_front(adjVal);
      currentMin = adjVal;
    }
  }
  
  Rcpp::DataFrame z = Rcpp::DataFrame::create(Rcpp::Named("rank")=Rcpp::wrap(indices),
                                              Rcpp::Named("pval")=Rcpp::wrap(pvals),
                                              Rcpp::Named("p.bh")=Rcpp::wrap(adj));
  return z;
}

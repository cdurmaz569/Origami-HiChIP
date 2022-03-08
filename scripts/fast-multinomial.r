library(Rcpp)
library(inline)

calc.zscore <- function(v) (v-mean(v))/sd(v)
bounded.prob <- function(v,low,high) pmin(pmax(v,low),high)
bounded.prob.matrix <- function(m,low,high) pmin(pmax(m/rowSums(m),low),high)


src <- '
int N = as<int>(n);
Rcpp::IntegerVector nSize(size);
Rcpp::NumericMatrix p(prob);
int ncol = p.ncol();
Rcpp::IntegerMatrix out(N,ncol);

if( N != nSize.size() || N != p.nrow() ) {
stop("currently, the size of each parameters must be the same");
}

RNGScope rngScope;

for( int i = 0; i < N; i++ ) {
Rcpp::NumericVector pv = p(i,_);
//std::sort(pv.begin(),pv.end());
Rcpp::NumericVector b(ncol);
double running = 0, total = 0;
double epsilon = .001;

for( int z = 0; z < ncol; z++ ) total += pv[z];

if( std::abs(total-1) > epsilon) {
for( int a = 0; a < ncol; a++ ) pv[a] /= total;
}

for( int x = 0; x < ncol; x++ ) {
running += pv[x];
b[x] = running;
}

int numElements = nSize[i];
Rcpp::NumericVector r = runif(numElements,0,1);
std::sort(r.begin(),r.end());
Rcpp::IntegerVector counts(ncol);

int idx = 0;
for( int j = 0; j < numElements; j++ ) {
while(r[j]>b[idx]) {
idx++;
if(idx>=ncol) break;
}

if(idx>=ncol) break; // sanity check

counts[idx]++;
}

out(i,_) = counts;
}

return out;
'

rmultinomial <- cxxfunction(signature(n = "integer", size="integer",prob="numeric"),
                            src, plugin="Rcpp")


\name{fclustIndex}
\title{Fuzzy Cluster Indexes (Validity/Performance Measures)}
\usage{fclustIndex(y, x, index = "all")}
\alias{fclustIndex}
\arguments{
  \item{y}{An object of a fuzzy clustering result of class \code{"fclust"}}
  \item{x}{Data matrix}
  \item{index}{The validity measures used: \code{"gath.geva"}, \code{"xie.beni"},
    \code{"fukuyama.sugeno"}, \code{"partition.coefficient"},
    \code{"partition.entropy"}, \code{"proportion.exponent"},
    \code{"separation.index"} and \code{"all"} for all the indexes.}}

    
\description{
  Calculates the values of several fuzzy validity measures. The values
  of the indexes can be independently used in order to evaluate and compare
  clustering partitions or even to determine the number of clusters
  existing in a data set.}

\details{
  The validity measures and a short description of them follows, where
  \eqn{N} is the number of data points, \eqn{u_{ij}} the values of the
  membership matrix, \eqn{v_j} the centers of the clusters and \eqn{k}
  te number of clusters.
  \describe{
    \item{\bold{gath.geva}:}{
      Gath and Geva introduced 2 main criteria for comparing and finding
      optimal partitions based on the heuristics that a better clustering
      assumes clear separation between the clusters, minimal volume of the
      clusters and maximal number of data points concentrated in the
      vicinity of the cluster centroids. These indexes are only for the
      cmeans clustering algorithm valid.
      For the first, the ``fuzzy hypervolume'' we have:
      \eqn{F_{HV}=\sum_{j=1}^{c}{[\det(F_j)]}^{1/2}}, where
      \eqn{F_j=\frac{\sum_{i=1}^N
	  u_{ij}(x_i-v_j)(x_i-v_j)^T}{\sum_{i=1}^{N}u_{ij}}}, for the
      case when the defuzzification parameter is 2.
      For the second, the ``average partition density'':
      \eqn{D_{PA}=\frac{1}{k}\sum_{j=1}^k\frac{S_j}{{[\det(F_j)]}^{1/2}}},
      where \eqn{S_j=\sum_{i=1}^N u_{ij}}.
      Moreover, the ``partition density'' which expresses the general
      partition density according to the physical definition of density
      is calculated by:
      \eqn{P_D=\frac{S}{F_{HV}}}, where \eqn{S=\sum_{j=1}^k\sum_{i=1}^N
	u_{ij}}.
    }
    \item{\bold{xie.beni}:}{
      This index is a function of the data set and the centroids of the
      clusters. Xie and Beni explained this index by writing it as a ratio
      of the total variation of the partition and the centroids $(U,V)$
      and the separation of the centroids vectors. The minimum values of
      this index under comparison support the best partitions.
      \eqn{u_{XB}(U,V;X)=\frac{\sum_{j=1}^k\sum_{i=1}^Nu_{ij}^2{||x_i-v_j||}^2}{N(\min_{j\neq l}\{{||v_j-v_l||}^2\})}}
    }
    \item{\bold{fukuyama.sugeno}:}{
      This index consists of the difference of two terms, the first
      combining the fuzziness in the membership matrix with the
      geometrical compactness of the representation of the data set via
      the prototypes, and the second the fuzziness in its row of the
      partition matrix with the distance from the $i$th prototype to the
      grand mean of the data. The minimum values of this index also
      propose a good partition.
      \eqn{u_{FS}(U,V;X)=\sum_{i=1}^{N}\sum_{j=1}^k
	(u_{ij}^2)^q(||x_i-v_j||^2-||v_j-\bar v||^2)}
    }
    \item{\bold{partition.coefficient}:}{
      An index which measures the fuzziness of the partition but without
      considering the data set itself. It is a heuristic measure since it
      has no connection to any property of the data. The maximum values of
      it imply a good partition in the meaning of a least fuzzy
      clustering.
      \eqn{F(U;k)=\frac{tr (UU^T)}{N}=\frac{<U,U>}{N}=\frac{||U||^2}{N}}
      \itemize{
	\item \eqn{F(U;k)} shows the fuzziness or the overlap of the partition
	and depends on \eqn{kN} elements. 
	\item \eqn{1/k\leq F(U;k)\leq 1}, where if \eqn{F(U;k)=1} then \eqn{U} is a hard
	partition and if \eqn{F(U;k)=1/k} then \eqn{U=[1/k]} is the centroid of
	the fuzzy partion space \eqn{P_{fk}}. The converse is also valid.
      }
    }
    \item{\bold{partition.entropy}:}{
      It is a measure that provides information about the membership
      matrix without also considering the data itself. The minimum values
      imply a good partition in the meaning of a more crisp partition.
      \eqn{H(U;k)=\sum_{i=1}^{N} h(u_i)/N}, where
      \eqn{h(u)=-\sum_{j=1}^{k} u_j\,\log _a (u_j)} the Shannon's entropy.
      \itemize{
	\item \eqn{H(U;k)} shows the uncertainty of a fuzzy partition and
	depends also on \eqn{kN} elements. Specifically, \eqn{h(u_i)} is
	interpreted as the amount of fuzzy information about the
	membership of \eqn{x_i} in \eqn{k} classes that is retained by column
	\eqn{u_j}. Thus, at \eqn{U=[1/k]} the most information is withheld since
	the membership is the fuzziest possible.
	\item \eqn{0\leq H(U;k)\leq \log_a(k)}, where for \eqn{H(U;k)=0} \eqn{U} is a
	hard partition and for \eqn{H(U;k)=\log_a(k)} \eqn{U=[1/k]}.
      }
    }
    \item{\bold{proportion.exponent}:}{
      It is a measure \eqn{P(U;k)} of fuzziness adept to detect structural variations
      in the partition matrix as it becomes more fuzzier. A crisp cluster
      in the partition matrix can drive it to infinity when the partition
      coefficient and the partition entropy are more sensitive to small
      changes when approaching a hard partition. Its evaluation does not also
      involve the data or the algorithm used to partition them and
      its maximum implies the optimal partition but without knowing what
      maximum is a statistically significant maximum.
      \itemize{
	\item \eqn{0\leq P(U;k)<\infty}, since the \eqn{[0,1]} values explode to
	\eqn{[0,\infty)} due to the natural logarithm. Specifically, \eqn{P=0}
	when and only when \eqn{U=[1/k]}, while \eqn{P\rightarrow\infty} when
	any column of \eqn{U} is crisp. 
	\item \eqn{P(U;k)} can easily explode and it is good for partitions
	with large column maximums and at detecting structural
	variations.
      }
    }
    \item{\bold{separation.index (known as CS Index)}:}{
      This index identifies unique cluster structure with well-defined
      properties that depend on the data and a measure of distance. It
      answers the question if the clusters are compact and separated, but
      it rather seems computationally infeasible for big data sets since a
      distance matrix between all the data membership values has to be
      calculated. It also presupposes that a hard partition is derived
      from the fuzzy one.\cr
      \eqn{D_1(U;k;X,d)=\min_{i+1\,\leq\,l\,\leq\,k-1}\left\{\min_{1\,\leq\,j\,\leq\,k}\left\{\frac{dis(u_j,u_l)}{\max_{1\leq m\leq k}\{dia(u_m)\}}\right\}\right\}}, where \eqn{dia}  is the diameter of the subset, \eqn{dis} the distance of
      two subsets, and \eqn{d} a metric.
      \eqn{U} is a CS partition of \eqn{X} \eqn{\Leftrightarrow D_1>1}. When this
      holds then \eqn{U} is unique.
    }
  }
}


\value{
  Returns a vector with the validity measures values.
}
\references{
  James C. Bezdek, \emph{Pattern Recognition with Fuzzy Objective
    Function Algorithms}, Plenum Press, 1981, NY.\cr
  L. X. Xie and G. Beni, \emph{Validity measure for fuzzy
    clustering}, IEEE Transactions on Pattern Analysis and Machine
  Intelligence, vol. \bold{3}, n. 8, p. 841-847, 1991.\cr
  I. Gath and A. B. Geva, \emph{Unsupervised Optimal Fuzzy
    Clustering}, IEEE Transactions on Pattern Analysis and Machine
  Intelligence, vol. \bold{11}, n. 7, p. 773-781, 1989.\cr
  Y. Fukuyama and M. Sugeno, \emph{A new method of choosing the
    number of clusters for the fuzzy $c$-means method}, Proc. 5th Fuzzy
  Syst. Symp., p. 247-250, 1989 (in japanese).}

\author{Evgenia Dimitriadou}

\seealso{\code{\link{cmeans}}}

\examples{
# a 2-dimensional example
x<-rbind(matrix(rnorm(100,sd=0.3),ncol=2),
         matrix(rnorm(100,mean=1,sd=0.3),ncol=2))
cl<-cmeans(x,2,20,verbose=TRUE,method="cmeans")
resultindexes <- fclustIndex(cl,x, index="all")
resultindexes   
}
\keyword{cluster}


\name{pathways}
\alias{pathways}
\docType{data}
\title{
Pathways List
}
\description{
Pathways IDs, Pathways Description, Number of Genes in the pathway 
}
\usage{data(pathways)}
\format{
Data file with pathway id, description, and number of genes in the pathway. PathId must correspond to the same ID of the pathways file names\cr
\preformatted{
      PathID     PathName                                   GenesCount
 [1,] "hsa00010" "Glycolysis / Gluconeogenesis"             "  66"    
 [2,] "hsa00020" "Citrate cycle (TCA cycle)"                "  31"    
 [3,] "hsa00030" "Pentose phosphate pathway"                "  27"    
 [4,] "hsa00040" "Pentose and glucuronate interconversions" "  32"    
 [5,] "hsa00051" "Fructose and mannose metabolism"          "  35"    
 [6,] "hsa00052" "Galactose metabolism"                     "  27"   
}
}
\source{
http://www.genome.jp/kegg/
}
\examples{
data(pathways)
}
\keyword{datasets}

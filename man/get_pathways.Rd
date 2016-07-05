\name{get_pathways}
\alias{get_pathways}
\title{
Pathways
}
\description{
Helper function to download pathways and their gene identifiers.
KEGG.db and reactome.db are used for pathway annotations.
}
\usage{
get_pathways(source="reactome",all_paths=TRUE,envir = "")
}
\arguments{
  \item{source}{
"reactome" or "kegg"
}
  \item{all_paths}{
TRUE or FALSE. If FALSE a subset will be asked by the function
}
  \item{envir}{
R environment to save Pathways to
}
}
\value{
Returns "Pathways description" All downloaded pathways are saved in the workspace
If reactome is selected as the source a prefix will be prompt to be set by user.
When kegg is selected the organism identifier is set automatically as the prefix (e.g."hsa").
}
\seealso{
\code{\link{read2_paths}}
}
\examples{
\dontrun{
# get pathways source = "kegg"
## library(KEGG.db)

# Create new environment to save data:
gper.env <- new.env()

# paths <- get_pathways(source="kegg",all_paths=FALSE,envir = gper.env)
	# when prompted introduce species as listed
#	hsa
	# if all_paths set to TRUE. All pathways are downloaded automatically
	# if all_paths set to FALSE. Introduce list of pathways separated by ","
	#hsa00010,hsa00020,hsa04670,hsa04672,hsa04710,hsa04720,hsa04722,hsa04730
	

# get pathways source = "reactome"
## library(reactome.db)
#paths <- get_pathways(source="reactome",all_paths=FALSE,envir=".GlobalEnv")
	# when prompted introduce species as listed
#	Homo sapiens
	# when prompted introduce prefix to be assigned to pathways
	#HSA
	# if all_paths set to TRUE. All pathways are downloaded automatically
	# IF all_paths set to FALSE, select a subset of pathway identifiers from 
	# list. Separated by ","
	1500931,1299503,...
}
}
\keyword{ annotation }
\keyword{ pathways }

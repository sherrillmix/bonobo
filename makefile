all: getData.md makeOtus.md

getData.md: getData.Rmd
	R -e 'knitr::knit("getData.Rmd")'

makeOtus.md: makeOtus.Rmd
	R -e 'knitr::knit("makeOtus.Rmd")'
	

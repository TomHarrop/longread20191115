cruft := slides.pdf slides.nav slides.log slides.aux slides.toc slides.snm
cache_directories := slides_cache slides_files

all: slides.pdf

clean:
ifneq ($(cruft),)
	rm -f $(cruft)
endif

clean_cache:
	$(RM) -rf $(cache_directories)

slides.pdf: slides.Rmd style/header.tex style/body.tex style/footer.tex
	# render the rmd to md using knitr
	R -e "rmarkdown::render('slides.Rmd')"

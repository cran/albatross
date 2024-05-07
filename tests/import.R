library(albatross)
library(tools)

ho_aq <- system.file('extdata/ho_aq.csv', package = 'albatross')
aq <- feem(ho_aq, 'table', sep = ',')

# 80% of values in the file are between -50 and 50, no NAs
stopifnot(
	sum(aq > -50 & aq < 50) / length(aq) > .80
)

pano1 <- as.data.frame(feem(
	system.file('extdata/panorama.dat', package = 'albatross'),
	'panorama'
))
with(pano1, stopifnot(
	# anti-Stokes + 10nm not measured
	emission >= excitation + 10,
	# intensity between .2 and .6
	intensity > .2, intensity < .6,
	min(intensity) < .25, max(intensity) > .55
))

# must understand all these parameters
if ('CP1251' %in% iconvlist()) {
	pano2 <- as.data.frame(feem(
		system.file('extdata/panorama.txt', package = 'albatross'),
		'table', fileEncoding = 'CP1251', transpose = TRUE, na = 0
	))
	with(pano2, stopifnot(
		# anti-Stokes not measured
		emission >= excitation,
		# 90% values between .25 and 3
		sum(intensity > .25 & intensity < 3) / length(intensity) > .9
	))
}

# pano2.txt should contain the same information as panorama.txt
pano3 <- as.data.frame(feem(
	system.file('extdata/pano2.txt', package = 'albatross'),
	'table', transpose = TRUE, na = 0
))
if ('CP1251' %in% iconvlist())
	stopifnot(all.equal(pano2, pano3))

# make sure that feem.connection works for both formats
f1 <- file(ho_aq)
feem(f1, 'table', sep = ',')

f2 <- file(system.file('extdata/panorama.dat', package = 'albatross'))
feem(f2, 'panorama')

# must stop on files in wrong format
assertError(feem(ho_aq, 'panorama'), verbose = TRUE)
assertError(feem(ho_aq, 'table'), verbose = TRUE)

# must handle empty Stokes area, too
pano4 <- as.data.frame(feem(
	system.file('extdata/pano-Stokes.dat', package = 'albatross'),
	'panorama'
))
with(pano4,
	stopifnot(range(emission - excitation) == c(10, 14))
)

# must handle flipped Panorama files
pano5 <- as.data.frame(feem(
	system.file('extdata/pano-flip.dat', package = 'albatross'),
	'panorama'
))
with(pano5, stopifnot(excitation >= emission + 10))

# must handle dec=','
dcs <- feem(
	system.file('extdata/decimal_comma.tsv', package = 'albatross'),
	'table', dec = ',', sep = '\t'
)
tools::assertError(feem( # must fail with wrong separator
	system.file('extdata/decimal_comma.tsv', package = 'albatross'),
	'table', dec = '.', sep = '\t'
), verbose = TRUE)
stopifnot(
	setequal(attr(dcs, 'excitation'), c(640, 645, 650)),
	dim(dcs) == c(5,3),
	!is.na(dcs),
	min(attr(dcs, 'emission')) < 245,
	max(attr(dcs, 'excitation')) > 247
)

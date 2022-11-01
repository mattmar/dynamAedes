w <- t(as.matrix(as.integer(seq(20,40,length.out=10))))*1000

test_that("To prove that seeded dynamAedes.m results are reproducible", {
	expect_equal(
		aeg.ws <- dynamAedes.m(species="aegypti", scale="ws", temps.matrix=w, startd="2022-06-01", n.clusters=1, iter=2, intro.eggs=100, jhwv=1, seeding=TRUE, compressed=TRUE, verbose=FALSE)
,
		aeg.ws1 <- dynamAedes.m(species="aegypti", scale="ws", temps.matrix=w, startd="2022-06-01", n.clusters=1, iter=2, intro.eggs=100, jhwv=1, seeding=TRUE, compressed=TRUE, verbose=FALSE)

		)
})
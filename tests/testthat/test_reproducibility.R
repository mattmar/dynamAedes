w <- t(as.matrix(as.integer(seq(20,40,length.out=10))))*1000

test_that("To prove that seeded dynamAedes results are reproducible", {
	expect_equal(
		aeg.ws <- dynamAedes(species="aegypti", scale="ws", temps.matrix=w, startd=2, endd=6, n.clusters=1, iter=2, intro.eggs=100, ihwv=1, seeding=TRUE, compressed=TRUE, verbose=FALSE)
,
		aeg.ws1 <- dynamAedes(species="aegypti", scale="ws", temps.matrix=w, startd=2, endd=6, n.clusters=1, iter=2, intro.eggs=100, ihwv=1, seeding=TRUE, compressed=TRUE, verbose=FALSE)

		)
})
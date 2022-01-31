context('Testing options')

pops <- options('psd.ops')[[1]]

test_that('psd options load',{
  expect_is(pops, 'list')
  expect_length(pops, 3)
  expect_named(pops, c('tapmin','tapcap','names'))
})

test_that('psd taper options are set appropriately',{
  expect_gte(pops[['tapmin']], 1)
  expect_gte(pops[['tapcap']], 1)
})

test_that('psd env names are set appropriately',{
  pnms <- pops[['names']]
  expect_is(pnms, 'list')
  expect_length(pnms, 11)
})
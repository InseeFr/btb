library('RUnit')

testsuite.btb <- defineTestSuite("btb",
                                 dirs = paste0(getwd(), "/btb/inst/testsUnitaires"),
                                 testFileRegexp = "^test_.+\\.R",
                                 testFuncRegexp = "^test.+",
                                 rngKind = "Marsaglia-Multicarry",
                                 rngNormalKind = "Kinderman-Ramage")

testResult <- runTestSuite(testsuite.btb)
printTextProtocol(testResult)

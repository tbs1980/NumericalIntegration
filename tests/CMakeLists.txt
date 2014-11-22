ADD_EXECUTABLE(outputKronrodNodesWeights outputKronrodNodesWeights.cpp)
#ADD_EXECUTABLE(IntegratorTest IntegratorTest.cpp)
ADD_EXECUTABLE(testLaurieGautschi testLaurieGautschi.cpp)
ADD_EXECUTABLE(testIntegrator_LogPow testIntegrator_LogPow.cpp)
ADD_EXECUTABLE(testIntegrator_Sine testIntegrator_Sine.cpp)
ADD_EXECUTABLE(testIntegrator_Peak testIntegrator_Peak.cpp)

TARGET_LINK_LIBRARIES(outputKronrodNodesWeights ${ALL_LIBS})
#TARGET_LINK_LIBRARIES(IntegratorTest ${ALL_LIBS})
TARGET_LINK_LIBRARIES(testLaurieGautschi ${ALL_LIBS})
TARGET_LINK_LIBRARIES(testIntegrator_LogPow ${ALL_LIBS})
TARGET_LINK_LIBRARIES(testIntegrator_Sine ${ALL_LIBS})
TARGET_LINK_LIBRARIES(testIntegrator_Peak ${ALL_LIBS})



#ADD_EXECUTABLE(test
#  testIntegrator_LogPow.cpp
#  testIntegrator_Sine.cpp
#  testIntegrator_Peak.cpp
#  ../source/
#  )

#TARGET_LINK_LIBRARIES(test
#  GTestMain
#  )

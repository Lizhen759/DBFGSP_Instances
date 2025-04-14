cmake_minimum_required(VERSION 3.22)
project(DBFGSP_new)

set(CMAKE_CXX_STANDARD 17)

add_executable(DBFGSP_new main.cpp Problem.cpp Problem.h NOperator.cpp NOperator.h Heu.cpp Heu.h Base.cpp Base.h rand.cpp rand.h TSCCEA.cpp TSCCEA.h Individual.cpp Individual.h ParetoReader.h ParetoReader.cpp ParetoOptimizer.cpp ParetoOptimizer.h GCCEA.cpp GCCEA.h MOSA.cpp MOSA.h CMOEA.cpp CMOEA.h
        QCCEA.cpp QCCEA.h QCCEA.h QCCEA.cpp NOperator_New.cpp NOperator_New.h QCCEANoQ.cpp QCCEANoQ.h QCCEANoRapid.h QCCEANoRapid.cpp)
